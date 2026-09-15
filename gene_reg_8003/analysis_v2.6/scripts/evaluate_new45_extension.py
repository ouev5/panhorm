#!/usr/bin/env python3
"""Score the fresh 45-pair extension with the frozen v2.6 deterministic pipeline.

Identical to validation/scripts/evaluate_genereg_orthogonal_v10.py except that
the expected record count is read from the input file instead of being fixed at
370, so it can drive the 45-record extension set.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import random
import time
from datetime import datetime, timezone
from pathlib import Path

THRESHOLD = 37.0
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20260716


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metric_values(labels, predictions):
    tp = sum(y == 1 and p == 1 for y, p in zip(labels, predictions))
    tn = sum(y == 0 and p == 0 for y, p in zip(labels, predictions))
    fp = sum(y == 0 and p == 1 for y, p in zip(labels, predictions))
    fn = sum(y == 1 and p == 0 for y, p in zip(labels, predictions))
    sensitivity = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision + sensitivity else 0.0
    accuracy = (tp + tn) / len(labels) if labels else 0.0
    return {"TP": tp, "TN": tn, "FP": fp, "FN": fn, "sensitivity": sensitivity,
            "specificity": specificity, "precision": precision, "f1": f1, "accuracy": accuracy}


def bootstrap(labels, predictions):
    rng = random.Random(BOOTSTRAP_SEED)
    values = {name: [] for name in ("sensitivity", "specificity", "precision", "f1", "accuracy")}
    n = len(labels)
    for _ in range(BOOTSTRAP_REPLICATES):
        indices = [rng.randrange(n) for _ in range(n)]
        sample = metric_values([labels[i] for i in indices], [predictions[i] for i in indices])
        for name in values:
            values[name].append(float(sample[name]))
    return {name: [round(sorted(v)[24], 6), round(sorted(v)[974], 6)] for name, v in values.items()}


def load_app(path: Path):
    spec = importlib.util.spec_from_file_location("genereg_app_for_extension_eval", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import GeneReg application from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-set", required=True, type=Path)
    parser.add_argument("--app", default="/www/wwwroot/gene_reg/app.py", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.test_set.open(encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    total = len(records)

    app = load_app(args.app)
    rows = []
    for index, record in enumerate(records, 1):
        started = time.monotonic()
        row = {
            "record_id": record["record_id"], "label": int(record["label"]), "subset": record["subset"],
            "tf": record["tf"], "gene": record["gene"], "threshold": THRESHOLD,
        }
        try:
            original_fetch = app.gf.fetch
            expected_ensembl_id = (record.get("ensembl_id") or "").strip()

            def fetch_with_frozen_ensembl(symbol):
                info = original_fetch(symbol)
                if info and expected_ensembl_id and symbol.strip().upper() == record["gene"].strip().upper():
                    info = dict(info)
                    info["ensembl_id"] = expected_ensembl_id
                return info

            app.gf.fetch = fetch_with_frozen_ensembl
            try:
                result = app._do_analysis(record["gene"], record["tf"], include_ai=False)
            finally:
                app.gf.fetch = original_fetch
            raw_score = float(result["raw_score"])
            scores = result.get("scores", {})
            row.update({
                "raw_score": raw_score,
                "prediction": int(raw_score >= THRESHOLD),
                "classification": result.get("classification", ""),
                "chip_score": scores.get("chip_seq", {}).get("score", ""),
                "motif_score": scores.get("motif", {}).get("score", ""),
                "literature_score": scores.get("literature", {}).get("score", ""),
                "ai_bonus": result.get("ai_bonus", ""),
                "status": "ok", "error": "",
            })
        except Exception as exc:  # noqa: BLE001
            row.update({"raw_score": "", "prediction": "", "classification": "", "chip_score": "",
                        "motif_score": "", "literature_score": "", "ai_bonus": "",
                        "status": "error", "error": f"{type(exc).__name__}: {exc}"})
        row["runtime_seconds"] = round(time.monotonic() - started, 3)
        rows.append(row)
        (args.output_dir / "progress.json").write_text(json.dumps({
            "completed": index, "total": total, "last_record_id": record["record_id"],
            "updated_at_utc": datetime.now(timezone.utc).isoformat()}, indent=2) + "\n", encoding="utf-8")
        print(f"{index}/{total} {record['tf']}->{record['gene']} {row['status']} score={row['raw_score']}", flush=True)

    fields = ["record_id", "label", "subset", "tf", "gene", "threshold", "raw_score", "prediction",
              "classification", "chip_score", "motif_score", "literature_score", "ai_bonus", "status", "error", "runtime_seconds"]
    with (args.output_dir / "genereg_deterministic_predictions_45.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    successful = [row for row in rows if row["status"] == "ok"]
    if len(successful) != total:
        raise RuntimeError(f"{total - len(successful)} analyses failed; row-level audit was retained")
    labels = [int(row["label"]) for row in successful]
    predictions = [int(row["prediction"]) for row in successful]
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "test_set": str(args.test_set), "test_set_sha256": sha256(args.test_set),
        "app_path": str(args.app), "app_sha256": sha256(args.app),
        "n_total": total, "n_successful": len(successful),
        "positive_count": sum(labels), "negative_count": total - sum(labels),
        "decision_rule": "raw_score >= 37.0 (v2.6 frozen rule, reported for reference only)",
        "metrics": metric_values(labels, predictions),
        "bootstrap_95ci": bootstrap(labels, predictions),
    }
    (args.output_dir / "genereg_deterministic_metrics.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary["metrics"], indent=2), flush=True)


if __name__ == "__main__":
    main()
