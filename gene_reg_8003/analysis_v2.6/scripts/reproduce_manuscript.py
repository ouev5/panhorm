"""Replay the published split and selected model without network or LLM calls.

By default, verifies all published row-level predictions and reports metrics.
Use --cross-validate to replay weight selection on development rows only.
Use --write-summary to refresh the machine-readable validation summary.
"""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scoring import MODEL_CONFIG, score_evidence
from rerun_415_leakfree import (
    BOOTSTRAP_SEED, bootstrap_metrics_ci, confusion, cv_evaluate_weights,
    select_best_weights,
)

CORPUS = ROOT / "analysis_v2.6/manuscript_benchmark_415.csv"
SUMMARY = ROOT / "analysis_v2.6/manuscript_metrics.json"


def read_corpus():
    with CORPUS.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["label"] = int(row["label"])
        for name in ("chip_score", "motif_score", "literature_score"):
            row[name] = float(row[name])
    return rows


def rank_metrics(labels, scores):
    """ROC AUC and trapezoidal PR AUC, grouped by unique score thresholds."""
    labels = np.asarray(labels)
    scores = np.asarray(scores)
    positives, negatives = int(labels.sum()), int((labels == 0).sum())
    tpr, fpr, recall, precision = [0.0], [0.0], [0.0], [1.0]
    for threshold in sorted(set(scores), reverse=True):
        prediction = scores >= threshold
        tp = int(np.sum(prediction & (labels == 1)))
        fp = int(np.sum(prediction & (labels == 0)))
        tpr.append(tp / positives)
        fpr.append(fp / negatives)
        recall.append(tp / positives)
        precision.append(tp / (tp + fp))
    return {"roc_auc": float(np.trapezoid(tpr, fpr)), "pr_auc": float(np.trapezoid(precision, recall))}


def build_summary(rows):
    assert len(rows) == len({row["record_id"] for row in rows}) == len({(row["tf"], row["gene"]) for row in rows}) == 415
    dev = [row for row in rows if row["partition"] == "development"]
    test = [row for row in rows if row["partition"] == "test"]
    assert (len(dev), sum(row["label"] for row in dev)) == (332, 160)
    assert (len(test), sum(row["label"] for row in test)) == (83, 40)
    assert not {(row["tf"], row["gene"]) for row in dev} & {(row["tf"], row["gene"]) for row in test}
    dev_scores = [score_evidence(row["chip_score"], row["motif_score"], row["literature_score"])["weighted_score"] for row in dev]
    assert min(dev_scores) == MODEL_CONFIG["normalization"]["min"]
    assert max(dev_scores) == MODEL_CONFIG["normalization"]["max"]
    for row in rows:
        result = score_evidence(row["chip_score"], row["motif_score"], row["literature_score"])
        assert abs(result["weighted_score"] - float(row["weighted_score"])) < 1e-6
        assert abs(result["normalized_score"] - float(row["weighted_score_normalized"])) < 1e-6
        assert result["prediction"] == int(row["prediction_final"])
    labels = [row["label"] for row in test]
    test_results = [score_evidence(row["chip_score"], row["motif_score"], row["literature_score"]) for row in test]
    scores = [result["normalized_score"] for result in test_results]
    predictions = [result["prediction"] for result in test_results]
    metrics = confusion(labels, predictions)
    metrics["recall"] = metrics.pop("recall_sensitivity")
    metrics["tn_proportion_evidence_absent"] = metrics.pop("specificity")
    metrics.update(rank_metrics(labels, scores))
    assert [metrics[name] for name in ("TP", "TN", "FP", "FN")] == [38, 36, 7, 2]
    ci = bootstrap_metrics_ci(labels, predictions, BOOTSTRAP_SEED)
    ci["recall"] = ci.pop("recall_sensitivity")
    ci["tn_proportion_evidence_absent"] = ci.pop("specificity")
    return {
        "model": MODEL_CONFIG,
        "scope": "Internal benchmark; evidence-absent labels are not experimentally validated biological negatives.",
        "decision_rule": "(0.10 * ChIP-seq + 0.70 * motif + 0.20 * literature - 24.50) / (67.00 - 24.50) >= 0.41",
        "corpus": {"path": "analysis_v2.6/manuscript_benchmark_415.csv", "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest(), "records": 415, "positives": 200, "evidence_absent": 215},
        "development_set": {"records": 332, "positives": 160, "evidence_absent": 172},
        "test_set": {"records": 83, "positives": 40, "evidence_absent": 43},
        "weight_selection": {"candidate_combinations": 231, "eligible_combinations": 171, "folds": 5, "one_se_cutoff": 0.9047631392051609, "selected_mean_cv_f1": 0.9083148708088927, "tie_break": "Equal development peak F1; higher mean F1 over the development threshold grid (0.6004 versus 0.5704)."},
        "metrics": metrics,
        "bootstrap": {"replicates": 1000, "seed": BOOTSTRAP_SEED, "ci_method": "percentile 2.5th-97.5th"},
        "bootstrap_95ci": ci,
        "ai_score_contribution": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cross-validate", action="store_true")
    parser.add_argument("--write-summary", action="store_true")
    args = parser.parse_args()
    rows = read_corpus()
    if args.cross_validate:
        dev = [row for row in rows if row["partition"] == "development"]
        best, _, eligible, diagnostics = select_best_weights(cv_evaluate_weights(dev), development_rows=dev)
        assert best["weights_chip_motif_literature"] == [0.10, 0.70, 0.20]
        print(f"Development-only CV verified: {len(eligible)} eligible combinations; 1-SE cutoff {diagnostics['one_se_cutoff']:.4f}")
    summary = build_summary(rows)
    if args.write_summary:
        SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
