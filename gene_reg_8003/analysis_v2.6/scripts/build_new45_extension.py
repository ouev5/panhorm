#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a fresh, non-overlapping 45-pair calibration extension (22 positive + 23 negative).

Rationale: the retired 45-pair calibration set (Table S2) shares record IDs with
the frozen 370-pair benchmark.  For the leakage-free 415-pair re-analysis we
replace those records with NEW pairs drawn under the SAME rules as the frozen
release, so that 370 + 45 = 415 unique TF-gene pairs:

- Positives: the next-ranked GTRD promoter targets (same eligibility filter and
  sort as ``select_gtrd_positives`` in the frozen build script), taken after the
  frozen top-20 per TF, round-robin across the eight TFs.
- Negatives: the same strict four-criterion screen (absent from the TF-specific
  GTRD list, ChIP-Atlas average = 0, no JASPAR PWM hit at relative threshold
  0.95 in the [-2000,+500] promoter, zero PubMed Title/Abstract co-occurrence),
  reshuffled with a new deterministic seed and excluding every pair already in
  the frozen 370 set.

The frozen release directory itself is never modified (caches are copied).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import random
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

FROZEN = Path("/www/wwwroot/gene_reg/validation/expanded_validation_v1_20260715_r7_orthogonal")
BUILD_SCRIPT = Path("/www/wwwroot/gene_reg/validation/scripts/build_expanded_genereg_validation.py")
OUT = Path("/www/wwwroot/gene_reg/validation/calibration_extension_45_20260910")
NEW_RUN_ID = "calib45ext_v1_20260910"
N_POS = 22
N_NEG = 23

spec = importlib.util.spec_from_file_location("frozen_build", BUILD_SCRIPT)
fb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fb)

TF_ORDER = list(fb.TF_CONFIG)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_final_370() -> tuple[list[dict], set[tuple[str, str]]]:
    path = FROZEN / "derived" / "combined_test_set_178_positive_192_negative.csv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 370, f"expected 370 frozen rows, got {len(rows)}"
    return rows, {(r["tf"], r["gene"]) for r in rows}


def load_sources() -> tuple[dict[str, dict], dict[str, dict]]:
    gtrd_by_tf: dict[str, dict] = {}
    chip_by_tf: dict[str, dict] = {}
    for tf, config in fb.TF_CONFIG.items():
        gtrd_by_tf[tf] = fb.load_gtrd_targets(FROZEN / "raw/gtrd/target_genes" / f"{tf}_{config['uniprot']}.tsv")
        chip_by_tf[tf] = fb.load_chip_atlas_targets(
            FROZEN / "raw/chip_atlas" / f"{tf}_hg38_tss1kb.tsv",
            fb.local_chip_atlas_experiment_ids(tf),
        )
    return gtrd_by_tf, chip_by_tf


def eligible_positive_ranking(tf, gtrd_by_tf, chip_by_tf, excluded_pairs):
    """Same filter and sort as frozen select_gtrd_positives, minus used pairs."""
    candidates = []
    for gene, source in gtrd_by_tf[tf].items():
        chip = chip_by_tf[tf].get(gene, {
            "chip_atlas_average": 0.0,
            "chip_atlas_nonzero_experiments": 0,
            "chip_atlas_independent_nonzero_experiments": 0,
            "chip_atlas_independent_experiment_universe": 0,
        })
        if (
            (tf, gene) in excluded_pairs
            or not fb.valid_gene_symbol(gene)
            or source["gtrd_site_count"] < 2
            or chip["chip_atlas_average"] <= 0
            or chip["chip_atlas_independent_nonzero_experiments"] < 2
        ):
            continue
        candidates.append({"gene": gene, **source, **chip})
    candidates.sort(key=lambda item: (-item["gtrd_site_count"], -item["chip_atlas_average"], item["gene"]))
    return candidates


def pick_positives(gtrd_by_tf, chip_by_tf, used_pairs) -> tuple[list[dict], list[dict]]:
    rankings = {tf: eligible_positive_ranking(tf, gtrd_by_tf, chip_by_tf, used_pairs) for tf in TF_ORDER}
    audit: list[dict] = []
    for tf, ranking in rankings.items():
        for rank, item in enumerate(ranking[:40], start=1):
            audit.append({
                "tf": tf, "rank_after_frozen_top20": rank, "gene": item["gene"],
                "ensembl_id": item["ensembl_id"], "gtrd_site_count": item["gtrd_site_count"],
                "chip_atlas_average": item["chip_atlas_average"],
                "chip_atlas_independent_nonzero_experiments": item["chip_atlas_independent_nonzero_experiments"],
                "picked": False,
            })
    quotas = {}
    base, extra = divmod(N_POS, len(TF_ORDER))
    for i, tf in enumerate(TF_ORDER):
        quotas[tf] = base + (1 if i < extra else 0)  # first TFs take the remainder
    picked: list[dict] = []
    cursor = {tf: 0 for tf in TF_ORDER}
    for tf in TF_ORDER:
        while len([p for p in picked if p["tf"] == tf]) < quotas[tf]:
            item = rankings[tf][cursor[tf]]
            cursor[tf] += 1
            picked.append({
                "label": 1,
                "subset": "calib45ext_GTRD_positive",
                "tf": tf,
                "gene": item["gene"],
                "ensembl_id": item["ensembl_id"],
                "gtrd_site_count": item["gtrd_site_count"],
                "chip_atlas_average": item["chip_atlas_average"],
                "chip_atlas_nonzero_experiments": item["chip_atlas_nonzero_experiments"],
                "chip_atlas_independent_nonzero_experiments": item["chip_atlas_independent_nonzero_experiments"],
                "chip_atlas_independent_experiment_universe": item["chip_atlas_independent_experiment_universe"],
                "evidence_summary": "Extension positive: GTRD promoter target with >=2 meta-clusters plus >=2 non-zero ChIP-Atlas SRX experiments absent from the local GeneReg peak corpus; next-ranked after the frozen top-20 for its TF, disjoint from the frozen 370-pair set.",
                "selection_rule": "gtrd_site_count>=2 AND chip_atlas_average>0 AND chip_atlas_independent_nonzero_experiments>=2 AND pair_not_in_frozen370; rank order inherited from the frozen sort",
            })
    for row in audit:
        row["picked"] = any(p["tf"] == row["tf"] and p["gene"] == row["gene"] for p in picked)
    assert len(picked) == N_POS, f"expected {N_POS} positives, got {len(picked)}"
    return picked, audit


def screen_negatives(session, gtrd_by_tf, chip_by_tf, used_pairs, new_positives, promoter_cache, jaspar_cache):
    all_genes: dict[str, str] = {}
    for targets in gtrd_by_tf.values():
        for gene, source in targets.items():
            if fb.valid_gene_symbol(gene):
                all_genes.setdefault(gene, source["ensembl_id"])
    manual_pairs = {(item["tf"], item["gene"]) for item in fb.RECENT_LITERATURE_PAIRS}
    new_positive_pairs = {(item["tf"], item["gene"]) for item in new_positives}

    quotas = {}
    base, extra = divmod(N_NEG, len(TF_ORDER))
    for i, tf in enumerate(TF_ORDER):
        quotas[tf] = base + (1 if i < extra else 0)

    negatives: list[dict] = []
    audit: list[dict] = []
    for tf in TF_ORDER:
        matrix = fb.fetch_jaspar_matrix(session, tf, fb.TF_CONFIG[tf]["jaspar_matrix_id"], jaspar_cache)
        candidates = [
            (gene, ensembl_id)
            for gene, ensembl_id in all_genes.items()
            if gene not in gtrd_by_tf[tf]
            and chip_by_tf[tf].get(gene, {}).get("chip_atlas_average", 0.0) <= 0
            and fb.valid_gene_symbol(gene)
            and (tf, gene) not in manual_pairs
            and (tf, gene) not in used_pairs
            and (tf, gene) not in new_positive_pairs
        ]
        rng = random.Random(f"{NEW_RUN_ID}:{tf}")
        rng.shuffle(candidates)
        picked = 0
        index = 0
        while picked < quotas[tf] and index < len(candidates):
            batch = candidates[index:index + 16]
            index += len(batch)
            sequence_by_id, errors_by_id = fb.fetch_promoters_bulk(session, [eid for _, eid in batch], promoter_cache)
            for gene, ensembl_id in sorted(batch):
                if picked >= quotas[tf]:
                    break
                row = {
                    "tf": tf, "gene": gene, "ensembl_id": ensembl_id,
                    "gtrd_absent": True, "chip_atlas_average": 0.0,
                    "chip_atlas_nonzero_experiments": 0,
                    "chip_atlas_independent_nonzero_experiments": 0,
                    "chip_atlas_independent_experiment_universe": matrix.get("chip_atlas_independent_experiment_universe", ""),
                    "jaspar_matrix_id": matrix["matrix_id"],
                    "jaspar_source_url": matrix["source_url"],
                    "jaspar_relative_threshold": fb.JASPAR_RELATIVE_THRESHOLD,
                    "pubmed_checked_at_utc": utc_now(),
                }
                if ensembl_id in errors_by_id or ensembl_id not in sequence_by_id:
                    row["screen_status"] = "rejected_source_error"
                    row["error"] = errors_by_id.get(ensembl_id, "No promoter sequence returned")
                    audit.append(row)
                    continue
                hits, raw_threshold = fb.pwm_hits(sequence_by_id[ensembl_id], matrix)
                row["jaspar_pwm_raw_threshold"] = raw_threshold
                row["jaspar_pwm_hits"] = hits
                if hits != 0:
                    row["screen_status"] = "rejected_jaspar_hit"
                    audit.append(row)
                    continue
                try:
                    count, query = fb.pubmed_count(session, tf, gene)
                except Exception as exc:  # noqa: BLE001
                    row["screen_status"] = "rejected_pubmed_error"
                    row["error"] = f"{type(exc).__name__}: {exc}"
                    audit.append(row)
                    continue
                row["pubmed_query"] = query
                row["pubmed_count"] = count
                if count != 0:
                    row["screen_status"] = "rejected_pubmed_cooccurrence"
                    audit.append(row)
                else:
                    row["screen_status"] = "accepted_negative"
                    audit.append(row)
                    negatives.append({
                        "label": 0,
                        "subset": "calib45ext_strict_negative",
                        "tf": tf,
                        "gene": gene,
                        "ensembl_id": ensembl_id,
                        "gtrd_site_count": 0,
                        "chip_atlas_average": 0.0,
                        "chip_atlas_nonzero_experiments": 0,
                        "chip_atlas_independent_nonzero_experiments": 0,
                        "chip_atlas_independent_experiment_universe": 0,
                        "jaspar_matrix_id": matrix["matrix_id"],
                        "jaspar_pwm_hits": 0,
                        "jaspar_relative_threshold": fb.JASPAR_RELATIVE_THRESHOLD,
                        "pubmed_query": query,
                        "pubmed_count": 0,
                        "evidence_summary": "Extension strict negative: absent from the TF-specific GTRD promoter target list, ChIP-Atlas score=0, no high-stringency JASPAR PWM hit in the promoter window, and zero PubMed title/abstract co-occurrence; disjoint from the frozen 370-pair set.",
                        "selection_rule": "gtrd_absent AND chip_atlas_average=0 AND jaspar_pwm_hits=0 AND pubmed_count=0 AND pair_not_in_frozen370",
                    })
                    picked += 1
                time.sleep(0.36)
            print(f"{tf}: extension negatives accepted {picked}/{quotas[tf]}; screened {index} candidates", flush=True)
        if picked != quotas[tf]:
            raise RuntimeError(f"{tf}: only {picked}/{quotas[tf]} extension negatives found")
    return negatives, audit


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"output already exists: {OUT}")
    OUT.mkdir(parents=True)
    (OUT / "raw").mkdir()

    final_rows, used_pairs = load_final_370()
    gtrd_by_tf, chip_by_tf = load_sources()

    positives, positive_audit = pick_positives(gtrd_by_tf, chip_by_tf, used_pairs)
    print(f"positives picked: {len(positives)}", flush=True)

    promoter_cache = OUT / "raw/ensembl_promoters"
    jaspar_cache = OUT / "raw/jaspar"
    shutil.copytree(FROZEN / "raw/ensembl_promoters", promoter_cache, copy_function=shutil.copy2)
    shutil.copytree(FROZEN / "raw/jaspar", jaspar_cache, copy_function=shutil.copy2)

    session = requests.Session()
    session.headers.update({"Accept": "application/json, text/plain, */*"})
    negatives, negative_audit = screen_negatives(session, gtrd_by_tf, chip_by_tf, used_pairs, positives, promoter_cache, jaspar_cache)
    print(f"negatives picked: {len(negatives)}", flush=True)

    combined = positives + negatives
    for i, row in enumerate(combined, 1):
        row["record_id"] = f"{NEW_RUN_ID}_{i:04d}"
        row["species"] = fb.SPECIES
        row["used_in_training"] = "false"
        row["pmid"] = ""
        row["publication_year"] = ""

    pair_keys = [(r["tf"], r["gene"]) for r in combined]
    assert len(pair_keys) == len(set(pair_keys)), "duplicate TF-gene pair inside the extension"
    assert not (set(pair_keys) & used_pairs), "extension overlaps the frozen 370 set"
    assert sum(int(r["label"]) for r in combined) == N_POS
    assert N_NEG == sum(1 for r in combined if int(r["label"]) == 0)

    fields = [
        "record_id", "label", "subset", "species", "tf", "gene", "ensembl_id", "used_in_training",
        "gtrd_site_count", "chip_atlas_average", "chip_atlas_nonzero_experiments", "jaspar_matrix_id",
        "chip_atlas_independent_nonzero_experiments", "chip_atlas_independent_experiment_universe",
        "jaspar_pwm_hits", "jaspar_relative_threshold", "pubmed_query", "pubmed_count", "pmid", "publication_year",
        "evidence_summary", "selection_rule",
    ]
    write_csv(OUT / "new45_extension_test_set.csv", combined, fields)
    write_csv(OUT / "positive_candidate_ranking_audit.csv", positive_audit,
              ["tf", "rank_after_frozen_top20", "gene", "ensembl_id", "gtrd_site_count",
               "chip_atlas_average", "chip_atlas_independent_nonzero_experiments", "picked"])
    write_csv(OUT / "negative_screening_audit_new45.csv", negative_audit,
              ["tf", "gene", "ensembl_id", "gtrd_absent", "chip_atlas_average", "chip_atlas_nonzero_experiments",
               "chip_atlas_independent_nonzero_experiments", "chip_atlas_independent_experiment_universe",
               "jaspar_matrix_id", "jaspar_source_url", "jaspar_relative_threshold", "jaspar_pwm_raw_threshold",
               "jaspar_pwm_hits", "pubmed_query", "pubmed_count", "pubmed_checked_at_utc", "screen_status", "error"])

    summary = {
        "created_at_utc": utc_now(),
        "new_run_id": NEW_RUN_ID,
        "purpose": "45 fresh pairs replacing the retired Table S2 calibration set for the leakage-free 415-pair re-analysis",
        "frozen_reference": str(FROZEN),
        "counts": {"total": len(combined), "positive": N_POS, "negative": N_NEG},
        "positive_rule": "next-ranked GTRD promoter targets after the frozen top-20 per TF (same eligibility filter and sort as the frozen build), round-robin quotas across the eight TFs",
        "negative_rule": "strict four-criterion screen identical to the frozen protocol, candidate reshuffle seed f'{NEW_RUN_ID}:<tf>', excluding every pair already inside the frozen 370 set",
        "positive_quotas": {tf: (N_POS // 8 + (1 if i < N_POS % 8 else 0)) for i, tf in enumerate(TF_ORDER)},
        "negative_quotas": {tf: (N_NEG // 8 + (1 if i < N_NEG % 8 else 0)) for i, tf in enumerate(TF_ORDER)},
        "overlap_with_frozen370": 0,
        "records": [{"record_id": r["record_id"], "label": int(r["label"]), "subset": r["subset"], "tf": r["tf"], "gene": r["gene"]} for r in combined],
    }
    (OUT / "build_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("new_run_id", "counts", "overlap_with_frozen370")}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
