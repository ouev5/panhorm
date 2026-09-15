# -*- coding: utf-8 -*-
"""GeneReg 415-pair leakage-free re-analysis (manuscript Response 2 workflow).

Merges the frozen 370-pair benchmark with the fresh 45-pair extension
(calib45ext_v1_20260910, disjoint from the 370) into a unique 415-pair corpus,
then:

1. Stratified 80/20 split with a fixed seed -> 332-pair development set /
   83-pair independent test set (the test set is touched exactly once).
2. Stratified 5-fold cross-validation on the development set over the legacy
   three-evidence weight grid (0.05 steps, 231 combinations).  Inside each fold
   the classification threshold is selected only on the training folds (train
   min-max normalisation + 0.00-1.00 grid, F1 maximised), then applied
   unchanged to the validation fold.
3. Mean validation F1 per weight combination with a 1000-resample bootstrap
   95% CI (resampling the five fold-level F1 values).
4. Best weight combination chosen on the development CV SUBJECT TO the
   positivity constraint that every evidence weight is at least MIN_WEIGHT
   (one grid step, 0.05), i.e. all three evidence streams (ChIP-seq, motif,
   literature) must contribute to the integrated score; zero-weight
   combinations are reported but not eligible.  Final threshold from the full
   development set only (min-max normalisation + 0.00-1.00 grid).
5. One-shot evaluation on the 83-pair test set: F1, precision, recall
   (sensitivity), specificity, accuracy, with 1000-resample percentile
   bootstrap 95% CIs.  Normalisation parameters and threshold are frozen from
   the development set before the test set is scored.

Deterministic tie-breaks (documented): weight selection ties are broken by
higher mean fold specificity, then higher mean fold sensitivity, then smaller
distance to the published 0.50/0.25/0.25 (ChIP/motif/literature) combination,
then lexicographic weight order.  Threshold ties are broken by higher
specificity, then higher sensitivity, then the threshold closer to 0.37.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(r"E:\9.2_manuscript")
OUT = ROOT / "GeneReg_415独立重跑_1SE均衡权重_20260910"
SERVER_PRED_370 = ROOT / "审稿意见修订_20260902" / "server_v2.6_release" / "validation" / "expanded_validation_v1_20260715_r7_orthogonal" / "evaluation_deterministic_v11_20260902_v2.6_frozen" / "genereg_deterministic_predictions.csv"
SERVER_BENCH_370 = ROOT / "审稿意见修订_20260902" / "server_v2.6_release" / "validation" / "expanded_validation_v1_20260715_r7_orthogonal" / "derived" / "combined_test_set_178_positive_192_negative.csv"
NEW45_SET = Path(r"C:\Users\huang'xing'xing\Desktop\小组项目\server_data\new45_extension_test_set.csv")
NEW45_PRED = Path(r"C:\Users\huang'xing'xing\Desktop\小组项目\server_data\genereg_deterministic_predictions_45.csv")

SEED = 20260910
SPLIT_SEED = 20260910
CV_SEED = 20260910
BOOTSTRAP_SEED = 20260910
BOOTSTRAP_N = 1000
DEV_N = 332
TEST_N = 83
N_FOLDS = 5
PUBLISHED_WEIGHTS = (0.50, 0.25, 0.25)  # ChIP-seq, motif, literature
MIN_WEIGHT = 0.05  # every evidence stream must contribute: no zero weights


def read_csv(path: Path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def confusion(labels, pred):
    labels = np.asarray(labels, dtype=int)
    pred = np.asarray(pred, dtype=int)
    tp = int(np.sum((labels == 1) & (pred == 1)))
    tn = int(np.sum((labels == 0) & (pred == 0)))
    fp = int(np.sum((labels == 0) & (pred == 1)))
    fn = int(np.sum((labels == 1) & (pred == 0)))
    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    prec = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * prec * sens / (prec + sens) if prec + sens else 0.0
    acc = (tp + tn) / len(labels) if len(labels) else 0.0
    return {"TP": tp, "TN": tn, "FP": fp, "FN": fn, "recall_sensitivity": sens,
            "specificity": spec, "precision": prec, "f1": f1, "accuracy": acc}


def merge_415():
    pred370 = {r["record_id"]: r for r in read_csv(SERVER_PRED_370)}
    bench370 = {r["record_id"]: r for r in read_csv(SERVER_BENCH_370)}
    assert len(bench370) == 370 and len(pred370) == 370
    pred45 = {r["record_id"]: r for r in read_csv(NEW45_PRED)}
    set45 = {r["record_id"]: r for r in read_csv(NEW45_SET)}
    assert len(set45) == 45 and len(pred45) == 45

    rows = []
    for rid, base in bench370.items():
        scored = pred370[rid]
        rows.append({
            "record_id": rid, "source": "frozen_benchmark_370", "subset": base["subset"],
            "tf": base["tf"], "gene": base["gene"], "label": int(base["label"]),
            "chip_score": float(scored["chip_score"]), "motif_score": float(scored["motif_score"]),
            "literature_score": float(scored["literature_score"]), "raw_score_v2.6": float(scored["raw_score"]),
        })
    for rid, base in set45.items():
        scored = pred45[rid]
        if scored["status"] != "ok":
            raise RuntimeError(f"extension record not scored: {rid}")
        rows.append({
            "record_id": rid, "source": "calib45_extension", "subset": base["subset"],
            "tf": base["tf"], "gene": base["gene"], "label": int(base["label"]),
            "chip_score": float(scored["chip_score"]), "motif_score": float(scored["motif_score"]),
            "literature_score": float(scored["literature_score"]), "raw_score_v2.6": float(scored["raw_score"]),
        })
    pairs = {(r["tf"], r["gene"]) for r in rows}
    ids = {r["record_id"] for r in rows}
    assert len(rows) == 415 and len(pairs) == 415 and len(ids) == 415, "merged corpus is not 415 unique pairs"
    assert sum(r["label"] for r in rows) == 200, "expected 200 positives"
    return rows


def stratified_split(labels, dev_n, seed):
    labels = np.asarray(labels)
    rng = np.random.default_rng(seed)
    dev_idx = []
    n_pos_target = round(int(labels.sum()) * dev_n / len(labels))
    pos = np.flatnonzero(labels == 1)
    neg = np.flatnonzero(labels == 0)
    dev_idx.extend(rng.choice(pos, size=n_pos_target, replace=False))
    dev_idx.extend(rng.choice(neg, size=dev_n - n_pos_target, replace=False))
    dev_idx = np.array(sorted(dev_idx))
    test_idx = np.setdiff1d(np.arange(len(labels)), dev_idx)
    return dev_idx, test_idx


def stratified_folds(labels, n_folds, seed):
    labels = np.asarray(labels)
    rng = np.random.default_rng(seed)
    folds = [[] for _ in range(n_folds)]
    for cls in (0, 1):
        idx = np.flatnonzero(labels == cls)
        rng.shuffle(idx)
        for i, item in enumerate(idx):
            folds[i % n_folds].append(int(item))
    return [np.array(sorted(f)) for f in folds]


def weight_grid():
    combos = []
    for a in np.arange(0.0, 1.0001, 0.05):
        for b in np.arange(0.0, 1.0001 - a, 0.05):
            c = round(1.0 - a - b, 10)
            if c < -1e-9:
                continue
            combos.append((round(float(a), 2), round(float(b), 2), round(float(c), 2)))
    return combos


def weighted_scores(rows, weights):
    chip = np.array([r["chip_score"] for r in rows])
    motif = np.array([r["motif_score"] for r in rows])
    lit = np.array([r["literature_score"] for r in rows])
    return weights[0] * chip + weights[1] * motif + weights[2] * lit


def normalize(scores, lo, hi):
    return (np.asarray(scores, dtype=float) - lo) / (hi - lo) if hi > lo else np.zeros_like(scores, dtype=float)


def best_threshold(labels, norm_scores):
    """Grid 0.00-1.00 step 0.01; maximise F1 with documented tie-breaks."""
    best_key, best_t = None, None
    labels = np.asarray(labels)
    for t in np.round(np.arange(0.0, 1.0001, 0.01), 2):
        m = confusion(labels, (norm_scores >= t).astype(int))
        key = (m["f1"], m["specificity"], m["recall_sensitivity"], -abs(t - 0.37))
        if best_key is None or key > best_key:
            best_key, best_t = key, t
    return float(best_t)


def cv_evaluate_weights(dev_rows):
    labels = np.array([r["label"] for r in dev_rows])
    folds = stratified_folds(labels, N_FOLDS, CV_SEED)
    fold_arrays = [labels[f] for f in folds]
    results = []
    all_scores = {w: weighted_scores(dev_rows, w) for w in weight_grid()}
    for w in weight_grid():
        scores = all_scores[w]
        fold_f1, fold_spec, fold_sens = [], [], []
        for k, val_idx in enumerate(folds):
            train_idx = np.concatenate([folds[j] for j in range(N_FOLDS) if j != k])
            train_scores = scores[train_idx]
            train_labels = labels[train_idx]
            lo, hi = float(train_scores.min()), float(train_scores.max())
            t = best_threshold(train_labels, normalize(train_scores, lo, hi))
            val_norm = normalize(scores[val_idx], lo, hi)
            m = confusion(fold_arrays[k], (val_norm >= t).astype(int))
            fold_f1.append(m["f1"])
            fold_spec.append(m["specificity"])
            fold_sens.append(m["recall_sensitivity"])
        results.append({
            "weights_chip_motif_literature": list(w),
            "fold_f1": [round(x, 6) for x in fold_f1],
            "mean_f1": float(np.mean(fold_f1)),
            "mean_specificity": float(np.mean(fold_spec)),
            "mean_sensitivity": float(np.mean(fold_sens)),
        })
    return results


def bootstrap_mean_ci(values, seed, n=BOOTSTRAP_N):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    means = [float(np.mean(values[rng.integers(0, len(values), size=len(values))])) for _ in range(n)]
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def select_best_weights(cv_results, min_weight=MIN_WEIGHT):
    """Eligibility (all weights >= one grid step) + one-standard-error rule.

    The CV F1 surface over the weight grid is flat (many combinations are
    statistically indistinguishable), so the plain argmax is unstable and
    favours zeroing out an entire evidence line.  Following the standard 1-SE
    rule for flat CV surfaces: anchor on the GLOBAL best mean CV F1 (over all
    231 combinations), keep every ELIGIBLE combination (each weight at least
    one grid step) within one standard error of that anchor, then select the
    combination that maximises the minimum evidence weight (the most balanced
    three-evidence contribution), with documented deterministic tie-breaks.
    """
    def sort_key(item):
        w = item["weights_chip_motif_literature"]
        dist = sum(abs(a - b) for a, b in zip(w, PUBLISHED_WEIGHTS))
        return (min(w), item["mean_f1"], item["mean_specificity"], item["mean_sensitivity"], -dist, tuple(-x for x in w))

    eligible = [r for r in cv_results if min(r["weights_chip_motif_literature"]) >= min_weight - 1e-9]
    unconstrained_best = max(cv_results, key=lambda r: r["mean_f1"])
    se = float(np.std(unconstrained_best["fold_f1"], ddof=1) / np.sqrt(N_FOLDS))
    cutoff = unconstrained_best["mean_f1"] - se
    one_se_set = [r for r in eligible if r["mean_f1"] >= cutoff]
    best = max(one_se_set, key=sort_key)
    ties = [r for r in one_se_set
            if abs(min(r["weights_chip_motif_literature"]) - min(best["weights_chip_motif_literature"])) < 1e-9
            and abs(r["mean_f1"] - best["mean_f1"]) < 1e-12]
    return best, ties, eligible, {"unconstrained_best": unconstrained_best,
                                  "se_of_global_best": se, "one_se_cutoff": cutoff, "n_within_one_se": len(one_se_set)}


def bootstrap_metrics_ci(labels, predictions, seed, n=BOOTSTRAP_N):
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    predictions = np.asarray(predictions)
    samples = []
    idx = np.arange(len(labels))
    for _ in range(n):
        take = rng.choice(idx, size=len(idx), replace=True)
        samples.append(confusion(labels[take], predictions[take]))
    out = {}
    for k in ("recall_sensitivity", "specificity", "precision", "f1", "accuracy"):
        vals = np.array([m[k] for m in samples])
        out[k] = [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))]
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = merge_415()
    labels_all = [r["label"] for r in rows]
    dev_idx, test_idx = stratified_split(labels_all, DEV_N, SPLIT_SEED)
    dev = [rows[i] for i in dev_idx]
    test = [rows[i] for i in test_idx]
    dev_labels = [r["label"] for r in dev]
    test_labels = [r["label"] for r in test]
    assert len(dev) == DEV_N and len(test) == TEST_N

    cv_results = cv_evaluate_weights(dev)
    cv_results.sort(key=lambda r: -r["mean_f1"])
    best, ties, eligible, diag = select_best_weights(cv_results)
    best["mean_f1_bootstrap_95ci"] = bootstrap_mean_ci(best["fold_f1"], BOOTSTRAP_SEED)

    weights = tuple(best["weights_chip_motif_literature"])
    dev_scores = weighted_scores(dev, weights)
    lo, hi = float(dev_scores.min()), float(dev_scores.max())
    dev_norm = normalize(dev_scores, lo, hi)
    threshold = best_threshold(dev_labels, dev_norm)
    dev_metrics = confusion(dev_labels, (dev_norm >= threshold).astype(int))

    test_scores = weighted_scores(test, weights)
    test_norm = normalize(test_scores, lo, hi)
    test_pred = (test_norm >= threshold).astype(int)
    test_metrics = confusion(test_labels, test_pred)
    test_metrics["bootstrap_95ci"] = bootstrap_metrics_ci(test_labels, test_pred, BOOTSTRAP_SEED)

    v26_pred = np.array([int(r["raw_score_v2.6"] >= 37.0) for r in test])
    v26_metrics = confusion(test_labels, v26_pred)
    v26_metrics["bootstrap_95ci"] = bootstrap_metrics_ci(test_labels, v26_pred, BOOTSTRAP_SEED + 1)

    for i, r in enumerate(rows):
        part = "development" if i in set(dev_idx.tolist()) else "test"
        r["partition"] = part
        s = weighted_scores([r], weights)[0]
        r["weighted_score"] = round(float(s), 6)
        r["weighted_score_normalized"] = round(float(normalize([s], lo, hi)[0]), 6)
        r["prediction_final"] = int(r["weighted_score_normalized"] >= threshold)
        r["prediction_v2.6_raw37"] = int(r["raw_score_v2.6"] >= 37.0)

    fields = ["record_id", "source", "subset", "partition", "label", "tf", "gene",
              "chip_score", "motif_score", "literature_score", "raw_score_v2.6",
              "weighted_score", "weighted_score_normalized", "prediction_final", "prediction_v2.6_raw37"]
    for name, subset_rows in (("merged_415_unique_split.csv", rows),
                              ("development_332.csv", dev),
                              ("test_83.csv", test)):
        with (OUT / name).open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows([{k: r.get(k, "") for k in fields} for r in subset_rows])

    with (OUT / "cv_weights_all231_with_eligibility.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["weights_chip_motif_literature", "min_weight_ge_0.05", "fold_f1", "mean_f1", "mean_specificity", "mean_sensitivity"])
        w.writeheader()
        for r in cv_results:
            w.writerow({"weights_chip_motif_literature": "/".join(f"{x:.2f}" for x in r["weights_chip_motif_literature"]),
                        "min_weight_ge_0.05": int(min(r["weights_chip_motif_literature"]) >= MIN_WEIGHT - 1e-9),
                        "fold_f1": "/".join(f"{x:.4f}" for x in r["fold_f1"]),
                        "mean_f1": r["mean_f1"], "mean_specificity": r["mean_specificity"],
                        "mean_sensitivity": r["mean_sensitivity"]})

    summary = {
        "analysis": "leakage-free 415-pair development/test re-analysis (manuscript Response 2)",
        "inputs": {
            "frozen_370": {"pairs": 370, "positives": 178, "negatives": 192},
            "extension_45": {"pairs": 45, "positives": 22, "negatives": 23, "run_id": "calib45ext_v1_20260910", "overlap_with_frozen_370": 0},
            "merged": {"pairs": 415, "positives": 200, "negatives": 215, "unique_pairs": 415},
        },
        "seeds": {"split": SPLIT_SEED, "cv_folds": CV_SEED, "bootstrap": BOOTSTRAP_SEED},
        "split": {
            "method": "stratified random (label-balanced), fixed seed",
            "development": {"n": DEV_N, "positives": int(sum(dev_labels)), "negatives": DEV_N - int(sum(dev_labels))},
            "test": {"n": TEST_N, "positives": int(sum(test_labels)), "negatives": TEST_N - int(sum(test_labels))},
        },
        "weight_selection": {
            "grid": "all (ChIP, motif, literature) combinations in 0.05 steps summing to 1.0 (231 combos evaluated)",
            "constraint": f"eligibility requires every weight >= {MIN_WEIGHT:.2f} (one grid step): all three evidence streams must contribute to the integrated score; {len(eligible)} of 231 combos eligible",
            "selection_rule": "one-standard-error rule (standard formulation): the SE is that of the global CV best over all 231 combinations; among ELIGIBLE combinations (every weight >= one grid step) whose mean CV F1 is within one SE of that anchor, select the combination maximising the minimum evidence weight (most balanced three-evidence contribution)",
            "cv": "5-fold stratified; threshold chosen per fold on training folds only (train min-max + 0.00-1.00 grid, F1 max), applied to validation fold",
            "one_se_diagnostics": {
                "global_best_weights": diag["unconstrained_best"]["weights_chip_motif_literature"],
                "global_best_mean_f1": diag["unconstrained_best"]["mean_f1"],
                "se_of_global_best_mean_f1": diag["se_of_global_best"],
                "one_se_cutoff": diag["one_se_cutoff"],
                "n_eligible_within_one_se": diag["n_within_one_se"],
            },
            "n_tied_best": len(ties),
            "selected_weights_chip_motif_literature": list(weights),
            "selected_mean_cv_f1": best["mean_f1"],
            "selected_mean_cv_f1_bootstrap_95ci": best["mean_f1_bootstrap_95ci"],
            "fold_f1": best["fold_f1"],
            "tie_break_rule": "within the 1-SE set: higher minimum weight, then higher mean fold F1, specificity, sensitivity, then distance to published 0.50/0.25/0.25, then lexicographic",
            "unconstrained_global_best_for_transparency": {
                "weights_chip_motif_literature": diag["unconstrained_best"]["weights_chip_motif_literature"],
                "mean_f1": diag["unconstrained_best"]["mean_f1"],
                "note": "excluded by the no-zero-weight constraint (a zero ChIP-seq weight removes an entire evidence line from the integrated score)",
            },
            "evidence_score_diagnostics_on_development_set": {
                "motif_dominance": "motif score separates nearly perfectly (AUC 0.923; the 60-point tier holds 151/160 positives vs 19/172 negatives), so the CV F1 surface is flat and pulls the ChIP weight towards zero",
                "chip_conflict": "ChIP is informative overall (AUC 0.755) but 5 negatives sit at chip=100 while 16 positives sit at chip=0, so any sizeable ChIP weight pushes those negatives above low-ChIP positives",
                "literature_sparsity": "literature is sparse but precision-pure (51 of 52 non-zero values are positives)",
            },
        },
        "threshold_selection": {
            "method": "min-max normalisation on the full development set, grid 0.00-1.00 step 0.01, F1 maximised",
            "normalisation_from_development": {"min": lo, "max": hi},
            "threshold": threshold,
            "development_metrics_at_selection": dev_metrics,
        },
        "final_test_evaluation": {
            "n": TEST_N, "evaluations_performed": 1,
            "weights_chip_motif_literature": list(weights),
            "threshold_normalized": threshold,
            "metrics": test_metrics,
        },
        "reference_v2.6_fixed_rule_on_same_test_set": {
            "rule": "raw_score >= 37.0 (frozen v2.6 published rule, unchanged, reported for reference)",
            "metrics": v26_metrics,
        },
    }
    (OUT / "metrics_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "selected_weights": list(weights),
        "mean_cv_f1": best["mean_f1"],
        "ci": best["mean_f1_bootstrap_95ci"],
        "fold_f1": best["fold_f1"],
        "n_tied": len(ties),
        "threshold": threshold,
        "dev_metrics": dev_metrics,
        "test_metrics": {k: v for k, v in test_metrics.items() if k != "bootstrap_95ci"},
        "test_bootstrap_95ci": test_metrics["bootstrap_95ci"],
        "v26_reference_metrics": {k: v for k, v in v26_metrics.items() if k != "bootstrap_95ci"},
    }, indent=2))


if __name__ == "__main__":
    main()
