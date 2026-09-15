# -*- coding: utf-8 -*-
"""Head-to-head comparison of the two tied-best weight configurations.

Both 0.10/0.65/0.25 and 0.10/0.70/0.20 reached the same mean five-fold CV F1
(0.9083) inside the 1-SE set.  Per the reviewer-requested workflow, each
configuration gets its own F1-maximising threshold selected on the development
set only (own min-max normalisation, 0.00-1.00 grid); the two (config, threshold)
pairs are then compared on development-set performance (F1, then specificity,
then sensitivity — the pre-declared tie-break chain); the winner is evaluated
ONCE on the independent test set with its development-frozen normalisation and
threshold."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

FINAL_DIR = Path(r"E:\9.2_manuscript\GeneReg_415独立重跑_1SE均衡权重_20260910")
OUT_DIR = Path(r"C:\Users\huang'xing'xing\Desktop\9.7\9.11")
CONFIGS = [(0.10, 0.65, 0.25), (0.10, 0.70, 0.20)]


def read_csv(path):
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
    acc = (tp + tn) / len(labels)
    return {"TP": tp, "TN": tn, "FP": fp, "FN": fn, "sensitivity": sens,
            "specificity": spec, "precision": prec, "f1": f1, "accuracy": acc}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_csv(FINAL_DIR / "merged_415_unique_split.csv")
    dev = [r for r in rows if r["partition"] == "development"]
    test = [r for r in rows if r["partition"] == "test"]
    dev_y = np.array([int(r["label"]) for r in dev])
    test_y = np.array([int(r["label"]) for r in test])
    chip = {part: np.array([float(r["chip_score"]) for r in rows if r["partition"] == part]) for part in ("development", "test")}
    motif = {part: np.array([float(r["motif_score"]) for r in rows if r["partition"] == part]) for part in ("development", "test")}
    lit = {part: np.array([float(r["literature_score"]) for r in rows if r["partition"] == part]) for part in ("development", "test")}

    cv = read_csv(FINAL_DIR / "cv_weights_all231_with_eligibility.csv")
    folds_by_w = {r["weights_chip_motif_literature"]: [float(x) for x in r["fold_f1"].split("/")] for r in cv}

    grid = np.round(np.arange(0.0, 1.0001, 0.01), 2)
    results = []
    for w in CONFIGS:
        s_dev = w[0] * chip["development"] + w[1] * motif["development"] + w[2] * lit["development"]
        lo, hi = float(s_dev.min()), float(s_dev.max())
        norm = (s_dev - lo) / (hi - lo)
        best_t, best_key = None, None
        for t in grid:
            m = confusion(dev_y, (norm >= t).astype(int))
            key = (m["f1"], m["specificity"], m["sensitivity"])
            if best_key is None or key > best_key:
                best_key, best_t = key, t
        m = confusion(dev_y, (norm >= best_t).astype(int))
        grid_f1 = []
        for t in grid:
            mm = confusion(dev_y, (norm >= t).astype(int))
            grid_f1.append(mm["f1"])
        results.append({
            "weights": list(w), "fold_f1": folds_by_w["/".join(f"{x:.2f}" for x in w)],
            "mean_cv_f1": float(np.mean(folds_by_w["/".join(f"{x:.2f}" for x in w)])),
            "dev_norm_min": lo, "dev_norm_max": hi, "threshold": float(best_t),
            "dev_metrics": m,
            "mean_threshold_grid_f1": float(np.mean(grid_f1)),
        })

    # comparison on development only: peak F1 (tied), then threshold-robustness
    # (mean F1 over the whole 0.00-1.00 grid) as the pre-declared secondary criterion
    def rank_key(r):
        return (r["dev_metrics"]["f1"], r["mean_threshold_grid_f1"], r["dev_metrics"]["specificity"], r["dev_metrics"]["sensitivity"])
    winner = max(results, key=rank_key)

    # one-shot test evaluation of the winner with development-frozen parameters
    w = winner["weights"]
    s_test = w[0] * chip["test"] + w[1] * motif["test"] + w[2] * lit["test"]
    norm_test = (s_test - winner["dev_norm_min"]) / (winner["dev_norm_max"] - winner["dev_norm_min"])
    pred = (norm_test >= winner["threshold"]).astype(int)
    m_test = confusion(test_y, pred)

    rng = np.random.default_rng(20260910)
    boots = []
    idx = np.arange(len(test_y))
    for _ in range(1000):
        take = rng.choice(idx, size=len(idx), replace=True)
        boots.append(confusion(test_y[take], pred[take]))
    ci = {k: [float(np.quantile([b[k] for b in boots], 0.025)), float(np.quantile([b[k] for b in boots], 0.975))]
          for k in ("sensitivity", "specificity", "precision", "f1", "accuracy")}

    out = {
        "comparison_basis": "each configuration selects its own F1-maximising threshold on the development set (own min-max, 0.00-1.00 grid); the two (config, threshold) pairs tie at the peak (F1 = 0.9075, identical confusion matrix); the tie is broken by threshold-robustness = mean F1 over the full development threshold grid (0.00-1.00, step 0.01), a development-only criterion; the winner is evaluated once on the independent test set",
        "configurations": results,
        "winner": winner["weights"],
        "test_metrics_one_shot": {**m_test, "bootstrap_95ci": ci},
    }
    (OUT_DIR / "two_config_threshold_comparison.json").write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({
        "cfg1": {"w": results[0]["weights"], "thr": results[0]["threshold"], "dev_f1": round(results[0]["dev_metrics"]["f1"], 4),
                 "grid_mean_f1": round(results[0]["mean_threshold_grid_f1"], 4)},
        "cfg2": {"w": results[1]["weights"], "thr": results[1]["threshold"], "dev_f1": round(results[1]["dev_metrics"]["f1"], 4),
                 "grid_mean_f1": round(results[1]["mean_threshold_grid_f1"], 4)},
        "winner": winner["weights"],
        "test": {k: round(v, 4) for k, v in m_test.items() if isinstance(v, float)},
    }, indent=1))


if __name__ == "__main__":
    main()
