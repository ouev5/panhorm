# -*- coding: utf-8 -*-
"""Generate single-panel PDF figures for the 9.11 deliverable folder.

Every figure is its own PDF (no a/b/c multi-panel pages), per user request:
 01 weight-grid scatter      (Blues: darker = higher mean CV F1)
 02 1-SE selection curve
 03 threshold grid search    (winner configuration, development)
 04 score distribution       (development)
 05 score distribution       (independent test)
 06 ROC                      (final model + frozen v2.6 reference)
 07 precision-recall
 08 confusion matrix
 09 evidence-score boxplots
 10 univariate AUC           (LINE plot, per request)
 11 five-fold CV F1          (BOX plot for the two tied configurations, per request)
 12 two-configuration threshold comparison (development)
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import auc, confusion_matrix, precision_recall_curve, roc_curve

FINAL_DIR = Path(r"E:\9.2_manuscript\GeneReg_415独立重跑_1SE均衡权重_20260910")
OUT = Path(r"C:\Users\huang'xing'xing\Desktop\9.7\9.11")

WINNER = (0.10, 0.70, 0.20)
RUNNER = (0.10, 0.65, 0.25)
THRESHOLD_WINNER = 0.41
THRESHOLD_RUNNER = 0.36
GRID_MEAN_F1 = {0.70: 0.6004, 0.65: 0.5704}  # by motif weight: winner/runner
CUTTOFF_1SE = 0.9048
PUBLISHED = (0.50, 0.25, 0.25)

C_POS, C_NEG, C_SEL, C_REF, C_GREY = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#999999"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})
MM = 1 / 25.4


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def load():
    rows = read_csv(FINAL_DIR / "merged_415_unique_split.csv")
    for r in rows:
        r["label"] = int(r["label"])
        for k in ("chip_score", "motif_score", "literature_score", "raw_score_v2.6", "weighted_score", "weighted_score_normalized"):
            r[k] = float(r[k])
    cv = read_csv(FINAL_DIR / "cv_weights_all231_with_eligibility.csv")
    for r in cv:
        r["w"] = [float(x) for x in r["weights_chip_motif_literature"].split("/")]
        r["folds"] = [float(x) for x in r["fold_f1"].split("/")]
        r["mean_f1"] = float(r["mean_f1"])
        r["eligible"] = r["min_weight_ge_0.05"] == "1"
    return rows, cv


def fig_weight_grid(cv):
    fig, ax = plt.subplots(figsize=(120 * MM, 92 * MM))
    elig = [r for r in cv if r["eligible"]]
    inelig = [r for r in cv if not r["eligible"]]
    ax.scatter([r["w"][1] for r in inelig], [r["w"][0] for r in inelig], s=12, c="#CCCCCC", alpha=0.55, label="Ineligible (a zero weight)", linewidths=0)
    vals = [r["mean_f1"] for r in elig]
    sc = ax.scatter([r["w"][1] for r in elig], [r["w"][0] for r in elig], s=20, c=vals, cmap="Blues", vmin=min(vals), vmax=max(vals), label="Eligible (all weights \u2265 0.05)", linewidths=0.15, edgecolors="none")
    ax.scatter([WINNER[1]], [WINNER[0]], s=110, facecolors="none", edgecolors=C_NEG, linewidths=1.8, zorder=5, clip_on=False, label="Selected 0.10/0.65/0.25")
    ax.scatter([PUBLISHED[1]], [PUBLISHED[0]], marker="X", s=55, c=C_NEG, zorder=5, clip_on=False, label="Original 0.50/0.25/0.25 (0.756)")
    ax.set_xlabel("Motif weight")
    ax.set_ylabel("ChIP-seq weight")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="upper right", frameon=True, framealpha=0.95, edgecolor="0.8")
    cb = plt.colorbar(sc, ax=ax, shrink=0.9)
    cb.set_label("Mean five-fold CV F1")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_01_weight_grid.pdf")
    plt.close(fig)


def fig_1se(cv):
    fig, ax = plt.subplots(figsize=(120 * MM, 85 * MM))
    groups = {}
    for r in cv:
        groups.setdefault(round(r["w"][0], 2), []).append(r["mean_f1"])
    xs = sorted(groups)
    ax.plot(xs, [max(groups[x]) for x in xs], "-o", ms=3.5, color=C_POS)
    ax.axhline(CUTTOFF_1SE, ls="--", lw=1, color=C_GREY)
    ax.text(0.985, CUTTOFF_1SE - 0.012, "1-SE cut-off (0.905)", color=C_GREY, ha="right", va="top")
    ax.axvline(WINNER[0], ls=":", lw=1, color=C_SEL)
    ax.scatter([WINNER[0]], [0.9083], s=80, facecolors="none", edgecolors=C_SEL, linewidths=1.6, zorder=5, clip_on=False)
    ax.annotate("Selected 0.10/0.65/0.25 (0.908)", (WINNER[0], 0.9083), xytext=(0.42, 0.868), color=C_SEL,
                arrowprops=dict(arrowstyle="-", color=C_SEL, lw=0.7))
    ax.scatter([0.0], [0.9225], marker="v", s=32, color=C_REF, zorder=5, clip_on=False)
    ax.annotate("Global best 0.00/0.60/0.40 (0.923)\nzero ChIP weight \u2014 excluded", (0.0, 0.9225), xytext=(0.13, 0.938), color=C_REF, va="top",
                arrowprops=dict(arrowstyle="-", color=C_REF, lw=0.7))
    ax.scatter([PUBLISHED[0]], [0.7561], marker="X", s=50, color=C_NEG, zorder=5, clip_on=False)
    ax.annotate("Original 0.50/0.25/0.25 (0.756)", (PUBLISHED[0], 0.7561), xytext=(0.47, 0.782), color=C_NEG, ha="right",
                arrowprops=dict(arrowstyle="-", color=C_NEG, lw=0.7))
    ax.set_xlabel("ChIP-seq weight")
    ax.set_ylabel("Best mean five-fold CV F1")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.73, 0.945)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_02_1SE_selection.pdf")
    plt.close(fig)


def f1_curve(rows, w):
    y = np.array([r["label"] for r in rows])
    s = w[0] * np.array([r["chip_score"] for r in rows]) + w[1] * np.array([r["motif_score"] for r in rows]) + w[2] * np.array([r["literature_score"] for r in rows])
    lo, hi = s.min(), s.max()
    norm = (s - lo) / (hi - lo)
    ts = np.round(np.arange(0.0, 1.0001, 0.01), 2)
    f1s = []
    for t in ts:
        pred = (norm >= t).astype(int)
        tp = np.sum((y == 1) & (pred == 1)); fp = np.sum((y == 0) & (pred == 1)); fn = np.sum((y == 1) & (pred == 0))
        prec = tp / (tp + fp) if tp + fp else 0.0
        sens = tp / (tp + fn) if tp + fn else 0.0
        f1s.append(2 * prec * sens / (prec + sens) if prec + sens else 0.0)
    return ts, np.array(f1s), lo, hi


def fig_threshold_search(dev):
    ts, f1s, lo, hi = f1_curve(dev, WINNER)
    fig, ax = plt.subplots(figsize=(100 * MM, 75 * MM))
    ax.plot(ts, f1s, color=C_POS, lw=1.5)
    ax.axvline(THRESHOLD_WINNER, ls="--", lw=1, color=C_SEL)
    ax.annotate(f"F1-maximising threshold = {THRESHOLD_WINNER:.2f}\n(F1 = {f1s.max():.3f})", (THRESHOLD_WINNER, f1s.max()),
                textcoords="offset points", xytext=(12, -34), color=C_SEL)
    ax.set_xlabel("Normalised threshold")
    ax.set_ylabel("F1 score")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.0)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_03_threshold_search.pdf")
    plt.close(fig)


def winner_norm(rows, dev_rows):
    w = WINNER
    sd = w[0] * np.array([r["chip_score"] for r in dev_rows]) + w[1] * np.array([r["motif_score"] for r in dev_rows]) + w[2] * np.array([r["literature_score"] for r in dev_rows])
    lo, hi = float(sd.min()), float(sd.max())
    s = w[0] * np.array([r["chip_score"] for r in rows]) + w[1] * np.array([r["motif_score"] for r in rows]) + w[2] * np.array([r["literature_score"] for r in rows])
    return (s - lo) / (hi - lo)


def fig_dist(rows, name, fname, dev_rows):
    fig, ax = plt.subplots(figsize=(100 * MM, 70 * MM))
    s = winner_norm(rows, dev_rows)
    y = np.array([r["label"] for r in rows])
    bins = np.linspace(0, 1, 26)
    centers = (bins[:-1] + bins[1:]) / 2
    width = (bins[1] - bins[0]) * 0.42
    hn, _ = np.histogram(s[y == 0], bins=bins)
    hp, _ = np.histogram(s[y == 1], bins=bins)
    ax.bar(centers - width / 2, hn, width=width, color=C_NEG, edgecolor="white", linewidth=0.4, label=f"Negative (n = {(y == 0).sum()})")
    ax.bar(centers + width / 2, hp, width=width, color=C_POS, edgecolor="white", linewidth=0.4, label=f"Positive (n = {(y == 1).sum()})")
    ax.axvline(THRESHOLD_WINNER, ls="--", lw=1.2, color=C_SEL)
    ax.text(THRESHOLD_WINNER + 0.025, ax.get_ylim()[1] * 0.52, f"threshold\n{THRESHOLD_WINNER:.2f}", color=C_SEL)
    ax.set_xlabel("Normalised weighted score")
    ax.set_ylabel("Number of pairs")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper center", frameon=True, framealpha=0.95, edgecolor="0.8")
    ax.set_title(name, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def fig_roc_pr(test, dev_rows):
    y = np.array([r["label"] for r in test])
    s_final = winner_norm(test, dev_rows)
    s_ref = np.array([r["raw_score_v2.6"] for r in test])
    fpr1, tpr1, _ = roc_curve(y, s_final)
    fpr2, tpr2, _ = roc_curve(y, s_ref)
    fig, ax = plt.subplots(figsize=(90 * MM, 90 * MM))
    ax.plot(fpr1, tpr1, color=C_POS, lw=1.6, label=f"Final model (AUC = {auc(fpr1, tpr1):.3f})")
    ax.plot(fpr2, tpr2, color=C_REF, lw=1.3, ls="--", label=f"Frozen v2.6 raw score (AUC = {auc(fpr2, tpr2):.3f})")
    ax.plot([0, 1], [0, 1], ls=":", lw=0.8, color=C_GREY)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.legend(loc="lower right", frameon=True, framealpha=0.95, edgecolor="0.8")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_06_roc.pdf")
    plt.close(fig)

    prec1, rec1, _ = precision_recall_curve(y, s_final)
    prec2, rec2, _ = precision_recall_curve(y, s_ref)
    fig, ax = plt.subplots(figsize=(90 * MM, 90 * MM))
    ax.plot(rec1, prec1, color=C_POS, lw=1.6, label=f"Final model (AP = {auc(rec1, prec1):.3f})")
    ax.plot(rec2, prec2, color=C_REF, lw=1.3, ls="--", label=f"Frozen v2.6 raw score (AP = {auc(rec2, prec2):.3f})")
    ax.axhline(y.mean(), ls=":", lw=0.8, color=C_GREY)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower left", frameon=True, framealpha=0.95, edgecolor="0.8")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_07_pr.pdf")
    plt.close(fig)


def fig_confusion(test, dev_rows):
    y = np.array([r["label"] for r in test])
    pred = (winner_norm(test, dev_rows) >= THRESHOLD_WINNER).astype(int)
    cm = confusion_matrix(y, pred)
    fig, ax = plt.subplots(figsize=(80 * MM, 80 * MM))
    ax.imshow(cm, cmap="Blues", vmin=0, vmax=cm.max() * 1.15)
    for (i, j), v in np.ndenumerate(cm):
        ax.text(j, i, str(v), ha="center", va="center", color="white" if v > cm.max() * 0.6 else "black", fontsize=13)
    ax.set_xticks([0, 1], ["Predicted \u2212", "Predicted +"])
    ax.set_yticks([0, 1], ["Reference \u2212", "Reference +"])
    ax.spines[:].set_visible(True)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_08_confusion_matrix.pdf")
    plt.close(fig)


def fig_evidence_box(dev):
    y = np.array([r["label"] for r in dev])
    data = [np.array([r["chip_score"] for r in dev]), np.array([r["motif_score"] for r in dev]), np.array([r["literature_score"] for r in dev])]
    names = ["ChIP-seq", "Motif", "Literature"]
    fig, ax = plt.subplots(figsize=(100 * MM, 75 * MM))
    for i, (d, name) in enumerate(zip(data, names)):
        bp = ax.boxplot([d[y == 0], d[y == 1]], positions=[i - 0.18, i + 0.18], widths=0.3, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=0.9))
        for patch, color in zip(bp["boxes"], [C_NEG, C_POS]):
            patch.set_facecolor(color)
            patch.set_alpha(0.65)
        rng = np.random.default_rng(7 + i)
        for j, grp in ((0, y == 0), (1, y == 1)):
            xs = rng.normal(i - 0.18 + 0.18 * j, 0.025, grp.sum())
            ax.scatter(xs, d[grp], s=3, color=C_NEG if j == 0 else C_POS, alpha=0.35, linewidths=0)
    ax.set_xticks(range(3), names)
    ax.set_ylabel("Evidence score (0\u2013100)")
    ax.set_ylim(-6, 130)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=C_NEG, alpha=0.65, label="Negative"), Patch(facecolor=C_POS, alpha=0.65, label="Positive")],
              loc="upper right", frameon=True, framealpha=0.95, edgecolor="0.8")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_09_evidence_boxplot.pdf")
    plt.close(fig)


def auc_score(d, y):
    order = np.argsort(d)
    ranks = np.empty(len(d))
    ranks[order] = np.arange(1, len(d) + 1)
    for v in np.unique(d):
        m = d == v
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    n1, n0 = y.sum(), (y == 0).sum()
    return (ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def fig_auc_line(dev):
    y = np.array([r["label"] for r in dev])
    data = [np.array([r["chip_score"] for r in dev]), np.array([r["motif_score"] for r in dev]), np.array([r["literature_score"] for r in dev])]
    names = ["ChIP-seq", "Motif", "Literature"]
    aucs = [auc_score(d, y) for d in data]
    fig, ax = plt.subplots(figsize=(100 * MM, 75 * MM))
    ax.plot(range(3), aucs, "-o", color=C_POS, lw=1.6, ms=6)
    for i, v in enumerate(aucs):
        ax.annotate(f"{v:.3f}", (i, v), textcoords="offset points", xytext=(8, 4), color="black")
    ax.axhline(0.5, ls=":", lw=0.9, color=C_GREY)
    ax.text(2.02, 0.505, "chance (0.5)", color=C_GREY, fontsize=7, ha="right")
    ax.set_xticks(range(3), names)
    ax.set_xlim(-0.25, 2.25)
    ax.set_ylim(0.4, 1.05)
    ax.set_ylabel("Univariate AUC")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_10_univariate_auc_line.pdf")
    plt.close(fig)


def fig_cv_box(cv):
    folds_by_w = {"/".join(f"{x:.2f}" for x in r["w"]): r["folds"] for r in cv}
    a = folds_by_w["/".join(f"{x:.2f}" for x in WINNER)]
    b = folds_by_w["/".join(f"{x:.2f}" for x in RUNNER)]
    fig, ax = plt.subplots(figsize=(100 * MM, 80 * MM))
    bp = ax.boxplot([a, b], positions=[1, 2], widths=0.42, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.0))
    for patch, color in zip(bp["boxes"], [C_SEL, C_REF]):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
    rng = np.random.default_rng(11)
    for i, v in enumerate((a, b), 1):
        xs = rng.normal(i, 0.045, len(v))
        ax.scatter(xs, v, s=22, color=C_SEL if i == 1 else C_REF, alpha=0.85, zorder=5, linewidths=0)
        ax.scatter([i], [np.mean(v)], marker="D", s=30, color="white", edgecolor="black", zorder=6)
        ax.annotate(f"mean {np.mean(v):.3f}", (i, np.mean(v)), textcoords="offset points", xytext=(12, -3), fontsize=8)
    ax.set_xticks([1, 2], ["/".join(f"{x:.2f}" for x in WINNER) + "\n(selected)", "/".join(f"{x:.2f}" for x in RUNNER) + "\n(tied alternative)"])
    ax.set_ylabel("Validation-fold F1")
    ax.set_ylim(0.85, 0.97)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_11_cv_f1_boxplot.pdf")
    plt.close(fig)


def fig_two_config(dev):
    ts1, f1s1, _, _ = f1_curve(dev, WINNER)
    ts2, f1s2, _, _ = f1_curve(dev, RUNNER)
    fig, ax = plt.subplots(figsize=(120 * MM, 88 * MM))
    ax.plot(ts1, f1s1, color=C_SEL, lw=1.7, label="0.10/0.70/0.20 — selected (threshold 0.41)")
    ax.plot(ts2, f1s2, color=C_REF, lw=1.4, ls="--", label="0.10/0.65/0.25 — tied at peak (threshold 0.36)")
    ax.axvline(THRESHOLD_WINNER, ls=":", lw=1, color=C_SEL)
    ax.axvline(THRESHOLD_RUNNER, ls=":", lw=1, color=C_REF)
    ax.scatter([THRESHOLD_WINNER, THRESHOLD_RUNNER], [f1s1.max(), f1s2.max()], s=48, color=[C_SEL, C_REF], zorder=5, clip_on=False)
    ax.annotate("identical peak: F1 = 0.908", (0.385, f1s1.max()), textcoords="offset points", xytext=(-6, 10), ha="center", color="black")
    ax.text(0.985, 0.80, "Threshold-robustness (mean F1 over the full grid):\n0.600 for 0.10/0.70/0.20 vs 0.570 for 0.10/0.65/0.25\n→ 0.10/0.70/0.20 selected", ha="right", va="top", color="black")
    ax.set_xlabel("Normalised threshold")
    ax.set_ylabel("Development-set F1")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.0)
    ax.legend(loc="lower left", frameon=True, framealpha=0.95, edgecolor="0.8")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_12_two_config_threshold.pdf")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, cv = load()
    dev = [r for r in rows if r["partition"] == "development"]
    test = [r for r in rows if r["partition"] == "test"]
    fig_weight_grid(cv)
    fig_1se(cv)
    fig_threshold_search(dev)
    fig_dist(dev, "Development set (n = 332)", "Figure_04_score_dist_development.pdf", dev)
    fig_dist(test, "Independent test set (n = 83)", "Figure_05_score_dist_test.pdf", dev)
    fig_roc_pr(test, dev)
    fig_confusion(test, dev)
    fig_evidence_box(dev)
    fig_auc_line(dev)
    fig_cv_box(cv)
    fig_two_config(dev)
    for p in sorted(OUT.glob("Figure_*.pdf")):
        print("written:", p.name, p.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
