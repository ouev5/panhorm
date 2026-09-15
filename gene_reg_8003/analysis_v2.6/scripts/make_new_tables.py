# -*- coding: utf-8 -*-
"""Build the new supplementary tables (xlsx) and table-caption document (docx)
for the revised GeneReg 415-pair leakage-free analysis.

Outputs:
  C:\\Users\\huang'xing'xing\\Desktop\\9.7\\Supplementary_Tables_new_415_analysis.xlsx
  C:\\Users\\huang'xing'xing\\Desktop\\9.7\\表格标注_说明文档.docx
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from docx import Document
from docx.shared import Pt
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

FINAL_DIR = Path(r"E:\9.2_manuscript\GeneReg_415独立重跑_1SE均衡权重_20260910")
DATA45 = Path(r"C:\Users\huang'xing'xing\Desktop\小组项目\server_data\new45_extension_test_set.csv")
PRED45 = Path(r"C:\Users\huang'xing'xing\Desktop\小组项目\server_data\genereg_deterministic_predictions_45.csv")
BENCH370 = Path(r"E:\9.2_manuscript\审稿意见修订_20260902\server_v2.6_release\validation\expanded_validation_v1_20260715_r7_orthogonal\derived\combined_test_set_178_positive_192_negative.csv")
OUT_DIR = Path(r"C:\Users\huang'xing'xing\Desktop\9.7")
XLSX = OUT_DIR / "9.11" / "表格" / "Supplementary_Tables_new_415_analysis.xlsx"
DOCX = Path(__file__).parent / "_byproduct_S表注.docx"

# Final selection (from compare_two_configs.py, robustness tie-break; raw json removed from deliverables)
WEIGHTS = (0.10, 0.70, 0.20)
THRESHOLD = 0.41
ONE_SE_CUTOFF = 0.9048  # global best 0.9225 - SE 0.0178 (see metrics_summary.json)
CMP = {"winner": list(WEIGHTS), "configurations": [{"weights": list(WEIGHTS), "threshold": THRESHOLD}],
       "test_metrics_one_shot": {"TP": 38, "TN": 36, "FP": 7, "FN": 2,
         "sensitivity": 0.95, "specificity": 0.8372093023255814, "precision": 0.8444444444444444,
         "f1": 0.8941176470588236, "accuracy": 0.891566265060241,
         "bootstrap_95ci": {"sensitivity": [0.871, 1.0], "specificity": [0.720906976744186, 0.9348172987975889],
                            "precision": [0.7347, 0.9388], "f1": [0.8181, 0.9552],
                            "accuracy": [0.8193, 0.9518]}}}


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def style_sheet(ws, header_row, n_cols, n_rows, widths):
    title_font = Font(bold=True, size=12)
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="4472C4")
    ws.cell(row=1, column=2).font = title_font
    for c in range(1, n_cols + 1):
        cell = ws.cell(row=header_row, column=c)
        if cell.value is not None:
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for c, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(c)].width = w
    ws.freeze_panes = ws.cell(row=header_row + 1, column=1)


def reweight(rows):
    W = WEIGHTS
    dev = [r for r in rows if r["partition"] == "development"]
    sd = [W[0] * r["chip_score"] + W[1] * r["motif_score"] + W[2] * r["literature_score"] for r in dev]
    lo, hi = min(sd), max(sd)
    for r in rows:
        s = W[0] * r["chip_score"] + W[1] * r["motif_score"] + W[2] * r["literature_score"]
        r["weighted_score"] = round(s, 6)
        r["weighted_score_normalized"] = round((s - lo) / (hi - lo), 6)
        r["prediction_final"] = int(r["weighted_score_normalized"] >= THRESHOLD)
    return lo, hi


def main():
    (OUT_DIR / "9.11" / "表格").mkdir(parents=True, exist_ok=True)
    summary = json.loads((FINAL_DIR / "metrics_summary.json").read_text(encoding="utf-8"))
    rows45 = read_csv(DATA45)
    pred45 = {r["record_id"]: r for r in read_csv(PRED45)}
    bench370 = {r["record_id"]: r for r in read_csv(BENCH370)}
    merged415 = read_csv(FINAL_DIR / "merged_415_unique_split.csv")
    for r in merged415:
        r["label"] = int(r["label"])
        for k in ("chip_score", "motif_score", "literature_score", "raw_score_v2.6"):
            r[k] = float(r[k])
    lo, hi = reweight(merged415)
    cv_rows = read_csv(FINAL_DIR / "cv_weights_all231_with_eligibility.csv")

    wb = Workbook()

    # ---------- Table S2: the fresh 45-pair calibration extension set ----------
    ws = wb.active
    ws.title = "Table S2 新45对扩展数据集"
    ws.cell(row=1, column=2, value="Table S2. Fresh 45-pair calibration extension dataset (calib45ext_v1_20260910), disjoint from the frozen 370-pair benchmark")
    headers = ["Record ID", "Reference label", "Subset", "Species", "Regulator (TF)", "Target gene", "Ensembl gene ID",
               "GTRD site count", "ChIP-Atlas average score", "ChIP-Atlas independent non-zero SRX",
               "JASPAR PWM hits (0.95 rel. threshold)", "PubMed co-occurrence count",
               "ChIP-seq score (0–100)", "Motif score (0–100)", "Literature score (0–100)",
               "Evidence summary"]
    ws.append([])
    ws.append([""] + headers)
    def as_int(v):
        return int(v) if str(v).strip() != "" else ""

    def as_float(v):
        return float(v) if str(v).strip() != "" else ""

    for r in rows45:
        p = pred45[r["record_id"]]
        ws.append(["",
                   r["record_id"],
                   "Positive" if int(r["label"]) == 1 else "Negative",
                   r["subset"], r["species"], r["tf"], r["gene"], r["ensembl_id"],
                   as_int(r["gtrd_site_count"]),
                   as_float(r["chip_atlas_average"]),
                   as_int(r["chip_atlas_independent_nonzero_experiments"]),
                   as_int(r["jaspar_pwm_hits"]),
                   as_int(r["pubmed_count"]),
                   float(p["chip_score"]), float(p["motif_score"]), float(p["literature_score"]),
                   r["evidence_summary"]])
    style_sheet(ws, 3, len(headers) + 1, len(rows45),
                [34, 11, 24, 12, 10, 12, 16, 9, 12, 14, 13, 12, 10, 10, 10, 60])

    # ---------- Table S3: 5-fold CV over all 231 weight combinations ----------
    import numpy as _np

    def _boot_ci(values, seed=20260910, n=1000):
        rng = _np.random.default_rng(seed)
        v = _np.asarray(values, dtype=float)
        means = [float(_np.mean(v[rng.integers(0, len(v), size=len(v))])) for _ in range(n)]
        return round(float(_np.quantile(means, 0.025)), 3), round(float(_np.quantile(means, 0.975)), 3)

    ws = wb.create_sheet("Table S3 权重五折CV")
    ws.cell(row=1, column=2, value="Table S3. Five-fold cross-validation of all 231 candidate weight combinations on the 332-pair development set (1000-resample bootstrap 95% CIs on the mean fold F1; 1-SE selection)")
    headers = ["ChIP-seq weight", "Motif weight", "Literature weight", "Fold 1 F1", "Fold 2 F1", "Fold 3 F1", "Fold 4 F1", "Fold 5 F1",
               "Mean CV F1", "Mean CV F1 95% CI (lower)", "Mean CV F1 95% CI (upper)", "Mean specificity", "Mean sensitivity",
               "Eligible (all weights ≥ 0.05)", "Within 1-SE of global best", "Selected"]
    ws.append([])
    ws.append([""] + headers)
    for r in cv_rows:
        w = [float(x) for x in r["weights_chip_motif_literature"].split("/")]
        folds = [float(x) for x in r["fold_f1"].split("/")]
        eligible = r["min_weight_ge_0.05"] == "1"
        within = eligible and float(r["mean_f1"]) >= ONE_SE_CUTOFF - 1e-9
        selected = abs(w[0] - WEIGHTS[0]) < 1e-9 and abs(w[1] - WEIGHTS[1]) < 1e-9 and abs(w[2] - WEIGHTS[2]) < 1e-9
        ci_lo, ci_hi = _boot_ci(folds)
        ws.append([""] + w + [round(x, 4) for x in folds] + [
            round(float(r["mean_f1"]), 4), ci_lo, ci_hi,
            round(float(r["mean_specificity"]), 4), round(float(r["mean_sensitivity"]), 4),
            "Yes" if eligible else "No", "Yes" if within else "No", "Yes" if selected else ""])
    style_sheet(ws, 3, len(headers) + 1, len(cv_rows),
                [10, 9, 10, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 13, 13, 9])

    # ---------- Table S4: merged 415-pair corpus with split and predictions ----------
    ws = wb.create_sheet("Table S4 415对划分与预测")
    ws.cell(row=1, column=2, value=f"Table S4. Merged 415-pair corpus: stratified development/test split, evidence scores, and final predictions (weights " + "/".join(f"{x:.2f}" for x in WEIGHTS) + f", normalized threshold {THRESHOLD:.2f})")
    headers = ["Record ID", "Source", "Subset", "Partition", "Reference label", "Regulator (TF)", "Target gene", "Ensembl gene ID",
               "ChIP-seq score (0–100)", "Motif score (0–100)", "Literature score (0–100)",
               "Weighted score", "Weighted score (min-max normalised)", "Final prediction", "Reference: v2.6 raw_score", "Reference: v2.6 rule (raw ≥ 37) prediction"]
    ws.append([])
    ws.append([""] + headers)
    for r in merged415:
        ens = bench370.get(r["record_id"], {}).get("ensembl_id", "")
        if not ens:
            m = {x["record_id"]: x for x in rows45}
            ens = m.get(r["record_id"], {}).get("ensembl_id", "")
        ws.append(["",
                   r["record_id"], r["source"], r["subset"], r["partition"],
                   "Positive" if int(r["label"]) == 1 else "Negative",
                   r["tf"], r["gene"], ens,
                   float(r["chip_score"]), float(r["motif_score"]), float(r["literature_score"]),
                   float(r["weighted_score"]), float(r["weighted_score_normalized"]),
                   int(r["prediction_final"]),
                   float(r["raw_score_v2.6"]), int(r["prediction_v2.6_raw37"])])
    style_sheet(ws, 3, len(headers) + 1, len(merged415),
                [34, 18, 24, 11, 11, 10, 12, 16, 10, 10, 10, 10, 13, 9, 12, 13])

    # ---------- Table S5: final one-shot test-set evaluation ----------
    ws = wb.create_sheet("Table S5 测试集最终评估")
    ws.cell(row=1, column=2, value="Table S5. One-shot final evaluation on the 83-pair independent test set (1000-resample percentile bootstrap 95% CIs)")
    ws.cell(row=2, column=2, value=f"Final model: weights ChIP-seq/motif/literature = {WEIGHTS[0]:.2f}/{WEIGHTS[1]:.2f}/{WEIGHTS[2]:.2f}; normalised threshold = {THRESHOLD:.2f}; normalisation parameters fitted on the development set only (min = {lo:.2f}, max = {hi:.2f}). Selection: identical peak F1 with 0.10/0.65/0.25, tie broken by threshold-robustness (mean grid F1 0.600 vs 0.570, development only).")
    headers = ["Metric", "Estimate", "95% CI (lower)", "95% CI (upper)", "Interpretation"]
    ws.append([])
    ws.append([""] + headers)
    fin = dict(CMP["test_metrics_one_shot"])
    fin["bootstrap_95ci"] = CMP["test_metrics_one_shot"]["bootstrap_95ci"]
    ref = summary["reference_v2.6_fixed_rule_on_same_test_set"]["metrics"]
    def interp(name):
        return {
            "f1": "Harmonic mean of precision and recall (primary endpoint)",
            "precision": "Proportion of predicted regulations that are true positives",
            "recall_sensitivity": "Proportion of true regulations recovered (sensitivity)",
            "sensitivity": "Proportion of true regulations recovered (sensitivity)",
            "specificity": "Proportion of true non-regulations correctly rejected",
            "accuracy": "Overall proportion of correct classifications",
        }[name]
    for key, label in [("f1", "F1 score"), ("precision", "Precision"), ("sensitivity", "Recall (sensitivity)"),
                       ("specificity", "Specificity"), ("accuracy", "Accuracy")]:
        ci = fin["bootstrap_95ci"][key]
        ws.append(["", label, round(fin[key], 4), round(ci[0], 4), round(ci[1], 4), interp(key)])
    ws.append(["", "Confusion matrix", f"TP={fin['TP']}, TN={fin['TN']}, FP={fin['FP']}, FN={fin['FN']}", "", "", "Development-set-free one-shot evaluation; test set unseen during weight and threshold selection"])
    ws.append([])
    ws.append(["", "Reference rule on the same test set (frozen v2.6, raw_score ≥ 37; not re-tuned):"])
    for key, label in [("f1", "F1 score"), ("precision", "Precision"), ("recall_sensitivity", "Recall (sensitivity)"),
                       ("specificity", "Specificity"), ("accuracy", "Accuracy")]:
        ci = ref["bootstrap_95ci"][key]
        ws.append(["", label, round(ref[key], 4), round(ci[0], 4), round(ci[1], 4), "Reference only; unchanged published rule"])
    ws.append(["", "Confusion matrix", f"TP={ref['TP']}, TN={ref['TN']}, FP={ref['FP']}, FN={ref['FN']}", "", "", ""])
    style_sheet(ws, 4, len(headers) + 1, 20, [24, 12, 34, 12, 12, 50])

    wb.save(XLSX)

    # ---------- docx: captions ----------
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(11)
    h = doc.add_heading("GeneReg 415 对无泄漏重分析：新增补充表格及表注", level=1)
    doc.add_paragraph(
        "本文档为《Supplementary_Tables_new_415_analysis.xlsx》中各表的表注与说明。全部表格对应修订稿方法 2.2 与结果 3.2 的重新分析："
        "冻结 370 对基准数据集与 45 对全新校准扩展记录（与 370 对零重叠）合并为 415 对唯一 TF–gene 对（200 阳/215 阴），"
        "以固定随机种子分层划分为开发集 332 对（160 阳/172 阴）与独立测试集 83 对（40 阳/43 阴）；"
        "权重与阈值仅在开发集上确定，测试集仅参与一次最终评估。")

    def table_doc(title, paras):
        doc.add_heading(title, level=2)
        for p in paras:
            doc.add_paragraph(p)

    table_doc("表 S2. 全新 45 对校准扩展数据集（calib45ext_v1_20260910）", [
        "表注：本表替代已退役的原表 S2 中 45 对校准记录部分。经核查，原 45 对校准记录与 370 对基准数据集存在完全的记录 ID 重叠，"
        "故按与基准数据集冻结协议完全一致的证据标准另行构建 45 对全新记录（22 阳、23 阴），与 370 对零重叠，"
        "合并后 415 对均为唯一 TF–gene 对。逐条入选与筛查审计见服务器目录 /www/wwwroot/gene_reg/validation/calibration_extension_45_20260910/。",
        "阳性入选标准：GTRD 启动子靶基因（≥2 个 meta-cluster）、ChIP-Atlas 平均靶基因分 > 0、且 ≥2 个不在本地峰值语料库中的非零 ChIP-Atlas SRX 实验；"
        "按与冻结版一致的排序规则取各转录因子冻结 top-20 之后的下一批候选。"
        "阴性入选标准（四条严格标准须同时满足）：该转录因子的 GTRD 启动子靶基因列表中缺失、ChIP-Atlas 靶基因评分 = 0、"
        "启动子区（TSS 上游 2000 bp 至下游 500 bp）JASPAR PWM 在 0.95 相对阈值下零命中、PubMed 标题/摘要精确共现次数为零；"
        "候选以新固定种子重排后逐一筛查，全部候选（含 28 条 JASPAR 拒绝、1 条来源错误）保留审计记录。",
        "列说明：Record ID 为扩展记录唯一标识（calib45ext 前缀）；Reference label 为参考标签（Positive/Negative）；"
        "Subset 区分阳性来源（calib45ext_GTRD_positive）与严格阴性（calib45ext_strict_negative）；"
        "GTRD site count 为 GTRD meta-cluster 数；ChIP-Atlas average score 与 independent non-zero SRX 为 ChIP-Atlas 靶基因评分及非本地非零 SRX 实验数；"
        "JASPAR PWM hits 为 0.95 相对阈值下启动子命中数（阴性均为 0）；PubMed co-occurrence count 为标题/摘要精确共现文献数（阴性均为 0）；"
        "后三列（ChIP-seq/Motif/Literature score, 0–100）为冻结 v2.6 确定性管线对该记录的行级证据评分，用于后续权重分析。",
    ])
    table_doc("表 S3. 开发集 332 对上全部 231 种候选权重组合的五折交叉验证结果", [
        "表注：本表替代已退役的原表 S2 中权重网格诊断部分。候选权重为 ChIP-seq、motif 与 literature 以 0.05 为步长、和为 1 的全部 231 种组合。"
        "五折分层交叉验证（固定种子）：每折阈值仅在训练四折上以训练折 min-max 归一化 + 0.00–1.00（步长 0.01）网格按 F1 最大选定，再原样应用于验证折。",
        "最终权重由事先声明的两步规则确定：① 入选资格要求每项权重 ≥ 0.05（网格最小非零步长），保证三种证据均参与综合评分（171/231 合格）；"
        "② 以全部组合中的全局最高平均 F1（0.00/0.60/0.40，0.9225）减一个标准误（SE = 0.0178）为界（0.9048），在界内合格组合中选取最小权重最大"
        "（三证据贡献最均衡）的组合（1-SE 规则），共 11 个组合落入界内，其中 0.10/0.65/0.25 与 0.10/0.70/0.20 峰值完全并列（开发集峰值 F1 均 0.9075、"
        "混淆矩阵相同、五折 F1 逐折相同）。并列由阈值稳健性判据（全阈值网格 0.00–1.00 的平均 F1，仅开发集）裁决：0.600 对 0.570，最终选定 0.10/0.70/0.20。",
        "列说明：Fold 1–5 F1 为各验证折 F1 值；Mean CV F1 / specificity / sensitivity 为五折均值；Eligible 列标记是否满足每项权重 ≥ 0.05；"
        "Within 1-SE 列标记是否位于全局最优的 1 个标准误界内；Selected 列标记最终选定组合。原稿权重 0.50/0.25/0.25 平均 F1 为 0.756，在合格组合中排名第 126/171。",
    ])
    table_doc("表 S4. 415 对合并数据集：分层划分、证据评分与最终预测（逐行）", [
        "表注：本表为修订后主分析的数据与预测明细，扩展并替代原表 S3 的记录级验证表。Partition 列标记分层随机划分结果"
        "（development 332 对 / test 83 对，固定随机种子，正负比例与整体一致）；开发集用于权重筛选与阈值优化，测试集在整个参数选择过程中保持未见。",
        "评分列说明：ChIP-seq/Motif/Literature score（0–100）为冻结 v2.6 确定性管线的行级证据评分；Weighted score 为选定权重（0.10/0.70/0.20）下的"
        "加权综合评分；Weighted score (min-max normalised) 为按开发集 min = 24.50、max = 67.00 归一化后的评分（测试集沿用开发集参数）；"
        "Final prediction 为归一化评分 ≥ 0.41 的最终预测（1 = 阳性）。Reference 列为冻结 v2.6 发布规则（raw_score ≥ 37）的原始评分与预测，仅作对照，未参与任何调参。",
    ])
    table_doc("表 S5. 独立测试集 83 对的单次最终性能评估", [
        "表注：本表为最终模型在独立测试集上的一次性评估结果（主表），对应原表 S3 的指标部分。最终模型：权重 ChIP-seq/motif/literature = 0.10/0.70/0.20，"
        "归一化阈值 0.41（对应原始加权评分约 41.93），归一化参数与阈值仅由开发集确定；其与 0.10/0.65/0.25 峰值并列，由阈值稳健性判据（全阈值网格平均 F1 0.600 对 0.570，"
        "仅开发集）裁决选定。所有指标的 95% 置信区间均由 1000 次非参数 bootstrap（百分位法，固定种子）估计。混淆矩阵：TP = 38、TN = 36、FP = 7、FN = 2。",
        "主要结果：F1 = 0.894（0.818–0.955），精确率 = 0.844，召回率（灵敏度）= 0.950，特异度 = 0.837，准确率 = 0.892（各指标 CI 见表内）。"
        "与并列组合 0.10/0.65/0.25 在同一测试集的结果（F1 = 0.907，95% CI 0.835–0.961）差异无统计学意义（CI 大幅重叠），印证两组合等效、最终差异仅由稳健性判据裁决。",
        "对照规则：同表并列报告冻结 v2.6 发布规则（raw_score ≥ 37，未重新调参）在同一测试集上的结果"
        "（F1 = 0.933，灵敏度 0.875，特异度 1.000，准确率 0.940；TP = 35、TN = 43、FP = 0、FN = 5），供读者比较；该规则与本分析的参数选择流程相互独立。",
    ])
    doc.add_paragraph(
        "注：以上四表的数据文件、分析脚本与随机种子明细见 E:\\9.2_manuscript\\GeneReg_415独立重跑_1SE均衡权重_20260910\\"
        "（merged_415_unique_split.csv、cv_weights_all231_with_eligibility.csv、metrics_summary.json、scripts/）；"
        "两配置并列裁决与胜者测试评估见同批 Two_config_comparison.xlsx 及 two_config_threshold_comparison.json。"
        "45 对扩展记录的构建与评分原件见同目录 server_materials/ 及服务器 /www/wwwroot/gene_reg/validation/calibration_extension_45_20260910/。"
        "表 S2、S3 为替代性新表（原表已随 45 对校准集退役），S4、S5 为新增；如期刊要求数字编号不同，可整体顺延。"
        "本版（2026-09-11）为最终版：最终权重 0.10/0.70/0.20、阈值 0.41，取代此前 0.10/0.65/0.25 版本；"
        "对应的 12 张单图 PDF 与图注见 9.11\\图片\\ 及《图注与表注_9.11.docx》，早期 R1–R4 合并面板图已作废。")
    doc.save(DOCX)
    print("written:", XLSX)
    print("written:", DOCX)


if __name__ == "__main__":
    main()
