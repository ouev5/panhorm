# GeneReg v2.6 — weighted TF–target gene prediction

GeneReg integrates ChIP-seq binding evidence, JASPAR 2024 motif scanning, and exact PubMed title/abstract co-occurrence (2010–2025) to prioritize transcription factor–target gene relationships, including nuclear receptor TFs. It does not predict ligand–receptor binding or establish causal disease mechanisms.

## Final scoring rule

```text
weighted_score = 0.10 × ChIP-seq + 0.70 × motif + 0.20 × literature
normalized_score = (weighted_score − 24.50) / (67.00 − 24.50)
prediction = normalized_score ≥ 0.41
```

All three component scores are on a 0–100 scale. The min–max parameters are fitted on the 332-pair development set and frozen for subsequent predictions. The normalized threshold corresponds to a raw weighted score of **41.925**. New observations may lie outside the development range; scores are not clipped before classification. Full-precision values determine the label. Optional DeepSeek explanations contribute **zero points**.

The manuscript's local ChIP-seq reference corpus comprises **436 TFs and 3380 BED files** from ENCODE and ChIP-Atlas. The application displays the actual available directory count in `/api/health`; large peak files must be restored separately. Promoter motif scanning uses the **−2000/+500 bp TSS window**, both strands, and a **0.95 relative PWM threshold**. ChIP component penalties, signal weighting, and cofactor evidence are retained within the ChIP score; no extra learned classifier or post-combination penalty modifies the weighted rule.

## Development and internal test benchmark

The 415-pair corpus contains 200 positives and 215 operational evidence-absent pairs, with a stratified fixed-seed split:

| Partition | Total | Positive | Evidence-absent |
|---|---:|---:|---:|
| Development | 332 | 160 | 172 |
| Internal test | 83 | 40 | 43 |

Development-only five-fold stratified cross-validation evaluated 231 weight combinations in 0.05 steps. All evidence weights must be at least 0.05, leaving 171 eligible combinations. The one-standard-error cutoff is 0.9048. The tied combinations 0.10/0.65/0.25 and 0.10/0.70/0.20 have mean CV F1 0.9083; the final combination is selected by development threshold-grid robustness (mean F1 **0.600** versus **0.570**). Min–max parameters and the final threshold are fitted on development data only.

| Internal-test metric | Estimate | Bootstrap 95% CI |
|---|---:|---|
| F1 | 0.894 | 0.818–0.955 |
| Precision | 0.844 | 0.735–0.939 |
| Recall | 0.950 | 0.871–1.000 |
| Accuracy | 0.892 | 0.819–0.952 |
| TN proportion on evidence-absent pairs | 0.837 | 0.721–0.935 |
| ROC AUC | 0.923 | — |
| PR AUC | 0.924 | — |

Confusion matrix: **TP 38, TN 36, FP 7, FN 2**. Confidence intervals use 1000 percentile bootstrap resamples. These are internal-benchmark results, not external independent validation. Evidence-absent labels mean no support was retrieved under the specified screening criteria; they are not experimentally validated biological negatives.

## Run and reproduce

```bash
cd gene_reg_8003
python -m pip install -r requirements.txt
# Optional AI explanation: provide DEEPSEEK_API_KEY in the process environment.
python app.py
```

Runtime resources are configured through environment variables; see [`.env.example`](./.env.example). Merely copying `.env.example` does not load it into Python: export the variables or use a process manager's environment-file option. No real key is needed for offline benchmark reproduction.

```bash
python -m pip install -r requirements-analysis.txt
python analysis_v2.6/scripts/reproduce_manuscript.py --cross-validate
python -m pytest tests -q
```

The reproduction script uses the **published frozen component scores**, not fresh external API responses: it verifies all 415 workbook predictions, the non-overlapping split, development normalization, selected weights, confusion matrix, and bootstrap intervals. Live requests depend on restored peak resources and upstream evidence availability.

### API

- `GET /api/analyze?gene=TFF1&receptor=ESR1`: deterministic prediction.
- `POST /api/analyze`: `{"gene":"TFF1","receptor":"ESR1","include_ai":false}`.
- `include_ai=true`: request a qualitative explanation without changing the score.
- `GET /api/health`: runtime data coverage and scoring configuration.
- `GET /api/validation`: published 415-pair corpus, development/test counts, hash, and internal-test results.

Responses include component scores, weights, weighted and normalized scores, the threshold, and the binary `prediction`. The web interface and CSV export display the same decision rule.

## Reproducibility files

- [`model_config.json`](./model_config.json): final frozen parameters.
- [`analysis_v2.6/manuscript_benchmark_415.csv`](./analysis_v2.6/manuscript_benchmark_415.csv): portable published row-level split and predictions.
- [`analysis_v2.6/manuscript_metrics.json`](./analysis_v2.6/manuscript_metrics.json): recomputed internal-test summary.
- [`analysis_v2.6/results/Supplementary_Tables_S2_S5.xlsx`](./analysis_v2.6/results/Supplementary_Tables_S2_S5.xlsx): published source tables.
- [`release_manifest.json`](./release_manifest.json): artifact hashes.
