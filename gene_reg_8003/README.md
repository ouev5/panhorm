# GeneReg v2.6 — TF target gene prediction

GeneReg combines **ChIP-seq (0.10)**, **JASPAR 2024 motif (0.70)**, and **PubMed literature (0.20)** evidence. The normalized decision threshold is **0.41**, with development-frozen min–max parameters **24.50/67.00** (raw weighted threshold **41.925**). Optional AI explanations contribute zero points.

The manuscript benchmark contains **415 unique pairs**: **332 development** and **83 internal test**. Internal-test results are **F1 0.894** (95% CI 0.818–0.955), **AUC 0.923**, **AUPRC 0.924**, with **TP/TN/FP/FN = 38/36/7/2**. Operational evidence-absent labels are not experimentally validated biological negatives.

See **[the complete methods, configuration, API and reproduction guide](./README_v2.6.md)**.

```bash
python -m pip install -r requirements-analysis.txt
python analysis_v2.6/scripts/reproduce_manuscript.py --cross-validate
python -m pytest tests -q
python app.py  # serves port 8003; restore runtime evidence resources first
```

The [model configuration](./model_config.json), [portable 415-pair data](./analysis_v2.6/manuscript_benchmark_415.csv), [recomputed metrics](./analysis_v2.6/manuscript_metrics.json), and [release manifest](./release_manifest.json) provide the reproducibility materials.
