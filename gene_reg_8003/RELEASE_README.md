# AHormoneDB GeneReg v2.6

GeneReg v2.6 uses ChIP-seq / motif / literature weights **0.10/0.70/0.20** and a development-frozen normalized threshold of **0.41**. Normalization uses min **24.50**, max **67.00**; the corresponding raw weighted threshold is **41.925**. Optional AI explanations do not alter scores or labels.

The published 415-pair corpus is split into 332 development pairs and 83 internal test pairs. Development-only five-fold cross-validation and the threshold-robustness tie-break select the final weights. The internal-test confusion matrix is **38/36/7/2 (TP/TN/FP/FN)**, F1 **0.894** (95% CI **0.818–0.955**), AUC **0.923**, and AUPRC **0.924**. The evidence-absent reference set is not experimentally validated biological negatives.

See [README_v2.6.md](./README_v2.6.md) for methods, configuration, API usage, and offline reproduction; [release_manifest.json](./release_manifest.json) lists artifact checksums. The application entry points `app.py` and `app_v2.6.py` use the same model. Large peak directories and production secrets are excluded.
