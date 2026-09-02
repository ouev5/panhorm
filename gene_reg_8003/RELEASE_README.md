# AHormoneDB GeneReg v2.6 frozen release

- Release ID: `AHormoneDB-GeneReg-v2.6-20260902`
- Production service: `http://127.0.0.1:8003` (public host should be `https://genereg.panhorm.cn`)
- Deterministic scorer: 12-feature logistic-regression model with fixed raw-score threshold 37.
- Optional DeepSeek explanation: disabled for validation and contributes zero points to the score.
- The legacy three-evidence 0.50/0.25/0.25 grid is retained only as a diagnostic; it is not the deployed v2.6 model and its reconstructed F1 is 0.8148 (rank 147).
- Frozen benchmark: 370 pairs (178 positive, 192 strict negative), test-set SHA-256 `2f0e8cb737d5e91eda6393bf0432e059e290a41b6069af2bc192a5ee8b02fac2`.
- Frozen application SHA-256: `c0f98646487737319bf5d29c3d8b1eb6560ae6c8a3ed71a933174baa69735a30`.
- Metrics are reported as internal-benchmark diagnostics; no claim of external generalization is made because an independently documented development/test threshold-selection split is not included in this release.

The release package contains the scoring source, validation protocol, negative-screening audit, source manifest, input checksums, row-level predictions, metric summary, and the PTPN1 per-edge outputs used to revise the manuscript. Production secrets and large runtime peak directories are intentionally excluded.
