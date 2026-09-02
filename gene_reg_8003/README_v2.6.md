# GeneReg v2.6

This directory contains the manuscript-matched GeneReg v2.6 release used for the AHormoneDB revision.

- Deterministic scorer: 12-feature logistic regression with ChIP-Atlas integration.
- Fixed release threshold: raw score `>= 37.0`.
- Optional LLM explanations do not modify the deterministic score.
- Validation is reported as a frozen internal benchmark, not as independent external generalization.
- Frozen benchmark: 370 pairs (178 positive, 192 strict negative).
- Release ID: `AHormoneDB-GeneReg-v2.6-20260902`.

See [`RELEASE_README_v2.6.md`](./RELEASE_README_v2.6.md), [`release_manifest_v2.6.json`](./release_manifest_v2.6.json), and the [`validation/`](./validation/) directory for protocols, checksums, row-level predictions, and audit materials.

Large runtime resources such as ENCODE peak directories and secrets are intentionally excluded from this public repository.
