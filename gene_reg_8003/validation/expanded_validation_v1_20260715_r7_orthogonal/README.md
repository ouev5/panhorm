# Expanded GeneReg independent validation set

This directory is immutable input/output for the expanded validation set. It is kept outside the active GeneReg scoring corpus (`encode_peaks/`).

- `raw/`: downloaded GTRD, ChIP-Atlas, JASPAR, Ensembl, and PubMed source data.
- `derived/combined_test_set_178_positive_192_negative.csv`: frozen joint test set.
- `derived/negative_screening_audit.csv`: all examined negative candidates and rejection reasons.
- `derived/source_manifest.json`: retrieval URL, timestamp, byte count and SHA-256 for every downloaded source file.
- `derived/local_chipatlas_accessions_excluded.json`: SRX accessions already in GeneReg's peak corpus and excluded from the GTRD positive confirmation.

The current GeneReg application is rule-based and does not read this directory; the test set is therefore held out from its scoring corpus.
