# Expanded GeneReg validation set: frozen protocol

**Release:** `expanded_validation_v1_20260715_r7_orthogonal`  
**Size:** 370 human TF–gene records: 178 positive and 192 negative.

## Purpose and isolation

The GeneReg application on port 8003 is a rule-based evidence integrator. Its
`app.py` reads local peak files under `encode_peaks/`, JASPAR motifs, and PubMed
queries; it contains no training loop and does not load the validation
directory. The validation labels are therefore held out from the application
rule base. The frozen release is retained separately under
`/www/wwwroot/gene_reg/validation/expanded_validation_v1_20260715_r7_orthogonal`.

Because AR and ESR1 directories already contain many `chipatlas_SRX*.bed`
files, GTRD positives are additionally protected against direct use of those
local ChIP-Atlas datasets: an eligible GTRD positive must have support from at
least two *non-local* ChIP-Atlas SRX accessions. The exact excluded SRX lists
are saved in `derived/local_chipatlas_accessions_excluded.json`.

## Positive labels (n = 178)

- **GTRD-labelled external set (n = 160):** 20 records for each of AR, ESR1,
  PGR, NR3C1, PPARG, HNF4A, RARA, and GATA3. Each is a GTRD promoter target
  (`[-1000,+100]` from TSS) with at least two GTRD meta-clusters and at least
  two non-zero ChIP-Atlas SRX experiments not stored in GeneReg's local peak
  corpus.
- **Recent manual literature set (n = 18):** 11 HNF4A–gene records from PMID
  38909044 (2024) and 7 AR–gene records from PMID 40596315 (2025). These were
  included only when the publication reported binding/occupancy together with
  expression, perturbation, transactivation, or functional follow-up. Neither
  source was used to derive the GeneReg rule base.

## Strict negative labels (n = 192)

There are 24 records for each of the eight TFs above. Every retained negative
is absent from the TF-specific GTRD promoter target list, has a ChIP-Atlas
target score of zero, has no promoter hit for the configured JASPAR CORE PWM at
a relative threshold of 0.95, and has zero PubMed Title/Abstract
co-occurrences for the exact TF/gene query. `negative_screening_audit.csv`
retains accepted records and every rejected candidate with its rejection reason.

## Reproducibility materials

- `derived/combined_test_set_178_positive_192_negative.csv`: final labels and
  row-level evidence fields.
- `derived/negative_screening_audit.csv`: candidate-level audit trail.
- `derived/source_manifest.json`: retrieval times, URLs, file sizes, and SHA-256
  digests for downloaded GTRD and ChIP-Atlas source files.
- `raw/`: retained GTRD, ChIP-Atlas, JASPAR, Ensembl promoter, and PubMed XML
  inputs.
- `scripts/build_expanded_genereg_validation.py`: frozen build program.

Suggested manuscript wording: “We expanded the held-out validation corpus with
an external GTRD-labelled set and manually curated 2023–2025 literature
associations. The final frozen set contains 178 positive and 192 strict
negative TF–gene regulatory pairs. Positive labels require multiple binding
evidence, whereas strict negatives require concordant absence of GTRD,
ChIP-Atlas, high-stringency JASPAR promoter-motif, and PubMed co-occurrence
evidence.”
