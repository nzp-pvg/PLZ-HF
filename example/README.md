# Example inputs

This directory contains lightweight example inputs for trial runs and format inspection.

## Files

- `toy_counts_matrix.tsv`: minimal toy count matrix
- `toy_metadata.tsv`: minimal toy sample metadata
- `HF_RNASeq_demo/`: compact demo subset extracted from real intermediate HF RNA-seq objects

## HF_RNASeq_demo contents

The `HF_RNASeq_demo/` directory contains small subsets derived from the Zenodo-style intermediate objects used in the full analysis.

### Discovery-stage demo objects

- `Discovery_start_data_demo.RData`
  - objects: `classes_df_demo`, `expr_mat_demo`
- `Discovery_start_data_demo_canonical_names.RData`
  - objects: `classes_df`, `expr_mat`
  - dimensions: `800 genes × 36 samples`
  - class balance: `12 CTL`, `12 DCM`, `12 ICM`

### Processed-count demo objects

- `HF_processed_counts_demo.RData`
  - objects: `count_keep`, `samples_keep`
- `HF_processed_counts_demo_canonical_names.RData`
  - objects: `count`, `samples`
  - dimensions: `800 genes × 36 samples` plus annotation columns

### Merged processed-count demo objects

- `GSE_processed_counts_demo.RData`
  - objects: `exp_keep`, `sample_info_keep`
- `GSE_processed_counts_demo_canonical_names.RData`
  - objects: `exp_no_0`, `sample_info`

### Per-cohort matrix examples

- `GSE46224_matrix_demo.rds`
- `GSE141910_matrix_demo.rds`
- `GSE55296_matrix_demo.rds`

These retain the original per-cohort matrix layout but only contain a small subset of rows and gene columns.

### Manifest

- `MANIFEST.tsv`: file-level summary of demo object purposes

## Intended use

These example files are provided so that users can:

- inspect object structure
- test read-in and path logic
- run small-scale dry runs without exposing manuscript-critical result tables

They are not intended to reproduce the final manuscript results.
