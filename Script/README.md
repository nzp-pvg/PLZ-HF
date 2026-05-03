# PLZ-HF Submission Code Package

This directory contains the code package used for the revised submission of the PLZ-HF manuscript.

## Purpose

The package is organized to support:

- RNA-seq preprocessing and discovery-stage one-vs-rest Random Forest modeling for heart failure (HF) subtype classification
- parameter scanning and robustness analyses
- external validation and benchmarking
- plasticizer-related downstream analyses
- HF-plasticizer integrative analyses for pathway, overlap, and hub-gene interpretation

The code organization in this submission package prioritizes readability, task-level separation, and reproducibility.

## Recommended code entry point

The recommended script set for review and reuse is:

`organized_by_task/Script`

This directory contains the task-organized scripts used for the revised submission. Scripts are grouped into:

- `HF`: HF discovery, tuning, robustness, and validation analyses
- `PLZ`: plasticizer-specific downstream analyses
- `PLZ_HF_integration`: HF-plasticizer integration analyses

## Directory overview

```text
code
├── README.md
├── organized_by_task
│   ├── CODE_DIRECTORY.md
│   └── Script
│       ├── HF
│       ├── PLZ
│       └── PLZ_HF_integration
├── PLZ-HF-main
│   ├── README.md
│   ├── LICENSE
│   ├── HF
│   ├── PLZ
│   └── Script
├── 01_preprocess_raw_counts.R
├── 02_discovery_ovr_rf_pipeline.R
├── 03_parameter_scans_and_robustness.R
├── 04_external_validation_and_evaluation.R
└── 05_plasticizer_downstream_analysis.R
```

## Main task-organized scripts

### HF

Located in:

`organized_by_task/Script/HF`

Primary scripts:

- `01_prepare_raw_counts_and_metadata.R`: raw RNA-seq count import, cohort harmonization, and count-matrix assembly
- `02_discovery_pipeline_main.R`: main discovery-stage one-vs-rest Random Forest pipeline
- `03_one_vs_each_training_function.R`: reusable one-vs-rest training function
- `04_topnumber_scan.R`: TopNumber sensitivity scan
- `05_topnumber_refine_500_600.R`: refined TopNumber scan around the final selected range
- `06_nfreq_frequency_scan.R`: nFreq recurrence scan across repeated CV outputs
- `07_nfreq_balance_selection.R`: subtype-specific nFreq balance-point selection
- `08_internal_evaluation_plots_metrics.R`: internal PR, calibration, ROC, and metric summaries
- `09_validate_GSE141910.R`: external validation on GSE141910
- `10_validate_GSE55296.R`: external validation and benchmarking on GSE55296
- `11_loco_cohort_sensitivity.R`: leave-one-cohort-out sensitivity analysis
- `12_permutation_negative_control.R`: permutation-based negative control
- `13_noise_injection_sensitivity.R`: Gaussian noise-injection sensitivity analysis
- `14_plot_fig_S3A_permutation.R`: publication plotting for the permutation panel

### PLZ

Located in:

`organized_by_task/Script/PLZ`

Primary scripts:

- `01_exposure_and_cluster_downstream.R`: cleaned downstream plasticizer analyses
- `02_exposure_evidence_tiering_37_original_logic.R`: evidence-tiering logic for plasticizer prioritization
- `03_cluster_k_by_k_compare_original_logic.R`: cluster-number sensitivity comparison

### PLZ_HF_integration

Located in:

`organized_by_task/Script/PLZ_HF_integration`

Primary scripts:

- `01_legacy_geo_helper_GSE116959.R`: retained legacy helper script; not part of the main integration workflow
- `02_prepare_hf_deg_sets_and_hgnc.R`: HF DEG preparation and HGNC symbol mapping
- `03_build_plasticizer_clusters_and_hf_intersections.R`: plasticizer cluster unions and HF-plasticizer overlap tables
- `04_dcm_cluster_enrichment_and_sankey.R`: DCM-focused enrichment and Sankey-ready outputs
- `05_icm_cluster_enrichment_and_sankey.R`: ICM-focused enrichment and Sankey-ready outputs
- `06_hub_ppi_and_cluster_pathways.R`: hub-gene, PPI, cluster-pathway, and Cytoscape export analyses
- `07_upset_hf_and_plasticizer_sets.R`: UpSet-style overlap summaries
- `08_hf_specific_kegg_reactome_enrichment.R`: HF-specific KEGG and Reactome enrichment
- `09_msigdb_enrichment_and_cytoscape_exports.R`: MSigDB enrichment and Cytoscape-friendly export tables
- `10_hub_core_logfc_summary.R`: hub-core follow-up summaries and export tables

## Suggested execution order

For HF discovery and validation:

1. `HF/01_prepare_raw_counts_and_metadata.R`
2. `HF/02_discovery_pipeline_main.R`
3. `HF/04_topnumber_scan.R`
4. `HF/05_topnumber_refine_500_600.R`
5. `HF/06_nfreq_frequency_scan.R`
6. `HF/07_nfreq_balance_selection.R`
7. `HF/08_internal_evaluation_plots_metrics.R`
8. `HF/09_validate_GSE141910.R`
9. `HF/10_validate_GSE55296.R`
10. `HF/11_loco_cohort_sensitivity.R`
11. `HF/12_permutation_negative_control.R`
12. `HF/13_noise_injection_sensitivity.R`
13. `HF/14_plot_fig_S3A_permutation.R`

For plasticizer-only analyses:

1. `PLZ/01_exposure_and_cluster_downstream.R`
2. `PLZ/02_exposure_evidence_tiering_37_original_logic.R`
3. `PLZ/03_cluster_k_by_k_compare_original_logic.R`

For HF-plasticizer integration:

1. `PLZ_HF_integration/02_prepare_hf_deg_sets_and_hgnc.R`
2. `PLZ_HF_integration/03_build_plasticizer_clusters_and_hf_intersections.R`
3. `PLZ_HF_integration/04_dcm_cluster_enrichment_and_sankey.R`
4. `PLZ_HF_integration/05_icm_cluster_enrichment_and_sankey.R`
5. `PLZ_HF_integration/06_hub_ppi_and_cluster_pathways.R`
6. `PLZ_HF_integration/07_upset_hf_and_plasticizer_sets.R`
7. `PLZ_HF_integration/08_hf_specific_kegg_reactome_enrichment.R`
8. `PLZ_HF_integration/09_msigdb_enrichment_and_cytoscape_exports.R`
9. `PLZ_HF_integration/10_hub_core_logfc_summary.R`

## Mapping to manuscript figures and supplementary figures

The task-organized scripts include comments indicating the main manuscript figure or table they support where applicable. Key mappings include:

- `HF/04_topnumber_scan.R` → Fig. S1A
- `HF/06_nfreq_frequency_scan.R`, `HF/07_nfreq_balance_selection.R` → Fig. S1B
- `HF/08_internal_evaluation_plots_metrics.R` → Fig. S1D-H
- `HF/10_validate_GSE55296.R` → Fig. S1F-H
- `HF/11_loco_cohort_sensitivity.R` → Fig. S1I
- `HF/12_permutation_negative_control.R` and `HF/14_plot_fig_S3A_permutation.R` → Fig. S3A
- `HF/13_noise_injection_sensitivity.R` → Fig. S3B

The manuscript figure copy used for revision tracking is:

`ENVINT-D-25-07776_R2.pdf`

## Reproducibility notes

- The organized scripts were cleaned to remove personal absolute local filesystem paths.
- The submission package retains both the original repository-style code (`PLZ-HF-main`) and the task-organized submission version (`organized_by_task/Script`).
- The task-organized scripts are intended to be the primary reference set for review, parameter tracing, and manuscript figure support.

## Additional script index

For a compact script inventory, see:

`organized_by_task/CODE_DIRECTORY.md`

