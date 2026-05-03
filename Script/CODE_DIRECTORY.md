# organized_by_task
├── 00_CODE_DIRECTORY.md
├── HF
│   ├── 01_prepare_raw_counts_and_metadata.R  # raw RNA-seq count import, harmonization, merged matrix assembly
│   ├── 02_discovery_pipeline_main.R  # main discovery-stage one-vs-rest Random Forest pipeline; intermediate outputs for Fig. S1
│   ├── 03_one_vs_each_training_function.R  # one-vs-rest training function used by discovery and sensitivity scripts
│   ├── 04_topnumber_scan.R  # TopNumber sensitivity scan; manuscript: Fig. S1A
│   ├── 05_topnumber_refine_500_600.R  # refined TopNumber scan around 500-600
│   ├── 06_nfreq_frequency_scan.R  # nFreq recurrence scan across repeated CV outputs; inputs for Fig. S1B
│   ├── 07_nfreq_balance_selection.R  # subtype-specific nFreq balance-point selection; manuscript: Fig. S1B
│   ├── 08_internal_evaluation_plots_metrics.R  # internal PR, calibration, ROC, and metric summaries; manuscript: Fig. S1D-H
│   ├── 09_validate_GSE141910.R  # external validation on GSE141910
│   ├── 10_validate_GSE55296.R  # external validation and benchmarking on GSE55296; manuscript: Fig. S1F-H
│   ├── 11_loco_cohort_sensitivity.R  # leave-one-cohort-out sensitivity analysis; manuscript: Fig. S1I
│   ├── 12_permutation_negative_control.R  # permutation-based negative control; manuscript: Fig. S3A
│   ├── 13_noise_injection_sensitivity.R  # Gaussian noise-injection sensitivity analysis; manuscript: Fig. S3B
│   └── 14_plot_fig_S3A_permutation.R  # publication plotting for the Fig. S3A permutation panel
├── PLZ
│   ├── 01_exposure_and_cluster_downstream.R  # cleaned downstream plasticizer analyses
│   ├── 02_exposure_evidence_tiering_37_original_logic.R  # original logic for plasticizer exposure-evidence tiering
│   └── 03_cluster_k_by_k_compare_original_logic.R  # original logic for plasticizer cluster-k sensitivity comparison
└── PLZ_HF_integration
    ├── 01_legacy_geo_helper_GSE116959.R  # legacy GEO helper script retained for completeness; not part of the main PLZ-HF workflow
    ├── 02_prepare_hf_deg_sets_and_hgnc.R  # build HF DEG sets from one-vs-rest results and convert retained genes from Ensembl IDs to HGNC symbols
    ├── 03_build_plasticizer_clusters_and_hf_intersections.R  # define plasticizer cluster gene unions and compute HF-versus-plasticizer overlap tables
    ├── 04_dcm_cluster_enrichment_and_sankey.R  # DCM-focused cluster enrichment analyses and Sankey-ready outputs
    ├── 05_icm_cluster_enrichment_and_sankey.R  # ICM-focused cluster enrichment analyses and Sankey-ready outputs
    ├── 06_hub_ppi_and_cluster_pathways.R  # hub-gene, PPI, cluster-pathway, and Cytoscape export analyses for DCM and ICM
    ├── 07_upset_hf_and_plasticizer_sets.R  # UpSet-style overlap summaries for HF DEG sets and plasticizer cluster gene sets
    ├── 08_hf_specific_kegg_reactome_enrichment.R  # HF-specific KEGG and Reactome enrichment workflows
    ├── 09_msigdb_enrichment_and_cytoscape_exports.R  # MSigDB enrichment workflows and Cytoscape-friendly export tables
    └── 10_hub_core_logfc_summary.R  # hub-core follow-up analysis, per-batch logFC summarization, cutoff scanning, and export tables
