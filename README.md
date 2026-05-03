# PLZ-HF
Plasticizer-Responsive Molecular Axes in Heart Failure

# 🧬 PLZ–HF Reproducibility Hub

This repository provides a reproducible analysis workflow for integrating plasticizer target profiles with subtype-resolved HF (DCM, ICM, CTL) RNA-Seq signatures. It is organized for traceability, parameter transparency, and rerunnable figure-level outputs.

---

## 🎯 What you can reproduce here

- ✅ Multi-cohort HF RNA-Seq preprocessing and harmonization  
- ✅ One-vs-rest Random Forest subtype modeling with stability-based feature selection (TopNumber + nFreq)  
- ✅ Internal evaluation (ROC/PR, calibration) using out-of-fold predictions  
- ✅ External validation (GSE141910, GSE55296) and benchmark comparisons  
- ✅ Robustness checks (permutation negative control; noise-injection sensitivity)  
- ✅ Plasticizer clustering and PLZ×HF intersections, enrichment, and network export tables

---

## 🚀 Quick start

### 1) Environment
R (>= 4.5) with core packages:
`edgeR`, `limma`, `sva`, `caret`, `randomForest`, `pROC`, `clusterProfiler`

### 2) Recommended entry point
Run task-organized scripts from:

`organized_by_task/Script/`

Suggested high-level order:

1. `HF/01_prepare_raw_counts_and_metadata.R`
2. `HF/02_discovery_pipeline_main.R`
3. `HF/06_nfreq_frequency_scan.R`
4. `HF/08_internal_evaluation_plots_metrics.R`
5. `HF/09_validate_GSE141910.R`
6. `HF/10_validate_GSE55296.R`
7. `HF/12_permutation_negative_control.R`
8. `HF/13_noise_injection_sensitivity.R`

Plasticizer-only and integration scripts are under:
- `PLZ/`
- `PLZ_HF_integration/`

---

## 🧾 Key parameters (quick reference)

A compact parameter checklist aligned with the Methods is provided in:

`docs/PLZ_HF_key_parameters.md`

This includes: cohort definitions, nested CV design, TopNumber/nFreq settings, evaluation rules (OOF predictions, calibration bins), CI method, and robustness checks.

---

## 📦 Example output structure (minimal)

Example directory layout produced by the scripts:

```text
output/
├── HF/
│   ├── models/                 # fold-level models, selected features
│   ├── evaluation/             # ROC/PR/calibration summaries
│   ├── robustness/             # permutation + noise injection results
│   └── validation/             # external validation metrics and predictions
├── PLZ/
│   ├── annotations/            # toxicity/regulatory feature tables
│   └── clustering/             # Jaccard similarity + cluster sensitivity
└── PLZ_HF_integration/
    ├── intersections/          # overlap matrices and set tables
    ├── enrichment/             # GO/KEGG/Reactome/MSigDB ORA tables
    └── network_exports/        # Cytoscape/STRING-ready tables
