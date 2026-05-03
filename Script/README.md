# 🧬 PLZ-HF: Reproducible Analysis Platform (Plasticizers × Heart Failure)

This repository provides a reproducible analysis platform for the PLZ–HF project: integrating plasticizer target profiles with subtype-resolved heart failure (HF) transcriptomic signatures (DCM, ICM, CTL) to derive mixture-relevant molecular axes, enrichment programs, and network-level hypotheses.

The codebase is organized for traceability, reuse, and figure-level reproducibility, with task-scoped scripts, explicit parameter definitions (see manuscript Table 2), and robustness checks (permutation and noise-injection).

---

## 🧲 What this repo contains

### 🫀 HF transcriptomic module (RNA-Seq)
- Multi-cohort RNA-Seq preprocessing and harmonization  
- Discovery-stage one-vs-rest Random Forest (RF) modeling for HF subtype classification  
- Feature stabilization via recurrence (nFreq), with TopNumber and nFreq sensitivity scans  
- Internal evaluation (ROC/PR, calibration, metrics) using out-of-fold predictions  
- External validation (GSE141910, GSE55296) and benchmarking (RF/XGB/SVM)  
- Cohort-sensitivity checks (LOCO)  
- Robustness checks (permutation negative control; noise-injection sensitivity)

### 🧪 Plasticizer module (PLZ)
- Toxicity/regulatory annotation summaries and exposure-tiering logic  
- Target-profile similarity clustering and cluster-number sensitivity

### 🧩 Integration module (PLZ × HF)
- Cluster-level target unions intersected with HF subtype DEG sets  
- Enrichment analysis (GO/KEGG/Reactome, MSigDB ORA)  
- Network construction and hub prioritization (STRING/Cytoscape exports)  
- Sankey-ready tables for axis-level summaries

---

## 🚀 Repository entry points

If you only run one pathway, start here:

`organized_by_task/Script/`

Scripts are grouped into three modules:
- `HF/` — HF discovery, tuning, robustness, validation, and evaluation  
- `PLZ/` — plasticizer-only analyses (annotation, clustering, sensitivity)  
- `PLZ_HF_integration/` — intersection, enrichment, network, and Sankey outputs  

A compact script inventory is available at:  
`organized_by_task/CODE_DIRECTORY.md`

---

## 🗂️ Directory overview

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
