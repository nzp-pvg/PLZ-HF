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

> Scope: method parameters and evaluation rules for reproducibility.
> ⚠️ Note: this document lists settings and definitions only (no primary-result claims).

---

## 🧬 Data and cohorts

- 📚 **Discovery cohorts (RNA-Seq, LV tissue)**: GSE46224, GSE116250, GSE120852, GSE135055, GSE48166, GSE126569  
- 🧪 **Independent validation cohorts**: GSE141910 (DCM vs CTL), GSE55296 (DCM/ICM/CTL)
- 📥 **Input type**: raw RNA-Seq **count** matrices

---

## 🧹 Preprocessing (RNA-Seq)

- 🧯 **Batch correction**: `sva::ComBat_seq()`  
  - batch variable: cohort-of-origin  
  - input scale: counts
- ⚖️ **Normalization (within-training)**: `edgeR::DGEList()` + `calcNormFactors(method="TMM")`
- 🔁 **Transformation / modeling**: `limma::voom()` (logCPM + precision weights)

---

## 🧬 DEG ranking (training-only)

- ✅ **Training-only rule**: DEG ranking and feature selection are performed strictly **within training splits only** (no information from held-out folds used for DEG ranking)
- 🧪 **Pipeline**: `voom → lmFit → contrasts.fit → eBayes(trend=TRUE)`
- 📊 **Ranking statistic**: limma **moderated t-statistic** (one-vs-rest contrast)

---

## 🤖 Classification task

- 🎯 **Scheme**: one-vs-rest Random Forest (RF)
  - DCM vs (ICM + CTL)
  - ICM vs (DCM + CTL)
  - CTL vs (DCM + ICM)

---

## 🔄 Nested resampling design

- 🧱 **Outer evaluation**: stratified 10-fold CV × 100 outer repeats  
  - total fold-level selections: **100 × 10 = 1000**
  - evaluation uses **out-of-fold** predictions from the outer CV
- 🧰 **Inner fitting/tuning (caret)**: `trainControl(method="repeatedcv", number=10, repeats=3, classProbs=TRUE, savePredictions="final")`
- 🧭 **caret optimization metric**: Kappa  
  - reporting focuses on AUROC/AUPRC and complementary metrics

---

## 🧬 Feature selection (TopNumber)

- 🧷 **TopNumber**: 560 DEGs per one-vs-rest task
- ↕️ **Symmetric rule**: select 280 most positive + 280 most negative DEGs  
  - corresponds to opposite directions in the one-vs-rest contrast

---

## 🧲 Stability filter (nFreq)

- 🔁 **Definition**: fold-level recurrence across **1000** fold-level TopNumber lists  
  - nFreq = number of fold-level TopNumber lists (out of 1000) in which a gene appears  
- 🎚️ **Subtype-specific thresholds**:
  - CTL = 2
  - DCM = 8
  - ICM = 17
- 🧠 **Selection rationale**: thresholds guided by sensitivity curves (trade-off between performance stability and signature compactness)

---

## 🎛️ RF hyperparameters

- ✅ **Tuned parameter**: `mtry ≈ TopNumber/3`  
- 🧱 **Inherited randomForest defaults** (unless otherwise stated):
  - `ntree = 500`
  - `nodesize = 1` (classification)
  - `replace = TRUE`

---

## ⚖️ Class imbalance handling

- ✅ **Primary RF pipeline**: stratified resampling (no explicit `classwt`)
- 🧪 **Sensitivity check (weighted RF)**: `randomForest::classwt` (inverse prevalence; neg/pos)

---

## 📈 Evaluation rules and metrics

- 🏆 **Primary metric**: AUROC  
- 📌 **Complementary metrics**: AUPRC, macro-F1, balanced accuracy, accuracy
- 🧾 **Prediction source**: out-of-fold predictions only (outer CV)
- 🧠 **Multiclass decision rule**: argmax of one-vs-rest predicted probabilities
- 🎯 **Calibration**:
  - 10 quantile-based bins of predicted probabilities
  - ECE and Brier score computed from binned summaries
- 📉 **AUROC 95% CI**: `pROC::ci.auc()` (DeLong)

---

## 🌍 External evaluation (primary RF transfer)

- 🚚 **Transfer rule (primary RF external validation)**:
  - apply discovery-derived DEGs and fold-specific trained models directly to external cohorts
  - no feature reselection; no retraining during the transfer step
- 🏷️ **Label constraint rule**: external evaluation follows label-defined tasks (no subtype imputation or label redefinition)

---

## 🧪 Robustness checks

- 🔀 **Permutation negative control**: 500 label permutations under identical pipeline
- 🌫️ **Noise injection sensitivity**: Gaussian noise 0%, 5%, 10%; 12 repeats per condition

---

## ♻️ Reproducibility

- 🎲 **Random seed**: `set.seed(2025)` (primary discovery-stage scripts)
- 🧰 **Core software**: R, edgeR, limma, sva, caret, randomForest, pROC (see manuscript Table 2 for versions)

---


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
