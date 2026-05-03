# === HF_subtyping_v3.R ===
rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)
library(limma)
library(caret)
library(randomForest)
library(glmnet)

# Load discovery matrix (count_annotated, sample_info_discovery)

load("HF_processed_counts.RData")

table(samples$batch)

rownames(count) <- count$ensembl_gene_id
expr <- count[, -(1:3)]
rownames(samples) <- samples$sample_id
expr <- expr[, samples$sample_id]

validation_batches <- c("B10")
validation_samples <- samples$sample_id[
  samples$batch %in% validation_batches & samples$group %in% c("DCM", "ICM", "CTL")
]

meta_val <- samples[validation_samples, ]
expr_val <- expr[, validation_samples]


classes_val <- meta_val[, c("sample_id", "group", "batch")]
colnames(classes_val) <- c("ID", "Classes", "batch")

deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

validate_models_on_valset <- function(model_result, expr_val, classes_val, subclass) {
  all_preds <- list()
  
  for (repeat_name in names(model_result)) {
    for (fold_name in names(model_result[[repeat_name]][[subclass]])) {
      res <- model_result[[repeat_name]][[subclass]][[fold_name]]
      model <- res$Model
      genes <- res$DEGs
      
      if (!all(genes %in% rownames(expr_val))) next
      
      expr_input <- t(expr_val[genes, , drop = FALSE])
      pred_probs <- predict(model, newdata = expr_input, type = "prob")
      
      df <- data.frame(pred_probs)
      df$Repeat <- repeat_name
      df$Fold <- fold_name
      df$Subclass <- subclass
      df$SampleID <- rownames(expr_input)
      df$ActualClass <- classes_val$Classes[match(df$SampleID, classes_val$ID)]
      df$PredictedClass <- predict(model, newdata = expr_input, type = "raw")
      
      all_preds[[paste0(repeat_name, "_", fold_name)]] <- df
    }
  }
  
  do.call(rbind, all_preds)
}

dcm_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "DCM")
ctl_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "CTL")
val_preds <- bind_rows(dcm_pred, ctl_pred)


library(pROC)
library(ggplot2)
library(dplyr)
library(purrr)

dcm_only <- dcm_pred
library(pROC)
roc_dcm <- roc(dcm_only$ActualClass == "DCM", dcm_only$One)
plot(roc_dcm, main = "DCM vs Others (Validation B10)")
auc(roc_dcm)
ci.auc(roc_dcm)

roc_ctl <- roc(ctl_pred$ActualClass == "CTL", ctl_pred$One)
plot(roc_ctl, main = "CTL vs Others (Validation B10)")
auc(roc_ctl)
ci.auc(roc_ctl)


# Set colors
subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0")

roc_dcm <- roc(dcm_pred$ActualClass == "DCM", dcm_pred$One, quiet = TRUE)
roc_ctl <- roc(ctl_pred$ActualClass == "CTL", ctl_pred$One, quiet = TRUE)

roc_df <- rbind(
  data.frame(Spec = 1 - roc_dcm$specificities, Sens = roc_dcm$sensitivities, Class = "DCM"),
  data.frame(Spec = 1 - roc_ctl$specificities, Sens = roc_ctl$sensitivities, Class = "CTL")
)
roc_df$Class <- factor(roc_df$Class, levels = c("DCM", "CTL"))

auc_ci_list <- list(
  DCM = ci.auc(roc_dcm),
  CTL = ci.auc(roc_ctl)
)

ggplot(roc_df, aes(x = Spec, y = Sens, color = Class)) +
  geom_line(data = subset(roc_df, Class == "CTL"),
            aes(x = Spec, y = Sens, color = Class), size = 3) +
  geom_line(data = subset(roc_df, Class == "DCM"),
            aes(x = Spec, y = Sens, color = Class), size = 3) +
  scale_color_manual(values = subclass_colors) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80" , size = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 50),
    axis.title = element_text(color = "black", size = 50)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves (Validation-1, Raw)") 
  annotate("text", 
           x = 0.2, y = 0.1, 
           label = paste0("DCM: ", sprintf("%.3f", auc_ci_list$DCM[2]),
                          " (", sprintf("%.3f", auc_ci_list$DCM[1]), "–", sprintf("%.3f", auc_ci_list$DCM[3]), ")"),
           fontface = "italic", color = subclass_colors["DCM"], size = 18, hjust = 0.01) +
  annotate("text", 
           x = 0.23, y = 0.04, 
           label = paste0("CTL: ", sprintf("%.3f", auc_ci_list$CTL[2]),
                          " (", sprintf("%.3f", auc_ci_list$CTL[1]), "–", sprintf("%.3f", auc_ci_list$CTL[3]), ")"),
           fontface = "italic", color = subclass_colors["CTL"], size = 18, hjust = 0.01)



subclasses_present <- c("DCM", "CTL")

dcm_pred$Subclass <- "DCM"
ctl_pred$Subclass <- "CTL"
val_df <- bind_rows(dcm_pred, ctl_pred)

roc_list_smooth <- list()
auc_ci_list <- list()

for (sc in subclasses_present) {
  df <- val_df %>% filter(Subclass == sc)
  roc_obj <- roc(
    response = factor(ifelse(df$ActualClass == sc, sc, "Other")),
    predictor = df$One,
    ci = TRUE,
    smooth = TRUE,
    quiet = TRUE
  )
  roc_list_smooth[[sc]] <- roc_obj
  auc_ci_list[[sc]] <- ci.auc(roc_obj)
}

roc_df <- purrr::map2_df(roc_list_smooth, names(roc_list_smooth), function(roc, label) {
  coords <- coords(roc, ret = c("specificity", "sensitivity"))
  data.frame(Spec = 1 - coords$specificity, Sens = coords$sensitivity, Class = label)
})

auc_label <- paste0(
  "DCM: ", sprintf("%.3f", auc_ci_list$DCM[2]), " (", 
  sprintf("%.3f", auc_ci_list$DCM[1]), "–", sprintf("%.3f", auc_ci_list$DCM[3]), ")\n",
  "CTL: ", sprintf("%.3f", auc_ci_list$CTL[2]), " (", 
  sprintf("%.3f", auc_ci_list$CTL[1]), "–", sprintf("%.3f", auc_ci_list$CTL[3]), ")"
)

subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0")

roc_df$Class <- factor(roc_df$Class, levels = c("DCM", "CTL"))

ggplot(roc_df, aes(x = Spec, y = Sens, color = Class)) +
  geom_line(data = subset(roc_df, Class == "CTL"),
            aes(x = Spec, y = Sens, color = Class), size = 3) +
  
  geom_line(data = subset(roc_df, Class == "DCM"),
            aes(x = Spec, y = Sens, color = Class), size = 3) +
  scale_color_manual(values = subclass_colors) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80" , size = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 50),
    axis.title = element_text(color = "black", size = 50)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Smoothed ROC Curves_Validation-1") 
  annotate("text", 
           x = 0.2, y = 0.1, 
           label = paste0("DCM: ", sprintf("%.3f", auc_ci_list$DCM[2]),
                          "(", sprintf("%.3f", auc_ci_list$DCM[1]), "–", sprintf("%.3f", auc_ci_list$DCM[3]), ")"),
           fontface = "italic", color = subclass_colors["DCM"],
           size = 18, hjust = 0.01) +
  annotate("text", 
           x = 0.23, y = 0.04, 
           label = paste0("CTL: ", sprintf("%.3f", auc_ci_list$CTL[2]),
                          "(", sprintf("%.3f", auc_ci_list$CTL[1]), "–", sprintf("%.3f", auc_ci_list$CTL[3]), ")"),
           fontface = "italic", color = subclass_colors["CTL"],
           size = 18, hjust = 0.01)


library(caret)
library(dplyr)

dcm_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "DCM")
ctl_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "CTL")

library(dplyr)

dcm_pred <- dcm_pred %>%
  rename(DCM = One) %>%
  select(-Others)

ctl_pred <- ctl_pred %>%
  rename(CTL = One) %>%
  select(-Others)



val_preds <- bind_rows(dcm_pred, ctl_pred)

final_preds <- val_preds %>%
  group_by(SampleID, Repeat, Fold) %>%
  summarise(
    DCM = mean(DCM, na.rm = TRUE),
    CTL = mean(CTL, na.rm = TRUE),
    ActualClass = dplyr::first(ActualClass),
    .groups = "drop"
  ) %>%
  mutate(
    PredictedClass = colnames(across(DCM:CTL))[apply(across(DCM:CTL), 1, which.max)]
  )

subset_preds <- final_preds %>%
  filter(ActualClass %in% c("DCM", "CTL"),
         PredictedClass %in% c("DCM", "CTL")) %>%
  mutate(
    ActualClass = factor(ActualClass, levels = c("DCM", "CTL")),
    PredictedClass = factor(PredictedClass, levels = c("DCM", "CTL"))
  )

conf_matrix <- table(Prediction = subset_preds$PredictedClass,
                     Reference = subset_preds$ActualClass)

print(conf_matrix)

cm_dcm <- confusionMatrix(subset_preds$PredictedClass, subset_preds$ActualClass, positive = "DCM")

cm_ctl <- confusionMatrix(subset_preds$PredictedClass, subset_preds$ActualClass, positive = "CTL")

summary_df <- data.frame(
  Class = c("DCM", "CTL"),
  Sensitivity = c(cm_dcm$byClass["Sensitivity"], cm_ctl$byClass["Sensitivity"]),
  Specificity = c(cm_dcm$byClass["Specificity"], cm_ctl$byClass["Specificity"]),
  Pos_Pred_Value = c(cm_dcm$byClass["Pos Pred Value"], cm_ctl$byClass["Pos Pred Value"]),
  Balanced_Accuracy = c(cm_dcm$byClass["Balanced Accuracy"], cm_ctl$byClass["Balanced Accuracy"]),
  Accuracy = c(cm_dcm$overall["Accuracy"], cm_ctl$overall["Accuracy"])
)

summary_df

conf_matrix <- table(Prediction = subset_preds$PredictedClass,
                     Reference = subset_preds$ActualClass)

print(conf_matrix)

TP_DCM <- conf_matrix["DCM", "DCM"]
FP_DCM <- conf_matrix["DCM", "CTL"]
FN_DCM <- conf_matrix["CTL", "DCM"]
TN_DCM <- conf_matrix["CTL", "CTL"]

TP_CTL <- conf_matrix["CTL", "CTL"]
FP_CTL <- conf_matrix["CTL", "DCM"]
FN_CTL <- conf_matrix["DCM", "CTL"]
TN_CTL <- conf_matrix["DCM", "DCM"]

total <- sum(conf_matrix)

accuracy_DCM <- (TP_DCM + TN_DCM) / total
accuracy_CTL <- (TP_CTL + TN_CTL) / total

cat("Using DCM as the positive class:\n")
cat("TP:", TP_DCM, "FP:", FP_DCM, "FN:", FN_DCM, "TN:", TN_DCM, "\n")
cat("Accuracy:", round(accuracy_DCM, 4), "\n\n")

cat("Using CTL as the positive class:\n")
cat("TP:", TP_CTL, "FP:", FP_CTL, "FN:", FN_CTL, "TN:", TN_CTL, "\n")
cat("Accuracy:", round(accuracy_CTL, 4), "\n")

# === HF_subtyping_v3.R ===
rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)

# ===============================
# ===============================
library(edgeR)
library(FactoMineR)
library(factoextra)
library(dplyr)

# ===============================
# ===============================
load("HF_processed_counts.RData")
rownames(count) <- count$ensembl_gene_id
expr <- count[, -(1:3)]
rownames(samples) <- samples$sample_id
expr <- expr[, samples$sample_id]

# ===============================
# ===============================
validation_samples <- samples %>%
  filter(batch == "B10" & group %in% c("DCM", "CTL"))

classes_val <- validation_samples[, c("sample_id", "group", "batch")]
colnames(classes_val) <- c("ID", "Classes", "Batch")
expr_val <- expr[, classes_val$ID]

# ===============================
# ===============================
dge <- DGEList(counts = expr_val)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)

# ===============================
# ===============================
Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")

DEG_DCM <- names(Feature_list$DCM[Feature_list$DCM >= 100])
DEG_CTL <- names(Feature_list$CTL[Feature_list$CTL >= 100])
DEG_union <- union(DEG_DCM, DEG_CTL)

#DCM_specific <- setdiff(DEG_DCM, DEG_CTL)
#CTL_specific <- setdiff(DEG_CTL, DEG_DCM)
#pca_genes <- union(DCM_specific, CTL_specific)


logCPM_pca <- logCPM[rownames(logCPM) %in% DEG_union, ]
#logCPM_pca <- logCPM[rownames(logCPM) %in% pca_genes, ]
# ===============================
# ===============================
sel_samples <- classes_val$Classes %in% c("DCM", "CTL")
logCPM_sel <- logCPM_pca[, sel_samples]
labels <- classes_val$Classes[sel_samples]

logCPM_filtered <- logCPM_sel[rowSums(logCPM_sel != logCPM_sel[, 1]) > 0, ]


expr_pca_df <- as.data.frame(t(logCPM_filtered))
expr_pca_df$SampleID <- rownames(expr_pca_df)

expr_pca_df <- left_join(expr_pca_df, classes_val, by = c("SampleID" = "ID"))

gene_cols <- setdiff(colnames(expr_pca_df), c("SampleID", "Classes", "Batch"))
expr_matrix <- expr_pca_df[, gene_cols]
expr_matrix <- as.data.frame(lapply(expr_matrix, as.numeric))
rownames(expr_matrix) <- expr_pca_df$SampleID

# ========================================
# ========================================
res.pca <- PCA(expr_matrix, graph = FALSE)

# ========================================
# ========================================
library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)

pca_scores <- as.data.frame(res.pca$ind$coord)
pca_scores$Class <- expr_pca_df$Classes
pca_scores$SampleID <- expr_pca_df$SampleID

compute_hull <- function(df) {
  df[chull(df$Dim.1, df$Dim.2), ]
}

hulls <- pca_scores %>%
  group_by(Class) %>%
  group_split() %>%
  lapply(compute_hull) %>%
  bind_rows()

eig_vals <- res.pca$eig

pc1_var <- round(eig_vals[1, 2], 1)
pc2_var <- round(eig_vals[2, 2], 1)

cat("Dim1: ", pc1_var, "%\n")
cat("Dim2: ", pc2_var, "%\n")

subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")

hulls <- pca_scores %>%
  group_by(Class) %>%
  slice(chull(Dim.1, Dim.2))

centers <- pca_scores %>%
  group_by(Class) %>%
  summarise(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2), .groups = "drop")

x_lab <- paste0("UMAP1 (", round(pc1_var, 1), "%)")
y_lab <- paste0("UMAP2 (", round(pc2_var, 1), "%)")

ggplot(pca_scores, aes(x = Dim.1, y = Dim.2, fill = Class)) +
  stat_density_2d(aes(fill = Class), geom = "polygon", contour = TRUE, alpha = 0.15, color = NA) +
  geom_point(shape = 21, size = 7, stroke = 0.5, alpha = 0.8, color = "black") +
  geom_polygon(data = hulls, aes(group = Class), color = "black", linewidth = 0, alpha = 0.1, show.legend = FALSE) +
  geom_point(data = centers, shape = 21, size = 18, alpha = 0.5, stroke = 1.2, color = "black") +
  scale_fill_manual(values = subclass_colors) +
  labs(
    x = x_lab,
    y = y_lab,
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    panel.grid = element_blank(),
    axis.text = element_text(size = 46, color = "black"),
    axis.title = element_text(size = 46, color = "black"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = -10, b = 0, l = 0)),
    legend.title = element_blank(),
    legend.text = element_text(size = 46),
    plot.title = element_text(size = 48, face = "bold", hjust = 0.5)
  )






# ===============================
# ===============================
rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)
library(limma)


# ===============================
# ===============================
library(edgeR)
library(dplyr)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

load("HF_processed_counts.RData")

rownames(count) <- count$ensembl_gene_id
expr <- count[, -(1:3)]
rownames(samples) <- samples$sample_id
expr <- expr[, samples$sample_id]

# ===============================
# ===============================
validation_samples <- samples %>%
  filter(batch == "B10" & group %in% c("DCM", "CTL"))

classes_val <- data.frame(
  ID = validation_samples$sample_id,
  Classes = validation_samples$group,
  Batch   = validation_samples$batch,
  stringsAsFactors = FALSE
)

expr_val <- expr[, classes_val$ID]

dge <- DGEList(counts = expr_val)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)

# ===============================
# ===============================
icm_tbl <- read_csv("ICM_hub_top60_edges_nodes.csv", show_col_types = FALSE) %>% as.data.frame()
dcm_tbl <- read_csv("DCM_hub_top60_edges_nodes.csv", show_col_types = FALSE) %>% as.data.frame()

icm_hubs <- icm_tbl %>%
  dplyr::arrange(desc(Degree)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(Gene, Cluster, Degree) %>%
  dplyr::mutate(Disease = "ICM")

dcm_hubs <- dcm_tbl %>%
  dplyr::arrange(desc(Degree)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(Gene, Cluster, Degree) %>%
  dplyr::mutate(Disease = "DCM")

hub_nodes <- bind_rows(icm_hubs, dcm_hubs) %>%
  group_by(Gene, Disease) %>%
  summarise(
    Cluster = paste(unique(Cluster), collapse="|"),
    Degree = max(Degree),
    .groups = "drop"
  )

# ===============================
# ===============================
calc_logfc_val <- function(subclass, expr, classes){
  case_ids <- classes %>% filter(Classes == subclass) %>% pull(ID)
  ctrl_ids <- classes %>% filter(Classes == "CTL") %>% pull(ID)
  
  if(length(case_ids) > 1 && length(ctrl_ids) > 1){
    case_mean <- rowMeans(expr[, case_ids, drop = FALSE])
    ctrl_mean <- rowMeans(expr[, ctrl_ids, drop = FALSE])
    logFC <- case_mean - ctrl_mean
    return(data.frame(ENSEMBL = rownames(expr), Val_DCM_logFC = logFC))
  } else {
    return(NULL)
  }
}

# ===============================
# ===============================
lfc_val_dcm <- calc_logfc_val("DCM", logCPM, classes_val)

map_val <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = lfc_val_dcm$ENSEMBL,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL","SYMBOL")
) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

lfc_val_dcm <- lfc_val_dcm %>%
  dplyr::left_join(map_val, by = "ENSEMBL") %>%
  dplyr::select(Gene = SYMBOL, Val_DCM_logFC)

# ===============================
# ===============================
lfc_val_all <- hub_nodes %>%
  left_join(lfc_val_dcm, by = "Gene")

# ===============================
# ===============================
write.csv(lfc_val_all, "GSE141910_ICM_and_DCM_hubs_in_DCMvsCTL_logFC.csv", row.names = FALSE)
print(lfc_val_all)
