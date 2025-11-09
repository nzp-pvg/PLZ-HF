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

setwd("~/Science/CVD/CVD_MS_1/data/RNA-Seq/subtype/top560/validation/")

# Load discovery matrix (count_annotated, sample_info_discovery)

# Step 1: 加载完整数据
load("HF_processed_counts.RData")  # 包含 count, samples

table(samples$batch)

rownames(count) <- count$ensembl_gene_id
expr <- count[, -(1:3)]
rownames(samples) <- samples$sample_id
expr <- expr[, samples$sample_id]

#选取GSE141910矩阵, batch = B10, 没有ICM!!
validation_batches <- c("B10")
validation_samples <- samples$sample_id[
  samples$batch %in% validation_batches & samples$group %in% c("DCM", "ICM", "CTL")
]

meta_val <- samples[validation_samples, ]
expr_val <- expr[, validation_samples]


# 构建类标签表
classes_val <- meta_val[, c("sample_id", "group", "batch")]
colnames(classes_val) <- c("ID", "Classes", "batch")

# 加载建模结果
deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

# 验证函数框架（对每个 repeat 和 subclass 预测）
validate_models_on_valset <- function(model_result, expr_val, classes_val, subclass) {
  all_preds <- list()
  
  for (repeat_name in names(model_result)) {
    for (fold_name in names(model_result[[repeat_name]][[subclass]])) {
      res <- model_result[[repeat_name]][[subclass]][[fold_name]]
      model <- res$Model
      genes <- res$DEGs
      
      # 若验证集缺失某些基因，则跳过（保守策略）
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

# 合并两类模型的所有预测
dcm_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "DCM")
ctl_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "CTL")
val_preds <- bind_rows(dcm_pred, ctl_pred)


# ROC分析
library(pROC)
library(ggplot2)
library(dplyr)
library(purrr)

# 非平滑处理
dcm_only <- dcm_pred  # 确保这个是 validate_models_on_valset() 的输出
library(pROC)
roc_dcm <- roc(dcm_only$ActualClass == "DCM", dcm_only$One)
plot(roc_dcm, main = "DCM vs Others (Validation B10)")
auc(roc_dcm)
ci.auc(roc_dcm)

roc_ctl <- roc(ctl_pred$ActualClass == "CTL", ctl_pred$One)
plot(roc_ctl, main = "CTL vs Others (Validation B10)")
auc(roc_ctl)
ci.auc(roc_ctl)


设定颜色
subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0")

# 1. 构建 ROC 对象
roc_dcm <- roc(dcm_pred$ActualClass == "DCM", dcm_pred$One, quiet = TRUE)
roc_ctl <- roc(ctl_pred$ActualClass == "CTL", ctl_pred$One, quiet = TRUE)

# 2. 构建绘图数据框（非平滑）
roc_df <- rbind(
  data.frame(Spec = 1 - roc_dcm$specificities, Sens = roc_dcm$sensitivities, Class = "DCM"),
  data.frame(Spec = 1 - roc_ctl$specificities, Sens = roc_ctl$sensitivities, Class = "CTL")
)
roc_df$Class <- factor(roc_df$Class, levels = c("DCM", "CTL"))

# 3. AUC 和 CI
auc_ci_list <- list(
  DCM = ci.auc(roc_dcm),
  CTL = ci.auc(roc_ctl)
)

# 4. 绘图
ggplot(roc_df, aes(x = Spec, y = Sens, color = Class)) +
  # 先画 CTL 曲线
  geom_line(data = subset(roc_df, Class == "CTL"),
            aes(x = Spec, y = Sens, color = Class), size = 3) +
  # 再画 DCM 曲线
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



# 平滑处理
# Step 1: 准备子类列表（仅包含验证集中存在的类别）
subclasses_present <- c("DCM", "CTL")

# Step 2: 构建用于 ROC 分析的数据框（分别来自 dcm_pred 和 ctl_pred）
dcm_pred$Subclass <- "DCM"
ctl_pred$Subclass <- "CTL"
val_df <- bind_rows(dcm_pred, ctl_pred)

# Step 3: 生成 smoothed ROC + AUC + 95% CI
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

# Step 4: 整理为 data.frame 供 ggplot 使用
roc_df <- purrr::map2_df(roc_list_smooth, names(roc_list_smooth), function(roc, label) {
  coords <- coords(roc, ret = c("specificity", "sensitivity"))
  data.frame(Spec = 1 - coords$specificity, Sens = coords$sensitivity, Class = label)
})

# Step 5: AUC 标签构建
auc_label <- paste0(
  "DCM: ", sprintf("%.3f", auc_ci_list$DCM[2]), " (", 
  sprintf("%.3f", auc_ci_list$DCM[1]), "–", sprintf("%.3f", auc_ci_list$DCM[3]), ")\n",
  "CTL: ", sprintf("%.3f", auc_ci_list$CTL[2]), " (", 
  sprintf("%.3f", auc_ci_list$CTL[1]), "–", sprintf("%.3f", auc_ci_list$CTL[3]), ")"
)

# Step 6: 绘图
subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0")

roc_df$Class <- factor(roc_df$Class, levels = c("DCM", "CTL"))

ggplot(roc_df, aes(x = Spec, y = Sens, color = Class)) +
  # 先画 CTL 曲线
  geom_line(data = subset(roc_df, Class == "CTL"),
            aes(x = Spec, y = Sens, color = Class), size = 3) +
  
  # 再画 DCM 曲线（在上层）
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
  # 单独添加 DCM 文字，蓝色
  annotate("text", 
           x = 0.2, y = 0.1, 
           label = paste0("DCM: ", sprintf("%.3f", auc_ci_list$DCM[2]),
                          "(", sprintf("%.3f", auc_ci_list$DCM[1]), "–", sprintf("%.3f", auc_ci_list$DCM[3]), ")"),
           fontface = "italic", color = subclass_colors["DCM"],
           size = 18, hjust = 0.01) +
  # 单独添加 CTL 文字，灰蓝色
  annotate("text", 
           x = 0.23, y = 0.04, 
           label = paste0("CTL: ", sprintf("%.3f", auc_ci_list$CTL[2]),
                          "(", sprintf("%.3f", auc_ci_list$CTL[1]), "–", sprintf("%.3f", auc_ci_list$CTL[3]), ")"),
           fontface = "italic", color = subclass_colors["CTL"],
           size = 18, hjust = 0.01)


# 混淆矩阵
library(caret)
library(dplyr)

# 1. 合并两个子模型预测
dcm_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "DCM")
ctl_pred <- validate_models_on_valset(deg_result, expr_val, classes_val, "CTL")

library(dplyr)

# 给每个子模型预测结果加上 subclass-specific 列名
dcm_pred <- dcm_pred %>%
  rename(DCM = One) %>%
  select(-Others)

ctl_pred <- ctl_pred %>%
  rename(CTL = One) %>%
  select(-Others)



val_preds <- bind_rows(dcm_pred, ctl_pred)

# 将每个样本在每个 repeat-fold 中的概率得分汇总成 wide 格式
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

# 筛选只含 DCM 和 CTL 的预测结果
subset_preds <- final_preds %>%
  filter(ActualClass %in% c("DCM", "CTL"),
         PredictedClass %in% c("DCM", "CTL")) %>%
  mutate(
    ActualClass = factor(ActualClass, levels = c("DCM", "CTL")),
    PredictedClass = factor(PredictedClass, levels = c("DCM", "CTL"))
  )

# 生成2x2混淆矩阵
conf_matrix <- table(Prediction = subset_preds$PredictedClass,
                     Reference = subset_preds$ActualClass)

# 查看混淆矩阵内容
print(conf_matrix)

# 正类为 DCM
cm_dcm <- confusionMatrix(subset_preds$PredictedClass, subset_preds$ActualClass, positive = "DCM")

# 正类为 CTL
cm_ctl <- confusionMatrix(subset_preds$PredictedClass, subset_preds$ActualClass, positive = "CTL")

# 提取核心指标
summary_df <- data.frame(
  Class = c("DCM", "CTL"),
  Sensitivity = c(cm_dcm$byClass["Sensitivity"], cm_ctl$byClass["Sensitivity"]),
  Specificity = c(cm_dcm$byClass["Specificity"], cm_ctl$byClass["Specificity"]),
  Pos_Pred_Value = c(cm_dcm$byClass["Pos Pred Value"], cm_ctl$byClass["Pos Pred Value"]),
  Balanced_Accuracy = c(cm_dcm$byClass["Balanced Accuracy"], cm_ctl$byClass["Balanced Accuracy"]),
  Accuracy = c(cm_dcm$overall["Accuracy"], cm_ctl$overall["Accuracy"])
)

summary_df

# Step 1: 生成混淆矩阵
conf_matrix <- table(Prediction = subset_preds$PredictedClass,
                     Reference = subset_preds$ActualClass)

# 查看矩阵结构确认
print(conf_matrix)

# Step 2: 提取预测计数
TP_DCM <- conf_matrix["DCM", "DCM"]
FP_DCM <- conf_matrix["DCM", "CTL"]
FN_DCM <- conf_matrix["CTL", "DCM"]
TN_DCM <- conf_matrix["CTL", "CTL"]

TP_CTL <- conf_matrix["CTL", "CTL"]
FP_CTL <- conf_matrix["CTL", "DCM"]
FN_CTL <- conf_matrix["DCM", "CTL"]
TN_CTL <- conf_matrix["DCM", "DCM"]

# Step 3: 手动计算每类为正类时的准确率（Accuracy）
total <- sum(conf_matrix)

accuracy_DCM <- (TP_DCM + TN_DCM) / total
accuracy_CTL <- (TP_CTL + TN_CTL) / total

# Step 4: 输出结果
cat("DCM为正类时:\n")
cat("TP:", TP_DCM, "FP:", FP_DCM, "FN:", FN_DCM, "TN:", TN_DCM, "\n")
cat("Accuracy:", round(accuracy_DCM, 4), "\n\n")

cat("CTL为正类时:\n")
cat("TP:", TP_CTL, "FP:", FP_CTL, "FN:", FN_CTL, "TN:", TN_CTL, "\n")
cat("Accuracy:", round(accuracy_CTL, 4), "\n")

# PCA分析
# === HF_subtyping_v3.R ===
rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)


setwd("~/Science/CVD/CVD_MS_1/data/RNA-Seq/subtype/top560/validation/")


# ===============================
# 0. 加载必要包
# ===============================
library(edgeR)
library(FactoMineR)
library(factoextra)
library(dplyr)

# ===============================
# 1. 加载表达矩阵与样本信息
# ===============================
load("HF_processed_counts.RData")  # 提供 count, samples
rownames(count) <- count$ensembl_gene_id
expr <- count[, -(1:3)]
rownames(samples) <- samples$sample_id
expr <- expr[, samples$sample_id]

# ===============================
# 2. 提取验证队列 B10 中的 DCM 和 CTL 样本
# ===============================
validation_samples <- samples %>%
  filter(batch == "B10" & group %in% c("DCM", "CTL"))

classes_val <- validation_samples[, c("sample_id", "group", "batch")]
colnames(classes_val) <- c("ID", "Classes", "Batch")
expr_val <- expr[, classes_val$ID]

# ===============================
# 3. TMM 标准化 + logCPM
# ===============================
dge <- DGEList(counts = expr_val)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)

# ===============================
# 4. 提取建模中用到的 DEG 特征
# ===============================
Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")

# 根据你优化的频次阈值进行筛选
DEG_DCM <- names(Feature_list$DCM[Feature_list$DCM >= 100])
DEG_CTL <- names(Feature_list$CTL[Feature_list$CTL >= 100])
DEG_union <- union(DEG_DCM, DEG_CTL)

#DCM_specific <- setdiff(DEG_DCM, DEG_CTL)
#CTL_specific <- setdiff(DEG_CTL, DEG_DCM)
#pca_genes <- union(DCM_specific, CTL_specific)
#pca_genes <- pca_genes[pca_genes %in% rownames(expr_val)]  # 保证这些基因在表达矩阵中存在


# 提取 DEG 表达矩阵（行为基因，列为样本）
logCPM_pca <- logCPM[rownames(logCPM) %in% DEG_union, ]
#logCPM_pca <- logCPM[rownames(logCPM) %in% pca_genes, ]
# ===============================
# 5. 选择 DCM 和 CTL 样本并筛掉恒定表达基因
# ===============================
sel_samples <- classes_val$Classes %in% c("DCM", "CTL")
logCPM_sel <- logCPM_pca[, sel_samples]
labels <- classes_val$Classes[sel_samples]

# 去除表达恒定的基因
logCPM_filtered <- logCPM_sel[rowSums(logCPM_sel != logCPM_sel[, 1]) > 0, ]


# 2. 转置表达矩阵为样本在行
expr_pca_df <- as.data.frame(t(logCPM_filtered))
expr_pca_df$SampleID <- rownames(expr_pca_df)

# 3. 合并类标签
expr_pca_df <- left_join(expr_pca_df, classes_val, by = c("SampleID" = "ID"))

# 4. 提取表达矩阵（纯 numeric），用于 PCA 分析
gene_cols <- setdiff(colnames(expr_pca_df), c("SampleID", "Classes", "Batch"))
expr_matrix <- expr_pca_df[, gene_cols]
expr_matrix <- as.data.frame(lapply(expr_matrix, as.numeric))
rownames(expr_matrix) <- expr_pca_df$SampleID

# ========================================
# 5. 执行 PCA 分析
# ========================================
res.pca <- PCA(expr_matrix, graph = FALSE)

# ========================================
# 6. 绘图（使用 convex hull 而不是 ellipse）
# ========================================
library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)

# 获取 PCA scores
pca_scores <- as.data.frame(res.pca$ind$coord)
pca_scores$Class <- expr_pca_df$Classes
pca_scores$SampleID <- expr_pca_df$SampleID

# 计算 convex hull 的函数
compute_hull <- function(df) {
  df[chull(df$Dim.1, df$Dim.2), ]
}

hulls <- pca_scores %>%
  group_by(Class) %>%
  group_split() %>%
  lapply(compute_hull) %>%
  bind_rows()

# 提取主成分方差解释比例
eig_vals <- res.pca$eig

# 取前两主成分的解释百分比
pc1_var <- round(eig_vals[1, 2], 1)
pc2_var <- round(eig_vals[2, 2], 1)

# 打印
cat("Dim1: ", pc1_var, "%\n")
cat("Dim2: ", pc2_var, "%\n")

# subclass 色彩设定
subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")

# 1. 计算 convex hull
hulls <- pca_scores %>%
  group_by(Class) %>%
  slice(chull(Dim.1, Dim.2))

# 2. 计算每组中心点
centers <- pca_scores %>%
  group_by(Class) %>%
  summarise(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2), .groups = "drop")

# 3. 提取解释比例（假设已经赋值）
x_lab <- paste0("UMAP1 (", round(pc1_var, 1), "%)")
y_lab <- paste0("UMAP2 (", round(pc2_var, 1), "%)")

# 4. 绘图
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
# 1. 加载 hub 基因表 (ICM + DCM)
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

setwd("~/Science/CVD/CVD_MS_1/data/RNA-Seq/subtype/top560/validation/")

# ===============================
# 1. 加载数据
# ===============================
library(edgeR)
library(dplyr)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

load("HF_processed_counts.RData")  # 包含 count, samples

# 设置行名
rownames(count) <- count$ensembl_gene_id
expr <- count[, -(1:3)]
rownames(samples) <- samples$sample_id
expr <- expr[, samples$sample_id]

# ===============================
# 2. 提取 GSE141910 验证集 (batch = B10, DCM vs CTL)
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

# TMM 标准化 + logCPM
dge <- DGEList(counts = expr_val)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)

# ===============================
# 3. 载入 hub 基因
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
# 4. 定义 logFC 计算函数
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
# 5. 计算验证集 logFC (DCM vs CTL)
# ===============================
lfc_val_dcm <- calc_logfc_val("DCM", logCPM, classes_val)

# ENSEMBL → SYMBOL 映射
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
# 6. 合并 hub 基因和 logFC
# ===============================
lfc_val_all <- hub_nodes %>%
  left_join(lfc_val_dcm, by = "Gene")

# ===============================
# 7. 保存结果
# ===============================
write.csv(lfc_val_all, "GSE141910_ICM_and_DCM_hubs_in_DCMvsCTL_logFC.csv", row.names = FALSE)
print(lfc_val_all)
