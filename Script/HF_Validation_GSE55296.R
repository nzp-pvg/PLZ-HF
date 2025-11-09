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



setwd("~/Science/CVD/CVD_MS_1/data/RNA-Seq/subtype/top560/validation/")

# 准备表达矩阵与样本信息
# 必要库
library(readxl)
library(biomaRt)
library(dplyr)

library(readxl)
library(biomaRt)

# Step 1: 读取样本信息
sample_info <- read.table("GSE55296_sample.txt", 
                          sep = "\t", 
                          header = TRUE, 
                          stringsAsFactors = FALSE)

# Step 2: 构建标签表 classes_val
classes_val <- data.frame(
  ID = sample_info$Sample_ID,
  Classes = sample_info$Group,
  batch = "validation",
  stringsAsFactors = FALSE
)

# Step 3: 读取表达矩阵（含 GeneID 列）
count_data <- read_excel("GSE55296_count_data.xlsx")
count_data <- as.data.frame(count_data)

# Step 4: 清洗 GeneID（确保无 NA、无重复）
count_data <- count_data[!is.na(count_data$GeneID), ]
count_data <- count_data[!duplicated(count_data$GeneID), ]

# Step 5: 设置行名为 GeneID，去除 GeneID 列
rownames(count_data) <- count_data$GeneID
count_data_mat <- count_data[, !(colnames(count_data) %in% "GeneID")]

# Step 6: 提取表达矩阵中 sample_info 对应的样本列
expr_val <- count_data_mat[, sample_info$Sample_ID]

# Step 7: 将验证集表达矩阵的 ENSEMBL ID（hg19）映射到 hg38 保留有效 ID
# 使用 hg19（GRCh37）biomaRt 查询当前仍在 hg38 中有效的 ENSEMBL ID
mart_hg19 <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)

gene_ids <- rownames(expr_val)

gene_map <- getBM(
  attributes = c("ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart_hg19
)

valid_genes <- unique(gene_map$ensembl_gene_id)
expr_val_hg38 <- expr_val[valid_genes, ]

# Step 8: 检查保留基因数量
cat("保留的验证集基因数：", nrow(expr_val_hg38), " / 原始基因数：", nrow(expr_val), "\n")

# 模型验证

# 加载发现队列构建的模型结果（按 One-vs-Eeach 策略）
deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

validate_models_on_valset <- function(model_result, expr_val, classes_val, subclass, min_gene_match = 0.8) {
  all_preds <- list()
  
  for (repeat_name in names(model_result)) {
    for (fold_name in names(model_result[[repeat_name]][[subclass]])) {
      res <- model_result[[repeat_name]][[subclass]][[fold_name]]
      model <- res$Model
      genes <- res$DEGs
      
      # 跳过无效模型或缺失基因
      if (is.null(model) || is.null(genes)) {
        message("Skipping ", repeat_name, "/", fold_name, "/", subclass, ": model or genes missing.")
        next
      }
      
      # 检查模型使用的变量
      model_vars <- names(model$forest$xlevels)
      if (is.null(model_vars) || length(model_vars) == 0) {
        message("Skipping ", repeat_name, "/", fold_name, "/", subclass, ": model_vars is NULL or empty.")
        next
      }
      
      # 匹配训练基因和验证集中的表达基因
      matched_genes <- intersect(genes, rownames(expr_val))
      match_ratio <- length(matched_genes) / length(genes)
      if (is.na(match_ratio) || match_ratio < min_gene_match) {
        message("Skipping ", repeat_name, "/", fold_name, "/", subclass, ": insufficient gene match (", round(match_ratio, 2), ")")
        next
      }
      
      # 构建表达矩阵输入
      expr_input <- t(expr_val[matched_genes, , drop = FALSE])
      colnames(expr_input) <- matched_genes
      
      # ✅ 容错处理：只选模型中实际用到的变量
      used_vars <- intersect(model_vars, colnames(expr_input))
      if (length(used_vars) == 0) {
        message("Skipping ", repeat_name, "/", fold_name, "/", subclass, ": no shared variables.")
        next
      }
      expr_input <- expr_input[, used_vars, drop = FALSE]
      
      # 预测
      pred_probs <- predict(model, newdata = expr_input, type = "prob")
      pred_labels <- predict(model, newdata = expr_input, type = "response")
      
      # 保存结果
      df <- data.frame(pred_probs)
      df$Repeat <- repeat_name
      df$Fold <- fold_name
      df$Subclass <- subclass
      df$SampleID <- rownames(expr_input)
      df$ActualClass <- classes_val$Classes[match(df$SampleID, classes_val$ID)]
      df$PredictedClass <- pred_labels
      
      all_preds[[paste0(repeat_name, "_", fold_name)]] <- df
    }
  }
  
  if (length(all_preds) == 0) {
    warning("No valid predictions generated for subclass ", subclass)
    return(NULL)
  }
  
  return(do.call(rbind, all_preds))
}

icm_pred <- validate_models_on_valset(deg_result, expr_val_hg38, classes_val, "ICM")
dcm_pred <- validate_models_on_valset(deg_result, expr_val_hg38, classes_val, "DCM")
ctl_pred <- validate_models_on_valset(deg_result, expr_val_hg38, classes_val, "CTL")

val_preds <- dplyr::bind_rows(icm_pred, dcm_pred, ctl_pred)



# 提取所有 DEGs（包含所有 repeat / subclass / fold）
all_deg_genes <- unique(unlist(lapply(deg_result, function(repeat_list) {
  lapply(repeat_list, function(subclass_list) {
    lapply(subclass_list, function(fold_list) {
      fold_list$DEGs
    })
  })
})))

# 检查与验证集表达矩阵中的基因是否重合
total_deg_genes <- length(all_deg_genes)
matched_genes <- sum(all_deg_genes %in% rownames(expr_val_hg38))

# 输出结果
cat("deg_result 中所有 DEGs 总数：", total_deg_genes, "\n")
cat("其中与验证集表达矩阵匹配上的基因数：", matched_genes, "\n")
cat("匹配比例：", round(matched_genes / total_deg_genes, 4) * 100, "%\n")

# PCA
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(dplyr)
library(umap)

# Step 1: logCPM 转换
log_expr <- log2(expr_val_hg38 + 1)

# Step 2: 提取 DEG
Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")
DEG_ICM <- names(Feature_list$DCM[Feature_list$ICM >= 22])
DEG_DCM <- names(Feature_list$DCM[Feature_list$DCM >= 7])
DEG_CTL <- names(Feature_list$CTL[Feature_list$CTL >= 70])

# （可选）Step 3: 进一步提取专属基因集
ICM_only <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM_only <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL_only <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))

three_class_DEGs <- unique(c(ICM_only, DCM_only, CTL_only))

# Step 3: 表达矩阵
logCPM_pca <- log_expr[rownames(log_expr) %in% three_class_DEGs, ]
sel_samples <- classes_val$Classes %in% c("ICM", "DCM", "CTL")
logCPM_sel <- logCPM_pca[, classes_val$ID[sel_samples]]
labels <- classes_val$Classes[sel_samples]

# Step 4: 去除恒定表达基因
logCPM_filtered <- logCPM_sel[rowSums(logCPM_sel != logCPM_sel[, 1]) > 0, ]

# Step 5: 转置并合并样本注释
expr_pca_df <- as.data.frame(t(logCPM_filtered))
expr_pca_df$SampleID <- rownames(expr_pca_df)
expr_pca_df <- left_join(expr_pca_df, classes_val, by = c("SampleID" = "ID"))

gene_cols <- setdiff(colnames(expr_pca_df), c("SampleID", "Classes", "batch"))
expr_matrix <- expr_pca_df[, gene_cols]
expr_matrix <- as.data.frame(lapply(expr_matrix, as.numeric))
rownames(expr_matrix) <- expr_pca_df$SampleID

# Step 6: 执行 UMAP
set.seed(2025)
umap_result <- umap(expr_matrix)
custom_config <- umap.defaults
custom_config$n_neighbors <- 10       # 默认15，越小越局部，强化类内聚类
custom_config$min_dist <- 0.2        # 默认0.1-0.5之间，越小越紧凑（增强类间间隔）
custom_config$spread <- 0.7         # 可微调类间散布程度
umap_result <- umap(expr_matrix, config = custom_config)

umap_plot_df <- data.frame(
  SampleID = rownames(expr_matrix),
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2]
) %>%
  left_join(expr_pca_df[, c("SampleID", "Classes")], by = "SampleID")

# Step 7: 中心点计算
centers <- umap_plot_df %>%
  group_by(Classes) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

# Step 8: 方差比例（近似）
umap_var <- apply(umap_result$layout, 2, var)
umap_var_exp <- round(100 * umap_var / sum(umap_var), 1)

# Step 9: 绘图
# 添加边界线（每类的凸包）
geom_polygon(data = umap_plot_df %>%
               group_by(Classes) %>%
               slice(chull(UMAP1, UMAP2)),
             aes(x = UMAP1, y = UMAP2, fill = Classes), 
             color = "black", linewidth = 0.1, alpha = 0, show.legend = FALSE)

ggplot(umap_plot_df, aes(x = UMAP1, y = UMAP2)) +
  stat_density_2d(aes(fill = Classes), geom = "polygon", alpha = 0.1, color = NA) +
  geom_point(aes(fill = Classes), shape = 21, size = 8, stroke = 1, alpha = 0.8, color = "black") +
  geom_point(data = centers, aes(x = UMAP1, y = UMAP2, fill = Classes),
             shape = 21, size = 18, stroke = 1, alpha = 0.8, color = "black") +
  geom_polygon(data = umap_plot_df %>%
                 group_by(Classes) %>%
                 slice(chull(UMAP1, UMAP2)),
               aes(x = UMAP1, y = UMAP2, fill = Classes),
               color = "black", linewidth = 0.5, alpha = 0.2, show.legend = FALSE) +

  scale_fill_manual(values = c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")) +
  labs(
    x = paste0("UMAP1 (", umap_var_exp[1], "%)"),
    y = paste0("UMAP2 (", umap_var_exp[2], "%)")
  ) +
  theme_minimal(base_size = 48) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    panel.grid = element_blank(),
    axis.text = element_text(size = 46, color = "black"),
    axis.title = element_text(size = 46, color = "black"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 0)),
    axis.title.y = element_text(margin = ggplot2::margin(r = -10)),
    legend.title = element_blank(),
    legend.text = element_text(size = 46),
    plot.title = element_text(size = 48, face = "bold", hjust = 0.5)
  )


ggplot(umap_plot_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(fill = logCPM_filtered["GENE_X", SampleID]), shape = 21) +
  scale_fill_gradient2() + ...

# 必要库
library(randomForest)
library(dplyr)
library(caret)
library(MLmetrics)

# Step 1: 读取表达矩阵（使用你之前处理好的 expr_val_hg38）
#         也确保 classes_val 已准备好

# Step 2: 读取 DEG 频次表
Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")

# Step 3: 筛选高频 DEG（例如频次 ≥ 100）
top_features <- names(Filter(function(x) x >= 100, Feature_list$DCM))
top_features <- intersect(top_features, rownames(expr_val_hg38))  # 仅保留在验证集中的

cat("最终用于验证的 DEG 数量：", length(top_features), "\n")

# Step 4: 构建表达矩阵（转置）
val_expr_input <- t(expr_val_hg38[top_features, ])
val_expr_input <- as.data.frame(val_expr_input)
val_expr_input$Class <- classes_val$Classes[match(rownames(val_expr_input), classes_val$ID)]
val_expr_input$Class <- as.factor(val_expr_input$Class)
# Step 5: 建立分类模型（使用交叉验证训练）
set.seed(2025)
ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = multiClassSummary)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,  # 因为你现在是三分类任务
  savePredictions = "all"
)

rf_model <- caret::train(
  Class ~ ., data = val_expr_input,
  method = "rf",
  trControl = ctrl,
  importance = TRUE
)

# Step 6: 输出模型性能
print(rf_model)


library(caret)
library(dplyr)

# 1. 初始化变量
set.seed(42)
n_repeat <- 100
rf_model_list <- list()

# 2. 训练多轮随机森林模型
for (i in 1:n_repeat) {
  cat("Training Repeat", i, "\n")
  
  set.seed(100 + i)
  
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 1,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = "final"
  )
  
  rf_model <- caret::train(
    Class ~ ., data = val_expr_input,
    method = "rf",
    trControl = ctrl,
    importance = TRUE
  )
  
  rf_model_list[[i]] <- rf_model
}

# 3. 提取每轮预测结果（使用最佳 mtry）
all_preds <- data.frame()

for (i in seq_along(rf_model_list)) {
  model <- rf_model_list[[i]]
  best_mtry <- model$bestTune$mtry
  
  preds <- model$pred %>%
    filter(mtry == best_mtry) %>%
    mutate(Repeat = paste0("Repeat", i))
  
  all_preds <- bind_rows(all_preds, preds)
}

# 4. 汇总混淆矩阵（真实数值）
conf_mat_total <- confusionMatrix(all_preds$pred, all_preds$obs)

cat("\n===== 总混淆矩阵（原始计数） =====\n")
print(conf_mat_total$table)

cat("\n===== 总体性能指标 =====\n")
print(conf_mat_total$overall)

# 可选：按类别精度汇总
cat("\n===== 每类指标 =====\n")
print(conf_mat_total$byClass)






library(caret)
library(dplyr)

# 参数设置
n_repeat <- 100
model_dir <- "RF_RepeatModels"
dir.create(model_dir, showWarnings = FALSE)

# 控制参数（每轮都是 10-fold CV）
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  savePredictions = "final"
)

# 主循环
rf_model_list <- list()
acc_table <- data.frame()

for (i in 1:n_repeat) {
  cat("==== Training Repeat", i, "====\n")
  set.seed(1000 + i)
  
  # 模型训练
  rf_model <- caret::train(
    Class ~ ., data = val_expr_input,
    method = "rf",
    trControl = ctrl,
    importance = TRUE
  )
  
  # 保存模型对象
  model_path <- file.path(model_dir, paste0("rf_model_repeat", i, ".rds"))
  saveRDS(rf_model, model_path)
  
  # 保存简要性能指标
  best_row <- rf_model$results[rf_model$results$mtry == rf_model$bestTune$mtry, ]
  acc_table <- rbind(acc_table, data.frame(
    Repeat = i,
    mtry = best_row$mtry,
    Accuracy = best_row$Accuracy,
    Kappa = best_row$Kappa
  ))
  
  # 保存进列表（可选）
  rf_model_list[[i]] <- rf_model
}

# 保存整体性能表
write.csv(acc_table, file = "RF_100Repeat_Accuracy_Kappa.csv", row.names = FALSE)
saveRDS(rf_model_list, file = "RF_100Repeat_ModelList.rds")

cat("\n✅ 已完成 100 轮训练，每轮模型和性能指标均已保存。\n")


library(caret)
library(dplyr)

# 1. 汇总所有重复的预测结果
all_preds <- lapply(rf_model_list, function(model) {
  model$pred
}) %>% bind_rows()

# 检查
dim(all_preds)
head(all_preds)

# 转为 factor，确保 levels 一致
all_preds$pred <- factor(all_preds$pred, levels = c("CTL", "DCM", "ICM"))
all_preds$obs <- factor(all_preds$obs, levels = c("CTL", "DCM", "ICM"))

# 混淆矩阵
conf_mat_total <- confusionMatrix(all_preds$pred, all_preds$obs)

cat("===== 总混淆矩阵（原始计数） =====\n")
print(conf_mat_total$table)

cat("\n===== 总体性能指标 =====\n")
print(conf_mat_total$overall)

cat("\n===== 每类指标 =====\n")
print(conf_mat_total$byClass)

library(pROC)
library(tidyr)
library(dplyr)

# 只保留需要列
pred_df <- all_preds[, c("obs", "CTL", "DCM", "ICM")]

# 转换为 long 格式
long_pred <- pred_df %>%
  pivot_longer(cols = -obs, names_to = "Class", values_to = "Prob") %>%
  mutate(Label = ifelse(obs == Class, 1, 0))

# 按 Class 计算 ROC 曲线
roc_list <- long_pred %>%
  group_by(Class) %>%
  group_split() %>%
  lapply(function(df) {
    roc(response = df$Label, predictor = df$Prob)
  })

# 输出 AUC + 95% CI
auc_ci_list <- lapply(roc_list, function(roc_obj) ci.auc(roc_obj))
# 手动为 ROC 对象命名
names(roc_list) <- unique(long_pred$Class)  # 确保 Class 是 CTL, DCM, ICM

# 计算每一类的 AUC 和 95% CI
auc_ci_df <- lapply(roc_list, function(r) {
  ci <- ci.auc(r)
  data.frame(
    AUC = as.numeric(ci[2]),
    CI_lower = as.numeric(ci[1]),
    CI_upper = as.numeric(ci[3])
  )
}) %>%
  bind_rows(.id = "Class")

# 输出清晰标注的表格
print(auc_ci_df)

library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)

# Step 1: 准备 long 格式数据
pred_df <- all_preds[, c("obs", "CTL", "DCM", "ICM")]

long_pred <- pred_df %>%
  pivot_longer(cols = -obs, names_to = "Class", values_to = "Prob") %>%
  mutate(Label = ifelse(obs == Class, 1, 0))

# Step 2: 计算 ROC 曲线和 AUC
roc_df_list <- list()
auc_ci_list <- list()

for (cls in unique(long_pred$Class)) {
  df <- filter(long_pred, Class == cls)
  roc_obj <- roc(response = df$Label, predictor = df$Prob)
  roc_df <- data.frame(
    Spec = 1 - roc_obj$specificities,
    Sens = roc_obj$sensitivities,
    Class = cls
  )
  roc_df_list[[cls]] <- roc_df
  auc_ci <- ci.auc(roc_obj)
  auc_ci_list[[cls]] <- c(lower = auc_ci[1], auc = auc_ci[2], upper = auc_ci[3])
}

roc_df <- bind_rows(roc_df_list)
auc_ci_list <- do.call(rbind, auc_ci_list)
rownames(auc_ci_list) <- c("CTL", "DCM", "ICM")

# Step 3: 设置颜色和顺序
subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")
roc_df$Class <- factor(roc_df$Class, levels = c("DCM", "CTL", "ICM"))

# Step 4: 绘图
ggplot(roc_df, aes(x = Spec, y = Sens, color = Class)) +
  geom_line(data = subset(roc_df, Class == "CTL"), size = 3) +
  geom_line(data = subset(roc_df, Class == "DCM"), size = 3) +
  geom_line(data = subset(roc_df, Class == "ICM"), size = 3) +
  scale_color_manual(values = subclass_colors) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80", size = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 50),
    axis.title = element_text(color = "black", size = 50)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Smoothed ROC Curves_Validation-2") 
  annotate("text", 
           x = 0.20, y = 0.12, 
           label = paste0("DCM: ", sprintf("%.3f", auc_ci_list["DCM", "auc"]),
                          " (", sprintf("%.3f", auc_ci_list["DCM", "lower"]),
                          "–", sprintf("%.3f", auc_ci_list["DCM", "upper"]), ")"),
           fontface = "italic", color = subclass_colors["DCM"],
           size = 18, hjust = 0.01) +
  annotate("text", 
           x = 0.23, y = 0.05, 
           label = paste0("CTL: ", sprintf("%.3f", auc_ci_list["CTL", "auc"]),
                          " (", sprintf("%.3f", auc_ci_list["CTL", "lower"]),
                          "–", sprintf("%.3f", auc_ci_list["CTL", "upper"]), ")"),
           fontface = "italic", color = subclass_colors["CTL"],
           size = 18, hjust = 0.01) +
  annotate("text", 
           x = 0.23, y = 0.18, 
           label = paste0("ICM: ", sprintf("%.3f", auc_ci_list["ICM", "auc"]),
                          " (", sprintf("%.3f", auc_ci_list["ICM", "lower"]),
                          "–", sprintf("%.3f", auc_ci_list["ICM", "upper"]), ")"),
           fontface = "italic", color = subclass_colors["ICM"],
           size = 18, hjust = 0.01)




## XGB algorithm
library(caret)
library(pROC)
library(dplyr)
library(tidyr)

# Step 0: 读取 DEG 频次表
Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")
# Step 1: 控制参数

set.seed(42)
n_repeat <- 100
xgb_model_list <- list()
all_preds_xgb <- data.frame()

# Step 2: 定义交叉验证参数
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 1,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  savePredictions = "final"
)

# Step 3: 开始重复训练
for (i in 1:n_repeat) {
  cat("Training Repeat", i, "\n")
  
  set.seed(100 + i)
  model <- caret::train(
    Class ~ ., data = val_expr_input,
    method = "xgbTree",
    trControl = ctrl,
    metric = "Accuracy"
  )
  
  xgb_model_list[[i]] <- model
  all_preds_xgb <- bind_rows(all_preds_xgb, model$pred)
}

# 总混淆矩阵
conf_mat_total_xgb <- confusionMatrix(
  factor(all_preds_xgb$pred, levels = c("CTL", "DCM", "ICM")),
  factor(all_preds_xgb$obs, levels = c("CTL", "DCM", "ICM"))
)

cat("===== XGBoost 总混淆矩阵（原始计数） =====\n")
print(conf_mat_total_xgb$table)

cat("\n===== XGBoost 总体性能指标 =====\n")
print(conf_mat_total_xgb$overall)

cat("\n===== XGBoost 每类指标 =====\n")
print(conf_mat_total_xgb$byClass)

saveRDS(all_preds_xgb, file = "GSE55296_all_preds_xgb.rds")

# ROC
library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)

# Step 1: 设定颜色
subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")

# Step 2: 预处理预测数据
pred_df_xgb <- all_preds_xgb[, c("obs", "CTL", "DCM", "ICM")]

long_pred_xgb <- pred_df_xgb %>%
  pivot_longer(cols = -obs, names_to = "Class", values_to = "Prob") %>%
  mutate(Label = ifelse(obs == Class, 1, 0))

# Step 3: 计算ROC对象
roc_list_xgb <- long_pred_xgb %>%
  group_by(Class) %>%
  group_split() %>%
  lapply(function(df) {
    roc_obj <- roc(response = df$Label, predictor = df$Prob, quiet = TRUE)
    return(roc_obj)
  })

names(roc_list_xgb) <- c("CTL", "DCM", "ICM")

# Step 4: 提取AUC和95%CI
auc_ci_list_xgb <- lapply(roc_list_xgb, function(roc_obj) {
  ci_obj <- ci.auc(roc_obj)
  c(CI_lower = ci_obj[1], AUC = ci_obj[2], CI_upper = ci_obj[3])
})

auc_ci_df_xgb <- do.call(rbind, auc_ci_list_xgb) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Class")

# Step 5: 构造绘图数据
roc_df_xgb <- lapply(names(roc_list_xgb), function(class) {
  data.frame(
    Spec = 1 - roc_list_xgb[[class]]$specificities,
    Sens = roc_list_xgb[[class]]$sensitivities,
    Class = class
  )
}) %>% bind_rows()

# Step 6: 绘图
ggplot(roc_df_xgb, aes(x = Spec, y = Sens, color = Class)) +
  geom_line(size = 3) +
  scale_color_manual(values = subclass_colors) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80", size = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 48),
    axis.title = element_text(color = "black", size = 48)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves (XGBoost, Validation Set)")+
  # 添加 AUC 注释
  annotate("text", x = 0.27, y = 0.18,
           label = paste0("ICM: ", sprintf("%.3f", auc_ci_df_xgb$AUC[auc_ci_df_xgb$Class == "ICM"]),
                          " (", sprintf("%.3f", auc_ci_df_xgb$CI_lower[auc_ci_df_xgb$Class == "ICM"]),
                          "–", sprintf("%.3f", auc_ci_df_xgb$CI_upper[auc_ci_df_xgb$Class == "ICM"]), ")"),
           fontface = "italic", color = subclass_colors["ICM"], size = 16, hjust = 0.01) +
  annotate("text", x = 0.18, y = 0.1,
           label = paste0("DCM: ", sprintf("%.3f", auc_ci_df_xgb$AUC[auc_ci_df_xgb$Class == "DCM"]),
                          " (", sprintf("%.3f", auc_ci_df_xgb$CI_lower[auc_ci_df_xgb$Class == "DCM"]),
                          "–", sprintf("%.3f", auc_ci_df_xgb$CI_upper[auc_ci_df_xgb$Class == "DCM"]), ")"),
           fontface = "italic", color = subclass_colors["DCM"], size = 16, hjust = 0.01) +
  annotate("text", x = 0.23, y = 0.04,
           label = paste0("CTL: ", sprintf("%.3f", auc_ci_df_xgb$AUC[auc_ci_df_xgb$Class == "CTL"]),
                          " (", sprintf("%.3f", auc_ci_df_xgb$CI_lower[auc_ci_df_xgb$Class == "CTL"]),
                          "–", sprintf("%.3f", auc_ci_df_xgb$CI_upper[auc_ci_df_xgb$Class == "CTL"]), ")"),
           fontface = "italic", color = subclass_colors["CTL"], size = 16, hjust = 0.01)
  
  
  ## SVM algorithm
  library(caret)
  library(dplyr)
  library(pROC)
  library(tidyr)
  
  # Step 0: 读取 DEG 频次表
  Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")
  top_ICM <- names(Filter(function(x) x >= 100, Feature_list$ICM))
  top_DCM <- names(Filter(function(x) x >= 100, Feature_list$DCM))
  top_CTL <- names(Filter(function(x) x >= 100, Feature_list$CTL))
  
  top_features <- unique(c(top_ICM, top_DCM, top_CTL))
  top_features <- intersect(top_features, rownames(expr_val_hg38))
  
  # Step 1: 构建输入表达矩阵
  val_expr_input <- t(expr_val_hg38[top_features, ]) %>% as.data.frame()
  val_expr_input$Class <- classes_val$Classes[match(rownames(val_expr_input), classes_val$ID)]
  val_expr_input$Class <- as.factor(val_expr_input$Class)
  
  
  # 参数设置
  n_repeat <- 100
  svm_model_list <- list()
  all_preds_svm <- data.frame()
  
  # 控制参数（和RF/XGB保持一致）
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 1,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = "final"
  )
  
  # 主循环
  set.seed(123)
  for (i in 1:n_repeat) {
    cat("Training SVM Repeat", i, "\n")
    
    set.seed(100 + i)
    
    svm_model <- caret::train(
      Class ~ ., data = val_expr_input,
      method = "svmRadial",  # 可替换为 svmLinear 试试线性核
      trControl = ctrl,
      metric = "Accuracy"
    )
    
    svm_model_list[[i]] <- svm_model
    all_preds_svm <- bind_rows(all_preds_svm, svm_model$pred)
  }
  
  # 混淆矩阵 + 性能
  conf_mat_svm <- confusionMatrix(
    factor(all_preds_svm$pred, levels = c("CTL", "DCM", "ICM")),
    factor(all_preds_svm$obs, levels = c("CTL", "DCM", "ICM"))
  )
  
  cat("===== SVM 总混淆矩阵 =====\n")
  print(conf_mat_svm$table)
  cat("\n===== SVM 总体性能指标 =====\n")
  print(conf_mat_svm$overall)
  cat("\n===== SVM 每类指标 =====\n")
  print(conf_mat_svm$byClass)

saveRDS(all_preds_svm, file = "GSE55296_all_preds_svm.rds")
  
  
  ## ROC
  library(ggplot2)
  library(pROC)
  library(dplyr)
  library(tidyr)
  subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")

  pred_df_svm <- all_preds_svm[, c("obs", "CTL", "DCM", "ICM")]
  
  long_pred_svm <- pred_df_svm %>%
    pivot_longer(cols = -obs, names_to = "Class", values_to = "Prob") %>%
    mutate(Label = ifelse(obs == Class, 1, 0))
  
  roc_list_svm <- long_pred_svm %>%
    group_by(Class) %>%
    group_split() %>%
    lapply(function(df) {
      roc_obj <- roc(response = df$Label, predictor = df$Prob, quiet = TRUE)
      return(roc_obj)
    })
  
  names(roc_list_svm) <- c("CTL", "DCM", "ICM")
  
  auc_ci_list_svm <- lapply(roc_list_svm, function(roc_obj) {
    ci_obj <- ci.auc(roc_obj)
    c(CI_lower = ci_obj[1], AUC = ci_obj[2], CI_upper = ci_obj[3])
  })
  
  auc_ci_df_svm <- do.call(rbind, auc_ci_list_svm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Class")

roc_df_svm <- lapply(names(roc_list_svm), function(class) {
    data.frame(
      Spec = 1 - roc_list_svm[[class]]$specificities,
      Sens = roc_list_svm[[class]]$sensitivities,
      Class = class
    )
  }) %>% bind_rows()

ggplot(roc_df_svm, aes(x = Spec, y = Sens, color = Class)) +
  geom_line(size = 3) +
  scale_color_manual(values = subclass_colors) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80", size = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 48),
    axis.title = element_text(color = "black", size = 48)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves (SVM, Validation Set)") 
  
  # AUC 注释（位置可调整）
  annotate("text", x = 0.27, y = 0.18,
           label = paste0("ICM: ", sprintf("%.3f", auc_ci_df_svm$AUC[auc_ci_df_svm$Class == "ICM"]),
                          " (", sprintf("%.3f", auc_ci_df_svm$CI_lower[auc_ci_df_svm$Class == "ICM"]),
                          "–", sprintf("%.3f", auc_ci_df_svm$CI_upper[auc_ci_df_svm$Class == "ICM"]), ")"),
           fontface = "italic", color = subclass_colors["ICM"], size = 16, hjust = 0.01) +
  
  annotate("text", x = 0.18, y = 0.1,
           label = paste0("DCM: ", sprintf("%.3f", auc_ci_df_svm$AUC[auc_ci_df_svm$Class == "DCM"]),
                          " (", sprintf("%.3f", auc_ci_df_svm$CI_lower[auc_ci_df_svm$Class == "DCM"]),
                          "–", sprintf("%.3f", auc_ci_df_svm$CI_upper[auc_ci_df_svm$Class == "DCM"]), ")"),
           fontface = "italic", color = subclass_colors["DCM"], size = 16, hjust = 0.01) +
  
  annotate("text", x = 0.23, y = 0.04,
           label = paste0("CTL: ", sprintf("%.3f", auc_ci_df_svm$AUC[auc_ci_df_svm$Class == "CTL"]),
                          " (", sprintf("%.3f", auc_ci_df_svm$CI_lower[auc_ci_df_svm$Class == "CTL"]),
                          "–", sprintf("%.3f", auc_ci_df_svm$CI_upper[auc_ci_df_svm$Class == "CTL"]), ")"),
           fontface = "italic", color = subclass_colors["CTL"], size = 16, hjust = 0.01)
  
  
# 雷达图比较
  # 安装包（如尚未安装）
  # install.packages("fmsb")
  
  library(fmsb)
  
  # 创建数据框（包含最大值、最小值以及三种模型的宏平均值）
  # 设置最大值和最小值行为 1 和 0.6
  radar_data <- data.frame(
    AUROC = c(0.9, 0.5, 0.8663, 0.8477, 0.8353),
    Accuracy = c(0.9, 0.5, 0.7890, 0.7740, 0.6840),
    Kappa = c(0.9, 0.5, 0.6840, 0.6600, 0.5240),
    Specificity = c(0.9, 0.5, 0.8930, 0.8843, 0.8407),
    Sensitivity = c(0.9, 0.5, 0.8050, 0.7850, 0.6913),
    Precision = c(0.9, 0.5, 0.7917, 0.7813, 0.6833),
    F1 = c(0.9, 0.5, 0.7960, 0.7827, 0.6863),
    Balanced_Accuracy = c(0.9, 0.6, 0.8490, 0.8267, 0.7660)
  )
  # 设置行名
  rownames(radar_data) <- c("Max", "Min", "RF", "XGBoost", "SVM")
  
  # 设置颜色
  colors_border <- c("#005493", "#DE4967","#7C8587")
  colors_in <- adjustcolor(colors_border, alpha.f = 0.3)
  colnames(radar_data) <- c("AUROC", 
                            "Accuracy", 
                            "Kappa", 
                            "Specificity", 
                            "Sensitivity", 
                            "Precision", 
                            "F1\nScore", 
                            "Accuracy\n(Balanced)")
  # 绘图
  radarchart(radar_data,
             axistype = 0,
             pcol = colors_border,
             pfcol = colors_in,
             plwd = 5,
             plty = 1,
             cglcol = "black",
             cglty = 5,
             axislabcol = "black",
             caxislabels = seq(0.5, 1, 0.1),
             cglwd = 2,
             vlcex = 2)
  

  
  legend("topright",
         legend = c("RF", "XGBoost", "SVM"),
         bty = "n",
         pch = 20,
         col = colors_border,
         text.col = "black",
         cex = 1.2,
         pt.cex = 3)
  

  
  
## hub 基因验证
# === HF_subtyping_v3.R ===
rm(list=ls())
gc()
setwd("~/Science/CVD/CVD_MS_1/data/RNA-Seq/subtype/top560/validation/")

## ================== 基础准备 ==================
library(readxl)
library(biomaRt)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readr)

# --- 读取样本信息 ---
sample_info <- read.table("GSE55296_sample.txt", 
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)

classes_val <- data.frame(
  ID = sample_info$Sample_ID,
  Classes = sample_info$Group,
  batch = "validation",
  stringsAsFactors = FALSE
)

# ================== 1. 表达矩阵 ==================
count_data <- read_excel("GSE55296_count_data.xlsx") %>% as.data.frame()
count_data <- count_data[!is.na(count_data$GeneID), ] %>% distinct(GeneID, .keep_all = TRUE)
rownames(count_data) <- count_data$GeneID
expr_val <- count_data[, sample_info$Sample_ID]

# ================== 2. Gene ID 映射 (hg19 → hg38, ENSEMBL → SYMBOL) ==================
mart_hg19 <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
gene_map <- getBM(
  attributes = c("ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = rownames(expr_val),
  mart = mart_hg19
)

valid_genes <- unique(gene_map$ensembl_gene_id)
expr_val <- expr_val[valid_genes, ]

# logCPM 转换
log_expr <- log2(expr_val + 1)

# --- ENSEMBL → SYMBOL 映射（按表达量取最显著的 SYMBOL） ---
map_val <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(log_expr),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL","SYMBOL")
) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(ENSEMBL) %>%
  mutate(mean_expr = rowMeans(log_expr[ENSEMBL, , drop=FALSE])) %>%
  arrange(desc(mean_expr)) %>%
  slice(1) %>%
  ungroup() %>%
  # 再保证每个 SYMBOL 唯一
  group_by(SYMBOL) %>%
  arrange(desc(mean_expr)) %>%
  slice(1) %>%
  ungroup()

log_expr_symbol <- log_expr[map_val$ENSEMBL, ]
rownames(log_expr_symbol) <- map_val$SYMBOL

# ================== 3. 定义 logFC 计算函数 ==================
calc_logfc_val <- function(subclass, expr, classes){
  case_ids <- classes %>% filter(Classes == subclass) %>% pull(ID)
  ctrl_ids <- classes %>% filter(Classes == "CTL") %>% pull(ID)
  
  if(length(case_ids) > 1 && length(ctrl_ids) > 1){
    case_mean <- rowMeans(expr[, case_ids, drop = FALSE])
    ctrl_mean <- rowMeans(expr[, ctrl_ids, drop = FALSE])
    logFC <- case_mean - ctrl_mean
    return(data.frame(Gene = rownames(expr), logFC = logFC))
  } else {
    return(NULL)
  }
}

# ================== 4. 计算 ICM vs CTL, DCM vs CTL ==================
lfc_val_icm <- calc_logfc_val("ICM", log_expr_symbol, classes_val)
lfc_val_dcm <- calc_logfc_val("DCM", log_expr_symbol, classes_val)

colnames(lfc_val_icm)[2] <- "Val_ICM_logFC"
colnames(lfc_val_dcm)[2] <- "Val_DCM_logFC"

# ================== 5. 读入 hub gene ==================
icm_tbl <- read_csv("ICM_hub_top60_edges_nodes.csv", show_col_types = FALSE)
dcm_tbl <- read_csv("DCM_hub_top60_edges_nodes.csv", show_col_types = FALSE)

library(dplyr)

icm_hubs <- icm_tbl %>%
  arrange(desc(Degree)) %>%
  slice_head(n = 15) %>%
  dplyr::select(Gene, Cluster, Degree) %>%   # ✅ 显式指定 dplyr
  mutate(Disease = "ICM")

dcm_hubs <- dcm_tbl %>%
  arrange(desc(Degree)) %>%
  slice_head(n = 15) %>%
  dplyr::select(Gene, Cluster, Degree) %>%
  mutate(Disease = "DCM")

hub_nodes <- bind_rows(icm_hubs, dcm_hubs)

# ================== 6. 提取验证队列 hub gene 的 logFC ==================
lfc_val_all <- hub_nodes %>%
  left_join(lfc_val_icm, by = "Gene") %>%
  left_join(lfc_val_dcm, by = "Gene")

# ================== 7. 输出 ==================
write.csv(lfc_val_all, "ValidationCohort_HubGenes_logFC.csv", row.names = FALSE)
print(lfc_val_all)
