#rm(list=ls())
#gc()
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)
library(randomForest)
library(glmnet)
library(limma)
library(doParallel)

setwd("~/Science/CVD/CVD_MS_1/data/RNA-Seq/subtype/Model_Tuning")

# 加载数据与函数
load("Discovery_start_data.RData")
source("one_vs_each_HF.R")

set.seed(2025)
k <- 10
repeat_n <- 30  # ⚠️ 已改为30次交叉重复
subclasses <- c("ICM", "DCM", "CTL")

cv_ctrl <- trainControl(
  method = "repeatedcv", number = 10, repeats = 3,
  classProbs = TRUE, verboseIter = TRUE, returnData = FALSE,
  savePredictions = "final"
)

# 启动并行环境
cl <- makeCluster(12)
registerDoParallel(cl)

# ======= 定义 topN 范围 =======
nDEG_range <- c(560, 580)
#nDEG_range <- c(600, 700, 1000)
#nDEG_range <- c(1200, 1500, 2000)

# 主循环：每个 topN 执行一次 sweep（不保存每个 repeat）
for (n in nDEG_range) {
  cat("\n===== Running topN =", n, "=====\n")
  All.Repeats <- list()
  
  for (r in 1:repeat_n) {
    cat("==== Repeat", r, "====\n")
    folds <- SplitkFold(expr_mat, classes_df$Classes, K = k)
    result_one_repeat <- list()
    
    for (sub in subclasses) {
      cat("  ---- Subclass:", sub, "----\n")
      sub_results <- list()
      
      for (i in 1:k) {
        cat("    Fold", i, "...\n")
        fold_result <- OnevsEach.HF(
          Mat = expr_mat,
          classes.df = classes_df,
          Indices = folds$samples[[i]],
          nDEG = n,
          subclass = sub,
          cv_control = cv_ctrl
        )
        sub_results[[paste0("Fold", i)]] <- fold_result
      }
      result_one_repeat[[sub]] <- sub_results
    }
    
    All.Repeats[[paste0("Repeat", r)]] <- result_one_repeat
  }
  
  # ✅ 仅保存每个 topN 的完整 sweep 结果
  saveRDS(All.Repeats, file = paste0("HF_CV_results_nDEG_top_", n, ".rds"))
}

# 关闭并行环境
stopCluster(cl)