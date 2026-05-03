# Manuscript outputs:
#   - Fig. S1A (TopNumber sensitivity scan)
# Inputs:
#   - Discovery_start_data.RData
#   - 03_one_vs_each_training_function.R
# Outputs:
#   - TopNumber scan result tables/plots in the working HF analysis directory
# Source basis:
#   - PLZ-HF-main/Script/HF_nDEG_topN.R

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

get_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

resolve_project_root <- function() {
  script_dir <- get_script_dir()
  cands <- unique(c(Sys.getenv("CVD_MS1_PROJECT_ROOT", unset = ""), script_dir, getwd()))
  cands <- cands[nzchar(cands)]
  for (seed in cands) {
    cur <- normalizePath(seed, winslash = "/", mustWork = FALSE)
    for (i in 0:8) {
      probe <- normalizePath(file.path(cur, paste(rep("..", i), collapse = "/")), winslash = "/", mustWork = FALSE)
      if (dir.exists(file.path(probe, "data", "RNA-Seq", "subtype", "top560"))) return(probe)
    }
  }
  stop("Unable to resolve project root.")
}

project_root <- resolve_project_root()
script_dir <- get_script_dir()
hf_dir <- file.path(project_root, "data", "RNA-Seq", "subtype", "top560")
setwd(hf_dir)

load("Discovery_start_data.RData")
source(file.path(script_dir, "03_one_vs_each_training_function.R"))

set.seed(2025)
k <- 10
repeat_n <- 50
subclasses <- c("ICM", "DCM", "CTL")

cv_ctrl <- trainControl(
  method = "repeatedcv", number = 10, repeats = 3,
  classProbs = TRUE, verboseIter = TRUE, returnData = FALSE,
  savePredictions = "final"
)

cl <- makeCluster(12)
registerDoParallel(cl)

nDEG_range <- c(100,200,300,400,500,600,700,1000,1200,1500,2000)

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
  
  saveRDS(All.Repeats, file = paste0("HF_CV_results_nDEG_top_", n, ".rds"))
}

stopCluster(cl)

# AUC analysis

nDEG_range <- c(100, 200, 300, 400, 500, 600, 700, 1000, 1200, 1500, 2000)
subclasses <- c("ICM", "DCM", "CTL")
extract_auc_from_test <- function(file_path, subclass) {
  res <- readRDS(file_path)
  aucs <- unlist(
    lapply(res, function(rep) {
      sapply(rep[[subclass]], function(fold) {
        prob <- fold$TestPred
        actual <- ifelse(prob$ActualClass == subclass, 1, 0)
        pred <- prob$One
        pROC::auc(pROC::roc(actual, pred))
      })
    })
  )
  return(aucs)
}

library(pROC)
auc_data <- data.frame()

for (n in nDEG_range) {
  file <- paste0("HF_CV_results_nDEG_top_", n, ".rds")
  for (sub in subclasses) {
    aucs <- extract_auc_from_test(file, sub)
    df <- data.frame(
      TopN = n,
      Subclass = sub,
      AUC = aucs
    )
    auc_data <- rbind(auc_data, df)
  }
}

library(dplyr)

auc_summary <- auc_data %>%
  group_by(TopN, Subclass) %>%
  summarise(
    AUC_mean = round(mean(AUC), 3),
    AUC_sd = round(sd(AUC), 3),
    .groups = "drop"
  )

write.csv(auc_summary, "AUC_summary_by_TopN.csv", row.names = FALSE)
library(openxlsx)

topn_list <- unique(auc_data$TopN)

topn_list <- unique(auc_data$TopN)

wb <- createWorkbook()

for (n in topn_list) {
  df_sub <- subset(auc_data, TopN == n)
  addWorksheet(wb, sheetName = paste0("TopN_", n))
  writeData(wb, sheet = paste0("TopN_", n), df_sub)
}

saveWorkbook(wb, file = "AUC_by_TopN_AllSheets.xlsx", overwrite = TRUE)



library(ggplot2)
ggplot(auc_data, aes(x = TopN, y = AUC, color = Subclass)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 50) +
  geom_point(position = position_jitter(width = 0, height = 0.003), alpha = 0.3) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  geom_vline(xintercept = 700, linetype = "dashed", color = "grey50") +
  annotate("text", x = 720, y = 0.985, label = "Inflection point", hjust = 0) +
  scale_x_continuous(breaks = c(100, 200, 300, 400, 500, 600, 700, 1000, 1200, 1500, 2000)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "AUC vs Top N DEGs (Zoomed View)",
    x = "Top N DEGs",
    y = "AUC (Test Fold)",
    color = "Subclass"
  )



auc_summary <- auc_data %>%
  group_by(Subclass, TopN) %>%
  summarise(
    MeanAUC = mean(AUC),
    SE = sd(AUC) / sqrt(n())
  )

ggplot(auc_summary, aes(x = TopN, y = MeanAUC, color = Subclass)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = MeanAUC - SE, ymax = MeanAUC + SE), width = 50) +
  labs(
    title = "AUC vs Top N DEGs for Each Subclass",
    x = "Top N DEGs",
    y = "Mean AUC (±SE)",
    color = "HF Subclass"
  ) +
  theme_minimal(base_size = 14)



library(dplyr)

auc_summary <- auc_data %>%
  group_by(TopN, Subclass) %>%
  summarise(
    AUC_mean = round(mean(AUC), 3),
    AUC_sd = round(sd(AUC), 3),
    .groups = "drop"
  )

write.csv(auc_summary, "AUC_summary_by_TopN.csv", row.names = FALSE)







