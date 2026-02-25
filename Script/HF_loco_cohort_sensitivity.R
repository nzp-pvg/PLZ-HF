#!/usr/bin/env Rscript
# Author: Zhaoxian Wang


suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(pROC)
})

set.seed(2026)

get_arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(default)
  sub(key, "", hit[[1]], fixed = TRUE)
}

resolve_project_root <- function() {
  from_arg <- get_arg_value("project-root", NA_character_)
  from_env <- Sys.getenv("CVD_MS1_PROJECT_ROOT", unset = "")
  cands <- unique(c(if (!is.na(from_arg)) from_arg else NULL, if (nzchar(from_env)) from_env else NULL, getwd()))
  for (s in cands) {
    p <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in 0:6) {
      r <- normalizePath(file.path(p, paste(rep("..", i), collapse = "/")), winslash = "/", mustWork = FALSE)
      if (file.exists(file.path(r, "data", "RNA-Seq", "subtype", "one_vs_each_HF.R"))) return(r)
    }
  }
  stop("Unable to resolve project root. Use --project-root=/path/to/repo")
}

project_root <- resolve_project_root()
review_dir <- file.path(project_root, "data", "RNA-Seq", "subtype", "reviewer_reply")
dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "data", "RNA-Seq", "subtype", "one_vs_each_HF.R"))
load(file.path(project_root, "data", "RNA-Seq", "subtype", "top560", "Discovery_start_data.RData")) # expr_mat, classes_df

subclasses <- c("ICM", "DCM", "CTL")
cohorts <- sort(unique(classes_df$batch))
nDEG_value <- as.integer(get_arg_value("nDEG", "560"))

safe_bal_acc <- function(actual, pred, labels = subclasses) {
  rec <- sapply(labels, function(lb) {
    idx <- actual == lb
    if (sum(idx) == 0) return(NA_real_)
    mean(pred[idx] == lb)
  })
  mean(rec, na.rm = TRUE)
}

safe_auc <- function(actual, prob, pos_label) {
  y <- ifelse(actual == pos_label, 1, 0)
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(y, prob, quiet = TRUE))
}

cv_ctrl <- trainControl(method = "none", classProbs = TRUE)

loco_summary <- list()
loco_binary <- list()
loco_preds <- list()
idx_s <- 1
idx_b <- 1
idx_p <- 1

for (h in cohorts) {
  cat("Running LOCO holdout cohort:", h, "\n")
  test_ids <- classes_df$ID[classes_df$batch == h]

  sub_pred_list <- list()
  sub_auc <- c()

  for (sc in subclasses) {
    test_idx <- which(classes_df$ID %in% test_ids)
    fit <- OnevsEach.HF(
      Mat = expr_mat,
      classes.df = classes_df,
      Indices = test_idx,
      nDEG = nDEG_value,
      subclass = sc,
      cv_control = cv_ctrl
    )

    pred_df <- fit$TestPred
    pred_df$SampleID <- rownames(pred_df)
    pred_df$Subclass <- sc
    sub_pred_list[[sc]] <- pred_df

    auc_sc <- safe_auc(pred_df$ActualClass, pred_df$One, sc)
    sub_auc[sc] <- auc_sc

    loco_binary[[idx_b]] <- data.frame(
      HoldoutCohort = h,
      Subclass = sc,
      n_test = nrow(pred_df),
      n_positive_in_test = sum(pred_df$ActualClass == sc),
      n_negative_in_test = sum(pred_df$ActualClass != sc),
      AUROC = auc_sc,
      stringsAsFactors = FALSE
    )
    idx_b <- idx_b + 1
  }

  merged <- NULL
  for (sc in subclasses) {
    tmp <- sub_pred_list[[sc]] %>%
      select(SampleID, ActualClass, One) %>%
      rename(!!paste0("Prob_", sc) := One)
    merged <- if (is.null(merged)) tmp else full_join(merged, tmp, by = c("SampleID", "ActualClass"))
  }

  prob_mat <- as.matrix(merged[, paste0("Prob_", subclasses), drop = FALSE])
  merged$PredictedClass <- subclasses[max.col(prob_mat, ties.method = "first")]

  acc <- mean(merged$PredictedClass == merged$ActualClass)
  bal <- safe_bal_acc(merged$ActualClass, merged$PredictedClass)

  loco_summary[[idx_s]] <- data.frame(
    HoldoutCohort = h,
    n_test = nrow(merged),
    n_CTL = sum(merged$ActualClass == "CTL"),
    n_DCM = sum(merged$ActualClass == "DCM"),
    n_ICM = sum(merged$ActualClass == "ICM"),
    Accuracy = acc,
    BalancedAccuracy = bal,
    MeanBinaryAUROC_Available = mean(sub_auc, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  idx_s <- idx_s + 1

  merged$HoldoutCohort <- h
  loco_preds[[idx_p]] <- merged
  idx_p <- idx_p + 1
}

loco_summary_df <- bind_rows(loco_summary)
loco_binary_df <- bind_rows(loco_binary)
loco_preds_df <- bind_rows(loco_preds)

overall <- data.frame(
  HoldoutCohort = "Overall_mean",
  n_test = sum(loco_summary_df$n_test),
  n_CTL = sum(loco_summary_df$n_CTL),
  n_DCM = sum(loco_summary_df$n_DCM),
  n_ICM = sum(loco_summary_df$n_ICM),
  Accuracy = mean(loco_summary_df$Accuracy, na.rm = TRUE),
  BalancedAccuracy = mean(loco_summary_df$BalancedAccuracy, na.rm = TRUE),
  MeanBinaryAUROC_Available = mean(loco_summary_df$MeanBinaryAUROC_Available, na.rm = TRUE),
  stringsAsFactors = FALSE
)

loco_summary_df_out <- bind_rows(loco_summary_df, overall)

out_summary <- file.path(review_dir, "Table_R3A_LOCO_cohort_summary.csv")
out_auc <- file.path(review_dir, "Table_R3A_LOCO_binary_auc_by_subclass.csv")
out_preds <- file.path(review_dir, "Table_R3A_LOCO_predictions.csv")

write.csv(loco_summary_df_out, out_summary, row.names = FALSE)
write.csv(loco_binary_df, out_auc, row.names = FALSE)
write.csv(loco_preds_df, out_preds, row.names = FALSE)

cat("Generated:\n")
cat(" -", out_summary, "\n")
cat(" -", out_auc, "\n")
cat(" -", out_preds, "\n")
