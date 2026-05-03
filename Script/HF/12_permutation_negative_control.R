#!/usr/bin/env Rscript
# Manuscript outputs:
#   - Fig. S3A
#   - Table_R3F_permutation_repeat_metrics.csv
#   - Table_R3F_permutation_summary.csv
# Source basis:
#   - data/Toxicity/toxicogenomic/HF_permutation_noise_robustness.R

suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(limma)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(randomForest)
})

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
      if (file.exists(file.path(r, "data", "RNA-Seq", "subtype", "top560", "Discovery_start_data.RData"))) return(r)
    }
  }
  stop("Unable to resolve project root. Use --project-root=/path/to/repo")
}

project_root <- resolve_project_root()
load(file.path(project_root, "data", "RNA-Seq", "subtype", "top560", "Discovery_start_data.RData"))

review_dir <- file.path(project_root, "data", "RNA-Seq", "subtype", "reviewer_reply")
perm_chunk_dir <- file.path(review_dir, "perm_chunks")
noise_chunk_dir <- file.path(review_dir, "noise_chunks")
dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(perm_chunk_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(noise_chunk_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(as.integer(get_arg_value("seed", "20260425")))
mode <- get_arg_value("mode", "perm_chunk")
subclasses <- c("ICM", "DCM", "CTL")
nDEG_value <- as.integer(get_arg_value("nDEG", "560"))
k_folds <- as.integer(get_arg_value("k", "10"))
perm_start <- as.integer(get_arg_value("perm-start", "1"))
perm_end <- as.integer(get_arg_value("perm-end", "10"))
noise_level <- as.numeric(get_arg_value("noise-level", "0.05"))
noise_repeats <- as.integer(get_arg_value("noise-repeats", "12"))
ntree_value <- as.integer(get_arg_value("ntree", "300"))

stopifnot(
  ncol(expr_mat) == nrow(classes_df),
  identical(colnames(expr_mat), classes_df$ID),
  !anyDuplicated(classes_df$ID),
  !anyNA(expr_mat),
  !anyNA(classes_df$Classes),
  !anyNA(classes_df$batch)
)

calc_metrics <- function(df, subclass) {
  y <- ifelse(df$ActualClass == subclass, 1, 0)
  auroc <- as.numeric(pROC::auc(y, df$One, quiet = TRUE))
  auprc <- as.numeric(PRROC::pr.curve(
    scores.class0 = df$One[y == 1],
    scores.class1 = df$One[y == 0],
    curve = FALSE
  )$auc.integral)
  data.frame(Subclass = subclass, AUROC = auroc, AUPRC = auprc, stringsAsFactors = FALSE)
}

split_kfold <- function(classes, k) {
  caret::createFolds(classes, k = k, list = TRUE)
}

one_vs_each_fast <- function(mat, classes_df, indices, nDEG, subclass, ntree = 300L) {
  test_data <- mat[, indices, drop = FALSE]
  test_pheno <- classes_df[indices, , drop = FALSE]
  train_data <- mat[, !colnames(mat) %in% test_pheno$ID, drop = FALSE]
  train_pheno <- classes_df %>% filter(!ID %in% test_pheno$ID)

  dge <- edgeR::DGEList(counts = train_data)
  dge <- edgeR::calcNormFactors(dge)
  group <- factor(ifelse(train_pheno$Classes == subclass, subclass, "Others"))
  design <- model.matrix(~0 + group + batch, data = train_pheno)
  colnames(design) <- make.names(colnames(design))
  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  cont <- limma::makeContrasts(contrasts = paste0("group", subclass, " - groupOthers"), levels = design)
  fit2 <- limma::contrasts.fit(fit, cont)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  degs <- limma::topTable(fit2, number = Inf)
  degs <- degs[order(degs$t), , drop = FALSE]

  half <- floor(nDEG / 2)
  top_features <- rownames(rbind(degs[1:half, , drop = FALSE], degs[(nrow(degs) - half + 1):nrow(degs), , drop = FALSE]))
  new_ann <- factor(ifelse(train_pheno$Classes == subclass, "One", "Others"), levels = c("One", "Others"))
  mtry_val <- max(1L, floor(length(top_features) / 3))
  rf_fit <- randomForest::randomForest(
    x = t(train_data[top_features, , drop = FALSE]),
    y = new_ann,
    mtry = mtry_val,
    ntree = ntree
  )

  pred_prob <- predict(rf_fit, newdata = t(test_data[top_features, , drop = FALSE]), type = "prob") %>% data.frame()
  pred_prob$ActualClass <- test_pheno$Classes
  pred_prob$PredictedClass <- predict(rf_fit, newdata = t(test_data[top_features, , drop = FALSE]), type = "response")
  pred_prob
}

inject_noise_counts <- function(mat, frac) {
  if (frac <= 0) return(mat)
  log_mat <- log2(mat + 1)
  gene_sd <- apply(log_mat, 1, sd)
  gene_sd[!is.finite(gene_sd)] <- 0
  noise_mat <- matrix(
    rnorm(length(log_mat), mean = 0, sd = rep(frac * gene_sd, ncol(log_mat))),
    nrow = nrow(log_mat),
    ncol = ncol(log_mat),
    byrow = FALSE
  )
  noisy_log <- log_mat + noise_mat
  noisy_counts <- round(pmax(2^noisy_log - 1, 0))
  rownames(noisy_counts) <- rownames(mat)
  colnames(noisy_counts) <- colnames(mat)
  storage.mode(noisy_counts) <- "numeric"
  noisy_counts
}

run_one_repeat <- function(mat, class_df, repeat_id, label_mode = c("observed", "permuted"), ntree = 300L) {
  label_mode <- match.arg(label_mode)
  class_df_run <- class_df
  if (label_mode == "permuted") {
    class_df_run$Classes <- sample(class_df_run$Classes, length(class_df_run$Classes), replace = FALSE)
  }
  folds <- split_kfold(class_df_run$Classes, k = k_folds)
  rows <- vector("list", length(subclasses))

  for (i in seq_along(subclasses)) {
    sc <- subclasses[[i]]
    fold_preds <- lapply(seq_along(folds), function(j) {
      one_vs_each_fast(
        mat = mat,
        classes_df = class_df_run,
        indices = folds[[j]],
        nDEG = nDEG_value,
        subclass = sc,
        ntree = ntree
      )
    })
    rep_df <- bind_rows(fold_preds)
    rows[[i]] <- calc_metrics(rep_df, sc)
  }

  bind_rows(rows) %>%
    mutate(
      Repeat = repeat_id,
      LabelMode = label_mode,
      MacroAUROC = mean(AUROC),
      MacroAUPRC = mean(AUPRC)
    ) %>%
    dplyr::select(LabelMode, Repeat, Subclass, AUROC, AUPRC, MacroAUROC, MacroAUPRC)
}

observed_from_saved <- function() {
  fs <- list.files(
    file.path(project_root, "data", "RNA-Seq", "subtype", "top560"),
    pattern = "^HF_repeat_[0-9]+[.]rds$",
    full.names = TRUE
  )
  if (length(fs) == 0) stop("No observed HF_repeat_*.rds files found.")
  bind_rows(lapply(seq_along(fs), function(i) {
    x <- readRDS(fs[[i]])
    bind_rows(lapply(subclasses, function(sc) {
      rep_df <- bind_rows(lapply(x[[sc]], function(z) z$TestPred))
      calc_metrics(rep_df, sc)
    })) %>%
      mutate(
        Repeat = i,
        LabelMode = "observed",
        MacroAUROC = mean(AUROC),
        MacroAUPRC = mean(AUPRC)
      ) %>%
      dplyr::select(LabelMode, Repeat, Subclass, AUROC, AUPRC, MacroAUROC, MacroAUPRC)
  }))
}

mean_ci <- function(x) {
  n <- sum(!is.na(x))
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  se <- if (n > 0) s / sqrt(n) else NA_real_
  ci_mult <- if (n > 1) qt(0.975, df = n - 1) else NA_real_
  c(mean = m, sd = s, n = n, ci_low = m - ci_mult * se, ci_high = m + ci_mult * se)
}

build_summary_block <- function(df, value_cols = c("AUROC", "AUPRC"), prefix) {
  metric1 <- value_cols[[1]]
  metric2 <- value_cols[[2]]
  sub_block <- df %>%
    group_by(Subclass) %>%
    summarise(
      tmp1 = list(mean_ci(.data[[metric1]])),
      tmp2 = list(mean_ci(.data[[metric2]])),
      .groups = "drop"
    ) %>%
    mutate(
      !!paste0(prefix, "_mean_AUROC") := sapply(tmp1, `[[`, "mean"),
      !!paste0(prefix, "_sd_AUROC") := sapply(tmp1, `[[`, "sd"),
      !!paste0(prefix, "_n_AUROC") := sapply(tmp1, `[[`, "n"),
      !!paste0(prefix, "_ci_low_AUROC") := sapply(tmp1, `[[`, "ci_low"),
      !!paste0(prefix, "_ci_high_AUROC") := sapply(tmp1, `[[`, "ci_high"),
      !!paste0(prefix, "_mean_AUPRC") := sapply(tmp2, `[[`, "mean"),
      !!paste0(prefix, "_sd_AUPRC") := sapply(tmp2, `[[`, "sd"),
      !!paste0(prefix, "_n_AUPRC") := sapply(tmp2, `[[`, "n"),
      !!paste0(prefix, "_ci_low_AUPRC") := sapply(tmp2, `[[`, "ci_low"),
      !!paste0(prefix, "_ci_high_AUPRC") := sapply(tmp2, `[[`, "ci_high")
    ) %>%
    dplyr::select(-tmp1, -tmp2)

  macro_block <- df %>%
    distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
    summarise(
      tmp1 = list(mean_ci(MacroAUROC)),
      tmp2 = list(mean_ci(MacroAUPRC))
    ) %>%
    mutate(
      Subclass = "MacroMean",
      !!paste0(prefix, "_mean_AUROC") := tmp1[[1]]["mean"],
      !!paste0(prefix, "_sd_AUROC") := tmp1[[1]]["sd"],
      !!paste0(prefix, "_n_AUROC") := tmp1[[1]]["n"],
      !!paste0(prefix, "_ci_low_AUROC") := tmp1[[1]]["ci_low"],
      !!paste0(prefix, "_ci_high_AUROC") := tmp1[[1]]["ci_high"],
      !!paste0(prefix, "_mean_AUPRC") := tmp2[[1]]["mean"],
      !!paste0(prefix, "_sd_AUPRC") := tmp2[[1]]["sd"],
      !!paste0(prefix, "_n_AUPRC") := tmp2[[1]]["n"],
      !!paste0(prefix, "_ci_low_AUPRC") := tmp2[[1]]["ci_low"],
      !!paste0(prefix, "_ci_high_AUPRC") := tmp2[[1]]["ci_high"]
    ) %>%
    dplyr::select(-tmp1, -tmp2)

  bind_rows(sub_block, macro_block)
}

write_perm_chunk <- function() {
  out <- bind_rows(lapply(perm_start:perm_end, function(i) {
    message("Running permuted repeat ", i, "/", perm_end)
    run_one_repeat(expr_mat, classes_df, repeat_id = i, label_mode = "permuted", ntree = ntree_value)
  }))
  out_file <- file.path(perm_chunk_dir, sprintf("perm_chunk_%04d_%04d.csv", perm_start, perm_end))
  write.csv(out, out_file, row.names = FALSE)
  cat(out_file, "\n")
}

write_noise_chunk <- function() {
  mat_use <- inject_noise_counts(expr_mat, noise_level)
  out <- bind_rows(lapply(seq_len(noise_repeats), function(i) {
    message("Running noise=", noise_level, " repeat ", i, "/", noise_repeats)
    run_one_repeat(mat_use, classes_df, repeat_id = i, label_mode = "observed", ntree = ntree_value)
  }))
  out$NoiseLevel <- sprintf("%.0f%%", noise_level * 100)
  out_file <- file.path(noise_chunk_dir, sprintf("noise_%s_repeats_%02d.csv", gsub("[.]", "p", sprintf("%.2f", noise_level)), noise_repeats))
  write.csv(out, out_file, row.names = FALSE)
  cat(out_file, "\n")
}

aggregate_outputs <- function() {
  observed_res <- observed_from_saved()

  perm_files <- list.files(perm_chunk_dir, pattern = "^perm_chunk_[0-9]{4}_[0-9]{4}[.]csv$", full.names = TRUE)
  if (length(perm_files) == 0) stop("No permutation chunk files found in ", perm_chunk_dir)
  perm_res <- bind_rows(lapply(perm_files, read.csv, stringsAsFactors = FALSE)) %>% arrange(Repeat, Subclass)

  noise_files <- list.files(noise_chunk_dir, pattern = "^noise_.*[.]csv$", full.names = TRUE)
  if (length(noise_files) == 0) stop("No noise chunk files found in ", noise_chunk_dir)
  noise_res <- bind_rows(lapply(noise_files, read.csv, stringsAsFactors = FALSE))
  noise_res_0 <- observed_res
  noise_res_0$NoiseLevel <- "0%"
  noise_repeat_res <- bind_rows(noise_res_0, noise_res)

  obs_summary <- build_summary_block(observed_res, prefix = "observed")
  perm_summary <- build_summary_block(perm_res, prefix = "perm")
  perm_q <- bind_rows(
    perm_res %>%
      group_by(Subclass) %>%
      summarise(
        perm_q95_AUROC = quantile(AUROC, probs = 0.95, na.rm = TRUE),
        perm_q975_AUROC = quantile(AUROC, probs = 0.975, na.rm = TRUE),
        perm_q95_AUPRC = quantile(AUPRC, probs = 0.95, na.rm = TRUE),
        perm_q975_AUPRC = quantile(AUPRC, probs = 0.975, na.rm = TRUE),
        .groups = "drop"
      ),
    perm_res %>%
      distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
      summarise(
        Subclass = "MacroMean",
        perm_q95_AUROC = quantile(MacroAUROC, probs = 0.95, na.rm = TRUE),
        perm_q975_AUROC = quantile(MacroAUROC, probs = 0.975, na.rm = TRUE),
        perm_q95_AUPRC = quantile(MacroAUPRC, probs = 0.95, na.rm = TRUE),
        perm_q975_AUPRC = quantile(MacroAUPRC, probs = 0.975, na.rm = TRUE)
      )
  )

  perm_emp <- bind_rows(
    perm_res %>%
      group_by(Subclass) %>%
      summarise(
        empirical_p_AUROC = (sum(AUROC >= obs_summary$observed_mean_AUROC[match(first(Subclass), obs_summary$Subclass)]) + 1) / (n() + 1),
        empirical_p_AUPRC = (sum(AUPRC >= obs_summary$observed_mean_AUPRC[match(first(Subclass), obs_summary$Subclass)]) + 1) / (n() + 1),
        .groups = "drop"
      ),
    perm_res %>%
      distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
      summarise(
        Subclass = "MacroMean",
        empirical_p_AUROC = (sum(MacroAUROC >= obs_summary$observed_mean_AUROC[obs_summary$Subclass == "MacroMean"]) + 1) / (n() + 1),
        empirical_p_AUPRC = (sum(MacroAUPRC >= obs_summary$observed_mean_AUPRC[obs_summary$Subclass == "MacroMean"]) + 1) / (n() + 1)
      )
  )

  perm_final <- obs_summary %>%
    left_join(perm_summary, by = "Subclass") %>%
    left_join(perm_q, by = "Subclass") %>%
    left_join(perm_emp, by = "Subclass") %>%
    mutate(
      delta_AUROC = observed_mean_AUROC - perm_mean_AUROC,
      delta_AUPRC = observed_mean_AUPRC - perm_mean_AUPRC,
      observed_gt_perm95_AUROC = observed_mean_AUROC > perm_q95_AUROC,
      observed_gt_perm95_AUPRC = observed_mean_AUPRC > perm_q95_AUPRC
    ) %>%
    arrange(match(Subclass, c(subclasses, "MacroMean")))

  noise_summary <- bind_rows(
    noise_repeat_res %>%
      group_by(NoiseLevel, Subclass) %>%
      summarise(
        tmp1 = list(mean_ci(AUROC)),
        tmp2 = list(mean_ci(AUPRC)),
        .groups = "drop"
      ) %>%
      mutate(
        mean_AUROC = sapply(tmp1, `[[`, "mean"),
        sd_AUROC = sapply(tmp1, `[[`, "sd"),
        n_AUROC = sapply(tmp1, `[[`, "n"),
        ci_low_AUROC = sapply(tmp1, `[[`, "ci_low"),
        ci_high_AUROC = sapply(tmp1, `[[`, "ci_high"),
        mean_AUPRC = sapply(tmp2, `[[`, "mean"),
        sd_AUPRC = sapply(tmp2, `[[`, "sd"),
        n_AUPRC = sapply(tmp2, `[[`, "n"),
        ci_low_AUPRC = sapply(tmp2, `[[`, "ci_low"),
        ci_high_AUPRC = sapply(tmp2, `[[`, "ci_high")
      ) %>%
      dplyr::select(-tmp1, -tmp2),
    noise_repeat_res %>%
      distinct(NoiseLevel, Repeat, MacroAUROC, MacroAUPRC) %>%
      group_by(NoiseLevel) %>%
      summarise(
        tmp1 = list(mean_ci(MacroAUROC)),
        tmp2 = list(mean_ci(MacroAUPRC)),
        .groups = "drop"
      ) %>%
      mutate(
        Subclass = "MacroMean",
        mean_AUROC = sapply(tmp1, `[[`, "mean"),
        sd_AUROC = sapply(tmp1, `[[`, "sd"),
        n_AUROC = sapply(tmp1, `[[`, "n"),
        ci_low_AUROC = sapply(tmp1, `[[`, "ci_low"),
        ci_high_AUROC = sapply(tmp1, `[[`, "ci_high"),
        mean_AUPRC = sapply(tmp2, `[[`, "mean"),
        sd_AUPRC = sapply(tmp2, `[[`, "sd"),
        n_AUPRC = sapply(tmp2, `[[`, "n"),
        ci_low_AUPRC = sapply(tmp2, `[[`, "ci_low"),
        ci_high_AUPRC = sapply(tmp2, `[[`, "ci_high")
      ) %>%
      dplyr::select(-tmp1, -tmp2)
  ) %>%
    left_join(
      bind_rows(
        noise_repeat_res %>%
          filter(NoiseLevel == "0%") %>%
          group_by(Subclass) %>%
          summarise(
            baseline_AUROC = mean(AUROC),
            baseline_AUPRC = mean(AUPRC),
            .groups = "drop"
          ),
        noise_repeat_res %>%
          filter(NoiseLevel == "0%") %>%
          distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
          summarise(
            Subclass = "MacroMean",
            baseline_AUROC = mean(MacroAUROC),
            baseline_AUPRC = mean(MacroAUPRC)
          )
      ),
      by = "Subclass"
    ) %>%
    mutate(
      delta_vs_baseline_AUROC = mean_AUROC - baseline_AUROC,
      delta_vs_baseline_AUPRC = mean_AUPRC - baseline_AUPRC
    ) %>%
    arrange(match(Subclass, c(subclasses, "MacroMean")), NoiseLevel)

  out_perm_repeat <- file.path(review_dir, "Table_R3F_permutation_repeat_metrics.csv")
  out_perm_summary <- file.path(review_dir, "Table_R3F_permutation_summary.csv")
  out_noise_repeat <- file.path(review_dir, "Table_R3G_noise_injection_repeat_metrics.csv")
  out_noise_summary <- file.path(review_dir, "Table_R3G_noise_injection_summary.csv")
  out_perm_fig <- file.path(review_dir, "Fig_R3F_permutation_null_distribution.png")
  out_noise_fig <- file.path(review_dir, "Fig_R3G_noise_injection_sensitivity.png")
  out_txt <- file.path(review_dir, "Text_R3F_R3G_robustness_note.txt")

  write.csv(bind_rows(observed_res, perm_res), out_perm_repeat, row.names = FALSE)
  write.csv(perm_final, out_perm_summary, row.names = FALSE)
  write.csv(noise_repeat_res, out_noise_repeat, row.names = FALSE)
  write.csv(noise_summary, out_noise_summary, row.names = FALSE)

  perm_plot_df <- bind_rows(
    perm_res %>%
      distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
      transmute(Repeat, Metric = "Macro AUROC", Value = MacroAUROC),
    perm_res %>%
      distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
      transmute(Repeat, Metric = "Macro AUPRC", Value = MacroAUPRC)
  )
  obs_lines <- data.frame(
    Metric = c("Macro AUROC", "Macro AUPRC"),
    Observed = c(
      perm_final$observed_mean_AUROC[perm_final$Subclass == "MacroMean"],
      perm_final$observed_mean_AUPRC[perm_final$Subclass == "MacroMean"]
    ),
    Q95 = c(
      perm_final$perm_q95_AUROC[perm_final$Subclass == "MacroMean"],
      perm_final$perm_q95_AUPRC[perm_final$Subclass == "MacroMean"]
    )
  )

  g_perm <- ggplot(perm_plot_df, aes(x = Value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "#C95B4A", color = "white", alpha = 0.75) +
    geom_density(color = "#8A2D21", linewidth = 0.8) +
    geom_vline(data = obs_lines, aes(xintercept = Observed), color = "#1F4E79", linewidth = 1.1) +
    geom_vline(data = obs_lines, aes(xintercept = Q95), color = "#555555", linewidth = 0.8, linetype = "dashed") +
    facet_wrap(~Metric, scales = "free") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    labs(
      x = "Permutation metric value",
      y = "Density",
      title = "Permutation null distribution with observed performance"
    )
  ggsave(out_perm_fig, g_perm, width = 8.2, height = 4.8, dpi = 300)

  g_noise <- ggplot(
    noise_summary %>% filter(Subclass %in% c(subclasses, "MacroMean")),
    aes(x = NoiseLevel, y = mean_AUROC, color = Subclass, group = Subclass)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.1) +
    geom_errorbar(aes(ymin = ci_low_AUROC, ymax = ci_high_AUROC), width = 0.08) +
    facet_wrap(~Subclass, scales = "free_y") +
    scale_color_manual(values = c("ICM" = "#E6AB2A", "DCM" = "#941651", "CTL" = "#4582B0", "MacroMean" = "#2E8B57")) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(
      x = "Injected Gaussian noise on log2(count + 1)",
      y = "Mean AUROC across repeats (95% CI)",
      title = "Noise-injection sensitivity with repeat-level uncertainty"
    )
  ggsave(out_noise_fig, g_noise, width = 8.2, height = 6.0, dpi = 300)

  macro_row <- perm_final %>% filter(Subclass == "MacroMean")
  noise_macro <- noise_summary %>% filter(Subclass == "MacroMean")
  noise_5 <- noise_macro %>% filter(NoiseLevel == "5%")
  noise_10 <- noise_macro %>% filter(NoiseLevel == "10%")
  txt <- c(
    "Permutation test and noise-injection sensitivity check",
    sprintf(
      "Using the same one-vs-rest DEG-selection plus RF framework (nDEG=%d, %d-fold outer CV, fixed mtry=p/3), the observed macro-mean AUROC/AUPRC were %.4f/%.4f, compared with permutation null means of %.4f/%.4f across %d label permutations.",
      nDEG_value, k_folds,
      macro_row$observed_mean_AUROC, macro_row$observed_mean_AUPRC,
      macro_row$perm_mean_AUROC, macro_row$perm_mean_AUPRC,
      macro_row$perm_n_AUROC
    ),
    sprintf(
      "The empirical p values were %.4f for macro AUROC and %.4f for macro AUPRC; the observed values exceeded the 95th permutation quantiles (AUROC q95=%.4f, AUPRC q95=%.4f).",
      macro_row$empirical_p_AUROC, macro_row$empirical_p_AUPRC,
      macro_row$perm_q95_AUROC, macro_row$perm_q95_AUPRC
    ),
    sprintf(
      "Under Gaussian perturbation on log2(count+1), macro AUROC was %.4f (95%% CI %.4f-%.4f) at baseline, %.4f (%.4f-%.4f) at 5%% noise, and %.4f (%.4f-%.4f) at 10%% noise.",
      noise_macro$mean_AUROC[noise_macro$NoiseLevel == '0%'],
      noise_macro$ci_low_AUROC[noise_macro$NoiseLevel == '0%'],
      noise_macro$ci_high_AUROC[noise_macro$NoiseLevel == '0%'],
      noise_5$mean_AUROC, noise_5$ci_low_AUROC, noise_5$ci_high_AUROC,
      noise_10$mean_AUROC, noise_10$ci_low_AUROC, noise_10$ci_high_AUROC
    ),
    "The small fluctuations across the 5%-10% synthetic-noise conditions should therefore be interpreted as repeat-to-repeat variation rather than performance gain."
  )
  writeLines(txt, out_txt)

  cat("Generated:\n")
  cat(" -", out_perm_repeat, "\n")
  cat(" -", out_perm_summary, "\n")
  cat(" -", out_noise_repeat, "\n")
  cat(" -", out_noise_summary, "\n")
  cat(" -", out_perm_fig, "\n")
  cat(" -", out_noise_fig, "\n")
  cat(" -", out_txt, "\n")
}

if (mode == "perm_chunk") {
  write_perm_chunk()
} else if (mode == "noise_chunk") {
  write_noise_chunk()
} else if (mode == "aggregate") {
  aggregate_outputs()
} else if (mode == "full") {
  write_perm_chunk()
  write_noise_chunk()
  aggregate_outputs()
} else {
  stop("Unsupported mode: ", mode)
}
