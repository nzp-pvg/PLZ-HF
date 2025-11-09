rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(readxl)
library(tibble)
library(purrr)
library(pheatmap)
library(readxl)

# 选定最适统一nFreq
setwd("~/Science/CVD/CVD_MS_1/data/RNA-Seq/subtype/nFreq/nFreq_results")

deg_counts <- sapply(1:42, function(i) {
  deg <- readRDS(paste0("nFreq_", i, ".rds"))
  length(deg)
})

plot(1:42, deg_counts, type = "b", pch = 16,
     xlab = "nFreq threshold", ylab = "Number of DEGs",
     main = "Number of DEGs vs. nFreq threshold")


# 例如选定 nFreq = 25
deg_final <- readRDS("nFreq_1.rds")
ICM_genes <- deg_final$ICM$genes
DCM_genes <- deg_final$DCM$genes
CTL_genes <- deg_final$CTL$genes
ICM_genes
DCM_genes
CTL_genes


deg_final$ICM$nGenes









# load HF subtyping classification result
res <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")






















# 标准曲线分析
# === 载入包 ===
library(dplyr)
library(ggplot2)
library(purrr)
library(dplyr)
library(ggplot2)
library(pROC)
library(scales)
library(boot)
library(rms)

# 提取预测概率 + 真值（假设结构为每个子模型结果存储在 res$ICM、res$DCM、res$CTL）
# 每个子模型是 One-vs-Rest，包含预测概率（prob）、实际标签（true）等
# === 提取每个类别的预测概率 ===
extract_probs <- function(res, subtype = c("ICM", "DCM", "CTL")) {
  subtype <- match.arg(subtype)
  prob_list <- lapply(res, function(rep) {
    lapply(rep[[subtype]], function(fold) {
      df <- fold$TestPred
      df$Sample_ID <- rownames(df)
      df[, c("Sample_ID", "One", "ActualClass")]
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  prob_summary <- prob_list %>%
    group_by(Sample_ID) %>%
    summarise(
      Predicted = mean(One),
      TrueLabel = dplyr::first(ActualClass)
    ) %>%
    ungroup()
  
  colnames(prob_summary)[2] <- paste0("prob_", subtype)
  return(prob_summary)
}

# 提取每个模型的预测概率与标签
prob_ICM <- extract_probs(res, "ICM")
prob_DCM <- extract_probs(res, "DCM")
prob_CTL <- extract_probs(res, "CTL")

# 合并成一个 df_pred
df_pred <- prob_ICM %>%
  inner_join(prob_DCM, by = c("Sample_ID", "TrueLabel")) %>%
  inner_join(prob_CTL, by = c("Sample_ID", "TrueLabel"))

# === 计算每一类的校准曲线数据 ===
get_calibration_df <- function(df, prob_col, class_name, n_bins = 10) {
  df_temp <- df %>%
    mutate(
      Observed = ifelse(TrueLabel == class_name, 1, 0),
      Predicted = .data[[prob_col]],
      bin = cut(Predicted, breaks = quantile(Predicted, probs = seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
    ) %>%
    group_by(bin) %>%
    summarise(
      bin_mean = mean(Predicted),
      obs_rate = mean(Observed),
      count = n(),
      .groups = "drop"
    ) %>%
    mutate(Class = class_name)
  
  return(df_temp)
}

# 分别计算每个类的校准数据
cal_ICM <- get_calibration_df(df_pred, "prob_ICM", "ICM")
cal_DCM <- get_calibration_df(df_pred, "prob_DCM", "DCM")
cal_CTL <- get_calibration_df(df_pred, "prob_CTL", "CTL")

# 合并所有
cal_all <- bind_rows(cal_ICM, cal_DCM, cal_CTL)
cal_all
# === 画合并校准曲线图 ===
ggplot(cal_all, aes(x = bin_mean, y = obs_rate, color = Class)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Calibration Curves for ICM, DCM, CTL",
       x = "Predicted Probability", y = "Observed Frequency") +
  theme_minimal() +
  scale_color_manual(values = c("ICM" = "#E6AB2A", "DCM" = "#941651", "CTL" = "#1F78B4"))

# === 计算每个 bin 的置信区间 ===
bootstrap_ci <- function(df, n_bins = 10, n_boot = 200, class_name = "ICM") {
  df$Observed <- ifelse(df$TrueLabel == class_name, 1, 0)
  df$Predicted <- df[[paste0("prob_", class_name)]]
  df$bin <- cut(df$Predicted,
                breaks = quantile(df$Predicted, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
                include.lowest = TRUE)
  
  bins <- levels(df$bin)
  ci_list <- list()
  
  for (b in bins) {
    bin_df <- df[df$bin == b, ]
    if (nrow(bin_df) < 2) next
    
    boot_fun <- function(d, i) {
      mean(d$Observed[i])
    }
    
    # tryCatch 处理 boot.ci 错误
    boot_out <- tryCatch({
      boot(data = bin_df, statistic = boot_fun, R = n_boot)
    }, error = function(e) return(NULL))
    
    if (is.null(boot_out)) next
    
    ci <- tryCatch({
      boot.ci(boot_out, type = "perc")
    }, error = function(e) return(NULL))
    
    if (is.null(ci) || is.null(ci$percent)) {
      lower_ci <- NA
      upper_ci <- NA
    } else {
      lower_ci <- ci$percent[4]
      upper_ci <- ci$percent[5]
    }
    
    ci_list[[b]] <- data.frame(
      bin = b,
      bin_mean = mean(bin_df$Predicted),
      obs_rate = mean(bin_df$Observed),
      count = nrow(bin_df),
      lower = lower_ci,
      upper = upper_ci,
      Class = class_name
    )
  }
  
  do.call(rbind, ci_list)
}

# === 主函数，合并绘图 ===
plot_calibration_with_ci <- function(df, class_names = c("ICM", "DCM", "CTL"), n_bins = 10, n_boot = 200) {
  all_df <- list()
  brier_scores <- c()
  
  for (class in class_names) {
    cal_df <- bootstrap_ci(df, n_bins = n_bins, n_boot = n_boot, class_name = class)
    all_df[[class]] <- cal_df
    
    # Brier Score
    brier <- mean((ifelse(df$TrueLabel == class, 1, 0) - df[[paste0("prob_", class)]])^2)
    brier_scores[class] <- round(brier, 3)
  }
  
  plot_df <- bind_rows(all_df)
  label_map <- paste0(class_names, " (Brier=", brier_scores, ")")
  names(label_map) <- class_names
  
  # 绘图
  ggplot(plot_df, aes(x = bin_mean, y = obs_rate, color = Class, fill = Class)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(size = 1.1) +
    geom_point(size = 2) +
    geom_text(aes(label = count), vjust = -1.3, size = 3, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("ICM" = "#E69F00", "DCM" = "#CC79A7", "CTL" = "#0072B2"),
                       labels = label_map) +
    scale_fill_manual(values = c("ICM" = "#E69F00", "DCM" = "#CC79A7", "CTL" = "#0072B2")) +
    labs(title = "Calibration Curve with CI Bands",
         x = "Predicted Probability",
         y = "Observed Frequency",
         color = "Class", fill = "Class") +
    ylim(0, 1) + xlim(0, 1) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# === 运行 ===
plot_calibration_with_ci(df_pred)

## 单独绘图
library(dplyr)
library(ggplot2)
library(tidyr)
bootstrap_ci <- function(df, n_bins = 10, n_boot = 200, class_name) {
  df_class <- df %>%
    mutate(Observed = ifelse(TrueLabel == class_name, 1, 0),
           Predicted = .[[paste0("prob_", class_name)]]) %>%
    select(Predicted, Observed)
  df_class$bin <- cut(df_class$Predicted,
                      breaks = quantile(df_class$Predicted, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
                      include.lowest = TRUE)
  bin_levels <- levels(df_class$bin)
  out_df <- data.frame(bin = bin_levels, bin_mean = NA, obs_rate = NA, count = NA, lower = NA, upper = NA)
  for (b in bin_levels) {
    bin_data <- df_class %>% filter(bin == b)
    obs_vals <- bin_data$Observed
    pred_vals <- bin_data$Predicted
    n <- length(obs_vals)
    out_df[out_df$bin == b, "bin_mean"] <- mean(pred_vals)
    out_df[out_df$bin == b, "obs_rate"] <- mean(obs_vals)
    out_df[out_df$bin == b, "count"] <- n
    if (n > 0) {
      # Bootstrap CI
      boots <- replicate(n_boot, mean(sample(obs_vals, replace = TRUE)))
      out_df[out_df$bin == b, "lower"] <- quantile(boots, 0.025)
      out_df[out_df$bin == b, "upper"] <- quantile(boots, 0.975)
    }
  }
  out_df$Class <- class_name
  return(out_df)
}

plot_single_calibration_with_ci <- function(df, class_name, bins = 10, boots = 200) {
  cal_df <- bootstrap_ci(df, n_bins = bins, n_boot = boots, class_name = class_name)
  color_map <- c("ICM" = "#E6ab2a", "DCM" = "#941751", "CTL" = "#0072B2")
  curve_color <- color_map[[class_name]]
  ggplot(cal_df, aes(x = bin_mean, y = obs_rate)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = curve_color, alpha = 0.4) +
    geom_line(color = curve_color, size = 4) +
    geom_point(shape = 21, size = 10, fill = "white", color = "black", alpha = 0.7, stroke = 1) +
    geom_text(aes(label = count), vjust = -1.5, size = 0, color = "black", position = position_nudge(y = -0.02)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 2) +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      title = paste("Calibration Curve with CI:", class_name),
      x = "Predicted Probability",
      y = "Observed Frequency"
    ) +
    theme_minimal(base_size = 36) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.text = element_text(size = 40),
      axis.title = element_text(size = 40),
      plot.title = element_text(size = 18, face = "bold")
    )
}
plot_single_calibration_with_ci(df_pred, "ICM")
plot_single_calibration_with_ci(df_pred, "DCM")
plot_single_calibration_with_ci(df_pred, "CTL")


df_long <- df_pred %>%
  select(Sample_ID, TrueLabel, prob_ICM, prob_DCM, prob_CTL) %>%
  pivot_longer(
    cols = starts_with("prob_"),
    names_to = "Class",
    values_to = "Predicted"
  ) %>%
  mutate(
    Class = gsub("prob_", "", Class),
    Observed = ifelse(TrueLabel == Class, 1, 0)
  )

plot_single_calibration <- function(class_name, df_long, bins = 10, return_data = FALSE) {
  df_class <- df_long %>% filter(Class == class_name)
  
  # 分箱
  df_class$bin <- cut(df_class$Predicted,
                      breaks = quantile(df_class$Predicted, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE),
                      include.lowest = TRUE)
  
  cal_df <- df_class %>%
    group_by(bin) %>%
    summarise(
      bin_mean = mean(Predicted),
      obs_rate = mean(Observed),
      count = n(),
      .groups = "drop"
    )
  
  # 如果需要返回数据而不是画图
  if (return_data) {
    return(cal_df)
  }
  
  # 设置颜色
  color_map <- c("ICM" = "#E6ab2a", "DCM" = "#941751", "CTL" = "#0072B2")
  curve_color <- color_map[[class_name]]
  
  ggplot(cal_df, aes(x = bin_mean, y = obs_rate)) +
    geom_line(color = curve_color, size = 3) +
    geom_point(shape = 21, size = 2, fill = curve_color, color = "black", alpha = 0.7, stroke = 1) +
    geom_text(aes(label = count), vjust = -1.5, size = 7, color = "black", position = position_nudge(y = -0.02)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 2) +
    xlim(0, 1) + ylim(0, 1) +
    scale_color_manual(
      name = "Class",
      values = setNames(curve_color, class_name)
    ) +
    labs(
      title = paste("Calibration Curve:", class_name),
      x = "Predicted Probability",
      y = "Observed Frequency"
    ) +
    theme_minimal(base_size = 36) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.text = element_text(size = 40),
      axis.title = element_text(size = 40),
      plot.title = element_text(size = 18, face = "bold")
    )
}
# Brier Score（越低越好）

plot_single_calibration("ICM", df_long)
plot_single_calibration("DCM", df_long)
plot_single_calibration("CTL", df_long)

brier_icm <- mean((df_long %>% filter(Class == "ICM")) %>%
                    with((Predicted - Observed)^2))

brier_dcm <- mean((df_long %>% filter(Class == "DCM")) %>%
                    with((Predicted - Observed)^2))

brier_ctl <- mean((df_long %>% filter(Class == "CTL")) %>%
                    with((Predicted - Observed)^2))

brier_icm
brier_dcm
brier_ctl


## 相线数目
# 查看每个类各自的分箱点数
cal_df_icm <- plot_single_calibration("ICM", df_long, return_data = TRUE)
cal_df_dcm <- plot_single_calibration("DCM", df_long, return_data = TRUE)
cal_df_ctl <- plot_single_calibration("CTL", df_long, return_data = TRUE)

# 每类的相线（分箱点）数目
nrow(cal_df_icm)  # 一般为10，除非某些bin没有数据
nrow(cal_df_dcm)
nrow(cal_df_ctl)

df_class <- df_long %>% filter(Class == "ICM" & Observed == 1)

# ICM
cal_df_icm <- plot_single_calibration("ICM", df_long, return_data = TRUE)
cat("ICM 校准曲线每个点的相线数目（样本数）:\n")
print(cal_df_icm$count)

# DCM
cal_df_dcm <- plot_single_calibration("DCM", df_long, return_data = TRUE)
cat("\nDCM 校准曲线每个点的相线数目（样本数）:\n")
print(cal_df_dcm$count)

# CTL
cal_df_ctl <- plot_single_calibration("CTL", df_long, return_data = TRUE)
cat("\nCTL 校准曲线每个点的相线数目（样本数）:\n")
print(cal_df_ctl$count)


## 计算ECE/MCE
compute_ece_mce <- function(df_long, bins = 10) {
  result_list <- list()
  
  for (cls in unique(df_long$Class)) {
    df_cls <- df_long %>% filter(Class == cls)
    
    # 分箱
    df_cls$bin <- cut(df_cls$Predicted,
                      breaks = quantile(df_cls$Predicted, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE),
                      include.lowest = TRUE)
    
    cal_df <- df_cls %>%
      group_by(bin) %>%
      summarise(
        bin_pred = mean(Predicted),
        bin_obs = mean(Observed),
        count = n(),
        .groups = "drop"
      ) %>%
      mutate(abs_error = abs(bin_pred - bin_obs),
             weight = count / sum(count))
    
    # 计算 ECE 和 MCE
    ece <- sum(cal_df$abs_error * cal_df$weight)
    mce <- max(cal_df$abs_error)
    
    result_list[[cls]] <- data.frame(
      Class = cls,
      ECE = round(ece, 4),
      MCE = round(mce, 4)
    )
  }
  
  result_df <- do.call(rbind, result_list)
  return(result_df)
}

df_long <- df_pred %>%
  select(Sample_ID, TrueLabel, prob_ICM, prob_DCM, prob_CTL) %>%
  pivot_longer(
    cols = starts_with("prob_"),
    names_to = "Class",
    values_to = "Predicted"
  ) %>%
  mutate(
    Class = gsub("prob_", "", Class),
    Observed = ifelse(TrueLabel == Class, 1, 0)
  )

ece_mce_result <- compute_ece_mce(df_long)
print(ece_mce_result)


## PR曲线+AUPRC
library(PRROC)
plot_pr_multi <- function(df) {
  library(PRROC)
  library(ggplot2)
  library(dplyr)
  
  classes <- c("ICM", "DCM", "CTL")
  colors <- c("ICM" = "#E6ab2a", "DCM" = "#941751", "CTL" = "#0072B2")
  pr_list <- list()
  
  for (cls in classes) {
    df_tmp <- df %>%
      mutate(Label = ifelse(TrueLabel == cls, 1, 0),
             Score = df[[paste0("prob_", cls)]])
    
    pr <- pr.curve(scores.class0 = df_tmp$Score[df_tmp$Label == 1],
                   scores.class1 = df_tmp$Score[df_tmp$Label == 0],
                   curve = TRUE)
    
    pr_list[[cls]] <- data.frame(
      Recall = pr$curve[, 1],
      Precision = pr$curve[, 2],
      Class = cls,
      AUPRC = round(pr$auc.integral, 4)
    )
  }
  
  pr_df <- bind_rows(pr_list)
  label_map <- paste0(names(pr_list), " (AUPRC=", sapply(pr_list, function(x) x$AUPRC[1]), ")")
  names(label_map) <- names(pr_list)
  
  ggplot(pr_df, aes(x = Recall, y = Precision, color = Class)) +
    geom_line(size = 2) +
    scale_color_manual(values = colors, labels = label_map) +
    labs(title = "Precision-Recall Curves", x = "Recall", y = "Precision") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 16) +  # 更大基础字体
    theme(
      axis.title = element_text(size = 42),
      axis.text = element_text(size = 42),
      axis.ticks = element_line(size = 0.8, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      panel.border = element_rect(color = "black", fill = NA, size = 2)
    )
}
plot_pr_multi(df_pred)





