#!/usr/bin/env Rscript
# Manuscript outputs:
#   - Fig. S3A (publication plot for permutation negative control)
# Source basis:
#   - Revision/R2/reviewer_reply/plot_r3f_permutation_smooth.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

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
      if (dir.exists(file.path(probe, "Revision", "R2", "reviewer_reply"))) return(probe)
    }
  }
  stop("Unable to resolve project root.")
}

col_observed <- "#7F7F7F"
col_permuted <- "#C34062"
alpha_pts <- 0.5
alpha_observed <- 0.8
col_q95 <- "#F5BD4D"
perm_fill <- grDevices::adjustcolor(col_permuted, alpha.f = 0.50)
perm_outline <- grDevices::adjustcolor(col_permuted, alpha.f = 0.48)

project_root <- resolve_project_root()
base_dir <- file.path(project_root, "Revision", "R2", "reviewer_reply")
src_dir <- base_dir

repeat_file <- file.path(src_dir, "Table_R3F_permutation_repeat_metrics.csv")
summary_file <- file.path(src_dir, "Table_R3F_permutation_summary.csv")

rep_df <- read.csv(repeat_file, stringsAsFactors = FALSE)
sum_df <- read.csv(summary_file, stringsAsFactors = FALSE)

macro_rep <- bind_rows(
  rep_df %>%
    filter(LabelMode == "observed") %>%
    distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
    transmute(Scenario = "Observed", Repeat, Metric = "AUROC", Value = MacroAUROC),
  rep_df %>%
    filter(LabelMode == "permuted") %>%
    distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
    transmute(Scenario = "Permuted (N=500)", Repeat, Metric = "AUROC", Value = MacroAUROC),
  rep_df %>%
    filter(LabelMode == "observed") %>%
    distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
    transmute(Scenario = "Observed", Repeat, Metric = "AUPRC", Value = MacroAUPRC),
  rep_df %>%
    filter(LabelMode == "permuted") %>%
    distinct(Repeat, MacroAUROC, MacroAUPRC) %>%
    transmute(Scenario = "Permuted (N=500)", Repeat, Metric = "AUPRC", Value = MacroAUPRC)
)

macro_row <- sum_df %>% filter(Subclass == "MacroMean")

metric_summary <- data.frame(
  Metric = c("AUROC", "AUPRC"),
  Observed = c(macro_row$observed_mean_AUROC, macro_row$observed_mean_AUPRC),
  Q95 = c(macro_row$perm_q95_AUROC, macro_row$perm_q95_AUPRC),
  P = c(macro_row$empirical_p_AUROC, macro_row$empirical_p_AUPRC),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    Q025 = quantile(macro_rep$Value[macro_rep$Metric == Metric & macro_rep$Scenario == "Permuted (N=500)"], 0.025),
    Q25 = quantile(macro_rep$Value[macro_rep$Metric == Metric & macro_rep$Scenario == "Permuted (N=500)"], 0.25),
    Q50 = quantile(macro_rep$Value[macro_rep$Metric == Metric & macro_rep$Scenario == "Permuted (N=500)"], 0.50),
    Q75 = quantile(macro_rep$Value[macro_rep$Metric == Metric & macro_rep$Scenario == "Permuted (N=500)"], 0.75),
    Q975 = quantile(macro_rep$Value[macro_rep$Metric == Metric & macro_rep$Scenario == "Permuted (N=500)"], 0.975)
  ) %>%
  ungroup()

build_panel <- function(metric_name, y_limits, q95_label_y_offset, labeled = FALSE) {
  panel_df <- macro_rep %>% filter(Metric == metric_name)
  sum_row <- metric_summary %>% filter(Metric == metric_name)

  perm_box <- data.frame(
    Scenario = "Permuted (N=500)",
    ymin = sum_row$Q025,
    lower = sum_row$Q25,
    middle = sum_row$Q50,
    upper = sum_row$Q75,
    ymax = sum_row$Q975
  )

  observed_df <- data.frame(
    Scenario = "Observed",
    Value = sum_row$Observed
  )

  x_q95 <- if (labeled) 1.62 else 1.62
  x_p <- if (labeled) 1.97 else 1.97
  p_label <- paste0("empirical p=", formatC(sum_row$P, format = "f", digits = 3))

  ggplot() +
    geom_hline(
      yintercept = sum_row$Q95,
      linetype = "dashed",
      color = col_q95,
      linewidth = 0.45
    ) +
    geom_jitter(
      data = panel_df %>% filter(Scenario == "Permuted (N=500)"),
      aes(x = Scenario, y = Value),
      width = 0.09,
      height = 0,
      size = 2.4,
      shape = 21,
      stroke = 0.28,
      alpha = alpha_pts,
      fill = perm_fill,
      color = perm_outline
    ) +
    geom_segment(
      data = observed_df,
      aes(x = 0.76, xend = 1.24, y = Value, yend = Value),
      linewidth = 2.0,
      color = col_observed,
      alpha = alpha_observed,
      lineend = "butt"
    ) +
    geom_point(
      data = observed_df,
      aes(x = Scenario, y = Value),
      size = 2.0,
      shape = 16,
      color = "black",
      alpha = 1
    ) +
    geom_boxplot(
      data = perm_box,
      aes(x = Scenario, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity",
      width = 0.28,
      fill = NA,
      color = "black",
      linewidth = 0.375
    ) +
    annotate(
      "text",
      x = x_q95,
      y = sum_row$Q95 + q95_label_y_offset,
      label = "q95",
      hjust = 0,
      vjust = -0.2,
      family = "Arial",
      size = if (labeled) 3.0 else 0
    ) +
    annotate(
      "text",
      x = x_p,
      y = y_limits[2] - 0.03,
      label = p_label,
      hjust = 1,
      vjust = 1,
      family = "Arial",
      size = if (labeled) 3.0 else 0
    ) +
    scale_x_discrete(
      limits = c("Observed", "Permuted (N=500)"),
      labels = if (labeled) c("Observed", "Permuted\n(n=500)") else c("", "")
    ) +
    scale_y_continuous(limits = c(0.2, 1.0), expand = c(0, 0)) +
    labs(
      title = if (labeled) metric_name else NULL,
      x = if (labeled) NULL else NULL,
      y = if (labeled) "Overall (macro-mean)" else NULL
    ) +
    theme_bw(base_family = "Arial", base_size = 10) +
    theme(
      panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.7),
      panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.7),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(color = "black", size = 18),
      axis.text.x = if (labeled) element_text(color = "black", size = 9) else element_blank(),
      axis.title.y = if (labeled) element_text(color = "black", size = 10, margin = margin(r = 6)) else element_blank(),
      axis.title.x = element_blank(),
      plot.title = if (labeled) element_text(color = "black", size = 11, hjust = 0.5, margin = margin(b = 4)) else element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.275),
      axis.line.x = if (labeled) element_line(color = "black", linewidth = 0.275) else element_blank(),
      axis.ticks.y = element_line(color = "black", linewidth = 0.275),
      axis.ticks.x = element_line(color = "black", linewidth = 0.275),
      axis.ticks.length.y = unit(0.12, "cm"),
      axis.ticks.length.x = unit(0.12, "cm"),
      plot.margin = margin(4, 6, 4, 4, "mm")
    )
}

p_auroc <- build_panel("AUROC", c(0.2, 1.00), 0.012)
p_auprc <- build_panel("AUPRC", c(0.2, 1.00), 0.018)
p_auroc_labeled <- build_panel("AUROC", c(0.2, 1.00), 0.012, labeled = TRUE)
p_auprc_labeled <- build_panel("AUPRC", c(0.2, 1.00), 0.018, labeled = TRUE)

p_all <- p_auroc + p_auprc + plot_layout(ncol = 2)
p_all_labeled <- p_auroc_labeled + p_auprc_labeled + plot_layout(ncol = 2)

out_png <- file.path(base_dir, "Fig_R3F_permutation_smooth_violin.png")
out_pdf <- file.path(base_dir, "Fig_R3F_permutation_smooth_violin.pdf")
out_png_labeled <- file.path(base_dir, "Fig_R3F_permutation_smooth_violin_labeled.png")
out_pdf_labeled <- file.path(base_dir, "Fig_R3F_permutation_smooth_violin_labeled.pdf")

ggsave(out_png, p_all, width = 5.28, height = 3.5, dpi = 600, bg = "white")
ggsave(out_pdf, p_all, width = 5.28, height = 3.5, device = cairo_pdf, bg = "white")
ggsave(out_png_labeled, p_all_labeled, width = 5.52, height = 4.1, dpi = 600, bg = "white")
ggsave(out_pdf_labeled, p_all_labeled, width = 5.52, height = 4.1, device = cairo_pdf, bg = "white")

cat(out_png, "\n", out_pdf, "\n", out_png_labeled, "\n", out_pdf_labeled, "\n")
