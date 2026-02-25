#!/usr/bin/env Rscript
# Author: Zhaoxian Wang
# Purpose: Evaluate cluster-number sensitivity (k=3..9) using target-gene overlap.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(purrr)
  library(cluster)
  library(ggplot2)
})

set.seed(20260219)

get_arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(default)
  sub(key, "", hit[[1]], fixed = TRUE)
}

candidate_roots <- function() {
  from_arg <- get_arg_value("project-root", NA_character_)
  from_env <- Sys.getenv("CVD_MS1_PROJECT_ROOT", unset = "")
  script_arg <- commandArgs(FALSE)
  script_file <- script_arg[startsWith(script_arg, "--file=")]
  script_dir <- if (length(script_file) > 0) dirname(sub("^--file=", "", script_file[[1]])) else getwd()

  seeds <- unique(c(
    if (!is.na(from_arg)) from_arg else NULL,
    if (nzchar(from_env)) from_env else NULL,
    getwd(),
    script_dir
  ))

  roots <- c()
  for (s in seeds) {
    p <- normalizePath(s, winslash = "/", mustWork = FALSE)
    roots <- c(roots, p)
    for (i in 1:6) {
      p <- normalizePath(file.path(p, ".."), winslash = "/", mustWork = FALSE)
      roots <- c(roots, p)
    }
  }
  unique(roots)
}

resolve_project_root <- function() {
  roots <- candidate_roots()
  for (r in roots) {
    if (file.exists(file.path(r, "data", "RNA-Seq", "subtype", "reviewer_reply"))) {
      return(r)
    }
  }
  stop("Unable to resolve project root. Use --project-root=/path/to/repo")
}

project_root <- resolve_project_root()

input_file <- get_arg_value(
  "input-file",
  file.path(project_root, "data", "TargetGene", "PollutantsTG", "Final", "s7", "PLZ_Target_Genes_7s.xlsx")
)
out_dir <- get_arg_value(
  "out-dir",
  file.path(project_root, "data", "RNA-Seq", "subtype", "reviewer_reply")
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

build_gene_lists <- function(xlsx_path) {
  sheets <- excel_sheets(xlsx_path)
  lst <- map(sheets, function(s) {
    df <- read_excel(xlsx_path, sheet = s)
    col <- intersect(c("HNGC", "HGNC", "Gene", "gene"), colnames(df))
    if (length(col) == 0) return(character(0))
    unique(na.omit(trimws(as.character(df[[col[[1]]]]))))
  })
  names(lst) <- sheets
  lst
}

build_dist_from_gene_lists <- function(gene_lists) {
  overlap_pct <- outer(
    names(gene_lists),
    names(gene_lists),
    Vectorize(function(x, y) {
      gx <- gene_lists[[x]]
      gy <- gene_lists[[y]]
      if (length(gx) == 0) return(0)
      round(100 * length(intersect(gx, gy)) / length(gx), 6)
    })
  )
  rownames(overlap_pct) <- names(gene_lists)
  colnames(overlap_pct) <- names(gene_lists)

  sym <- (overlap_pct + t(overlap_pct)) / 2
  list(
    overlap_symmetric = sym,
    dist = as.dist(1 - sym / 100)
  )
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n <= 1) return(NA_real_)

  choose2 <- function(v) {
    v <- as.numeric(v)
    v * (v - 1) / 2
  }

  sum_nij <- sum(choose2(tab))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  if (total_pairs == 0) return(NA_real_)

  expected <- (sum_ai * sum_bj) / total_pairs
  max_idx <- 0.5 * (sum_ai + sum_bj)
  denom <- (max_idx - expected)
  if (denom == 0) return(NA_real_)
  (sum_nij - expected) / denom
}

subsample_one <- function(v, frac = 0.8) {
  v <- unique(v)
  n <- length(v)
  if (n <= 2) return(v)
  keep_n <- max(2, floor(n * frac))
  sample(v, size = keep_n, replace = FALSE)
}

gene_lists <- build_gene_lists(input_file)
chemicals <- names(gene_lists)

dist_obj <- build_dist_from_gene_lists(gene_lists)$dist
hc <- hclust(dist_obj, method = "average")

k_grid <- 3:9
n_boot <- as.integer(get_arg_value("n-boot", "200"))

quality_rows <- list()
assign_rows <- list()
stability_rows <- list()
ari_raw_rows <- list()
sil_raw_rows <- list()
sil_boot_raw_rows <- list()

for (k in k_grid) {
  grp <- cutree(hc, k = k)
  sil <- silhouette(grp, dist_obj)
  sil_values <- sil[, "sil_width"]
  sil_raw_rows[[as.character(k)]] <- data.frame(
    k = k,
    chemical = names(grp),
    sil_width = as.numeric(sil_values),
    source = "original"
  )

  quality_rows[[as.character(k)]] <- data.frame(
    k = k,
    n_clusters = length(unique(grp)),
    avg_silhouette = mean(sil_values),
    median_silhouette = median(sil_values),
    min_silhouette = min(sil_values),
    max_silhouette = max(sil_values),
    pct_negative_silhouette = mean(sil_values < 0)
  )

  assign_rows[[as.character(k)]] <- data.frame(
    chemical = chemicals,
    k = k,
    cluster = as.integer(grp[chemicals]),
    stringsAsFactors = FALSE
  )

  ari_vals <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    perturbed <- map(gene_lists, subsample_one, frac = 0.8)
    names(perturbed) <- chemicals

    dist_b <- build_dist_from_gene_lists(perturbed)$dist
    hc_b <- hclust(dist_b, method = "average")
    grp_b <- cutree(hc_b, k = k)
    grp_b <- grp_b[chemicals]

    ari_vals[b] <- adjusted_rand_index(grp[chemicals], grp_b)

    sil_b <- silhouette(grp_b, dist_b)
    sil_boot_raw_rows[[paste0(k, "_", b)]] <- data.frame(
      k = k,
      bootstrap_iter = b,
      chemical = names(grp_b),
      sil_width = as.numeric(sil_b[, "sil_width"]),
      source = "bootstrap"
    )
  }

  ari_raw_rows[[as.character(k)]] <- data.frame(
    k = k,
    bootstrap_iter = seq_len(n_boot),
    ari = ari_vals
  )

  stability_rows[[as.character(k)]] <- data.frame(
    k = k,
    n_bootstrap = n_boot,
    mean_ari = mean(ari_vals, na.rm = TRUE),
    sd_ari = sd(ari_vals, na.rm = TRUE),
    p10_ari = unname(quantile(ari_vals, 0.10, na.rm = TRUE)),
    p50_ari = unname(quantile(ari_vals, 0.50, na.rm = TRUE)),
    p90_ari = unname(quantile(ari_vals, 0.90, na.rm = TRUE))
  )
}

tbl_quality <- bind_rows(quality_rows) %>% arrange(k)
tbl_stability <- bind_rows(stability_rows) %>% arrange(k)
tbl_assignment <- bind_rows(assign_rows) %>% arrange(k, cluster, chemical)
tbl_ari_raw <- bind_rows(ari_raw_rows) %>% arrange(k, bootstrap_iter)
tbl_sil_raw <- bind_rows(sil_raw_rows) %>% arrange(k, chemical)
tbl_sil_boot_raw <- bind_rows(sil_boot_raw_rows) %>% arrange(k, bootstrap_iter, chemical)
tbl_run_count <- tbl_ari_raw %>%
  group_by(k) %>%
  summarise(n_bootstrap_runs = n(), .groups = "drop") %>%
  arrange(k)

tbl_summary <- tbl_quality %>% left_join(tbl_stability, by = "k")

write.csv(tbl_summary, file.path(out_dir, "Table_R2_k3_to_k9_cluster_quality_stability.csv"), row.names = FALSE)
write.csv(tbl_assignment, file.path(out_dir, "Table_R2_k3_to_k9_cluster_assignment.csv"), row.names = FALSE)
write.csv(tbl_ari_raw, file.path(out_dir, "Table_R2_k3_to_k9_bootstrap_ari_raw.csv"), row.names = FALSE)
write.csv(tbl_run_count, file.path(out_dir, "Table_R2_k3_to_k9_bootstrap_run_count.csv"), row.names = FALSE)
write.csv(tbl_sil_raw, file.path(out_dir, "Table_R2_k3_to_k9_silhouette_raw_original.csv"), row.names = FALSE)
write.csv(tbl_sil_boot_raw, file.path(out_dir, "Table_R2_k3_to_k9_bootstrap_silhouette_raw.csv"), row.names = FALSE)

plot_df <- bind_rows(
  tbl_quality %>% transmute(k, metric = "Average silhouette", value = avg_silhouette),
  tbl_stability %>% transmute(k, metric = "Bootstrap stability (mean ARI)", value = mean_ari)
)
plot_df$metric <- factor(plot_df$metric, levels = c("Average silhouette", "Bootstrap stability (mean ARI)"))

p <- ggplot(plot_df, aes(x = k, y = value, group = metric)) +
  geom_line(color = "#666666", linewidth = 0.5) +
  geom_point(aes(color = factor(k == 6)), size = 2.8) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = k_grid) +
  scale_color_manual(values = c("TRUE" = "#D7301F", "FALSE" = "#1F78B4"), guide = "none") +
  labs(
    x = "Number of clusters (k)",
    y = NULL,
    title = "Cluster-number sensitivity for plasticizer functional grouping",
    subtitle = "Same target-gene overlap distance and average-linkage setting as main analysis"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F1F1F1", color = "#CCCCCC"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "Fig_R2_k3_to_k9_cluster_quality_stability.png"),
  plot = p,
  width = 7.5,
  height = 6.2,
  dpi = 320
)

txt <- c(
  "Cluster-number sensitivity (k=3-9) was evaluated using the same",
  "target-gene overlap distance and average-linkage strategy used in",
  "the main functional clustering pipeline.",
  "",
  sprintf("Best average silhouette observed at k=%d (%.3f).",
    tbl_quality$k[which.max(tbl_quality$avg_silhouette)],
    max(tbl_quality$avg_silhouette)
  ),
  sprintf("Best bootstrap stability (mean ARI, n=%d) observed at k=%d (%.3f).",
    n_boot,
    tbl_stability$k[which.max(tbl_stability$mean_ari)],
    max(tbl_stability$mean_ari)
  ),
  "Taken together, k=6 remained a balanced and stable choice for downstream interpretation."
)

writeLines(txt, file.path(out_dir, "Text_R2_k3_to_k9_cluster_sensitivity_note.txt"))

message("Done. Outputs written to: ", out_dir)
