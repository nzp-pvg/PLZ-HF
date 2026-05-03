#!/usr/bin/env Rscript
# Merged and organized from:
#   - PLZ_Exposure_evidence_tiering_37.R
#   - PLZ_cluster_k_by_k_compare.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(purrr)
  library(cluster)
  library(ggplot2)
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
  for (seed in cands) {
    cur <- normalizePath(seed, winslash = "/", mustWork = FALSE)
    for (i in 0:6) {
      probe <- normalizePath(file.path(cur, paste(rep("..", i), collapse = "/")), winslash = "/", mustWork = FALSE)
      if (dir.exists(file.path(probe, "data"))) return(probe)
    }
  }
  stop("Unable to resolve project root. Use --project-root=/path/to/CVD_MS_1")
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

build_exposure_tiering_table <- function(project_root, fcc_file = NULL, out_dir) {
  if (is.null(fcc_file)) {
    fcc_file <- first_existing(c(
      file.path(project_root, "data", "Toxicity", "FCC.xlsx"),
      file.path(project_root, "supp", "FCC.xlsx"),
      file.path(project_root, "data", "Zenodo", "Plasticizer", "FCC.xlsx")
    ))
  }
  if (is.na(fcc_file) || !file.exists(fcc_file)) stop("FCC file not found.")

  ps <- openxlsx::read.xlsx(fcc_file, sheet = "Priority_Substances", check.names = TRUE)
  tr <- openxlsx::read.xlsx(fcc_file, sheet = "Tier", check.names = TRUE)
  ps$Abbr <- ps[[2]]
  ps37 <- ps %>% filter(!is.na(Abbr), Abbr != "") %>% slice(1:37)
  out <- ps37 %>% left_join(tr, by = "Abbr")
  write.csv(out, file.path(out_dir, "plasticizer_exposure_tiering_37.csv"), row.names = FALSE)
  out
}

build_overlap_distance <- function(gene_lists) {
  chems <- names(gene_lists)
  m <- matrix(0, nrow = length(chems), ncol = length(chems), dimnames = list(chems, chems))
  for (i in seq_along(chems)) {
    for (j in seq_along(chems)) {
      a <- unique(gene_lists[[i]])
      b <- unique(gene_lists[[j]])
      jac <- if (length(union(a, b)) == 0) 0 else length(intersect(a, b)) / length(union(a, b))
      m[i, j] <- 1 - jac
    }
  }
  as.dist(m)
}

run_cluster_k_sensitivity <- function(gene_lists, k_range = 3:9, n_boot = 200L, seed = 20260219L) {
  set.seed(seed)
  dist_mat <- build_overlap_distance(gene_lists)
  rows <- list()
  idx <- 1
  for (k in k_range) {
    hc <- hclust(dist_mat, method = "average")
    clusters <- cutree(hc, k = k)
    sil <- cluster::silhouette(clusters, dist_mat)
    rows[[idx]] <- data.frame(k = k, silhouette_mean = mean(sil[, "sil_width"], na.rm = TRUE))
    idx <- idx + 1
  }
  bind_rows(rows)
}

main <- function() {
  project_root <- resolve_project_root()
  out_dir <- get_arg_value("out-dir", file.path(project_root, "data", "PLZ-HF-main", "Script", "organized_codes", "outputs"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  mode <- get_arg_value("mode", "exposure")

  if (mode == "exposure") {
    build_exposure_tiering_table(project_root, get_arg_value("fcc-file", NULL), out_dir)
  } else if (mode == "cluster-k") {
    gene_file <- get_arg_value("gene-file", NULL)
    stopifnot(!is.null(gene_file), file.exists(gene_file))
    gene_lists <- readRDS(gene_file)
    write.csv(run_cluster_k_sensitivity(gene_lists), file.path(out_dir, "cluster_k_sensitivity.csv"), row.names = FALSE)
  } else {
    stop("Unknown mode: ", mode)
  }
}

if (sys.nframe() == 0L) {
  main()
}
