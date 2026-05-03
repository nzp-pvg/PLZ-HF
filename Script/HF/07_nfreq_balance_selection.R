#!/usr/bin/env Rscript
# Author: Zhaoxian Wang


suppressPackageStartupMessages({
  library(dplyr)
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
      if (file.exists(file.path(r, "data", "RNA-Seq", "subtype", "nFreq"))) return(r)
    }
  }
  stop("Unable to resolve project root. Use --project-root=/path/to/repo")
}

project_root <- resolve_project_root()
res_dir <- file.path(project_root, "data", "RNA-Seq", "subtype", "nFreq", "nFreq_results")
out_curve <- file.path(project_root, "data", "RNA-Seq", "subtype", "reviewer_reply", "Table_R3D_nFreq_sensitivity_curve.csv")
out_plateau <- file.path(project_root, "data", "RNA-Seq", "subtype", "reviewer_reply", "Table_R3D_nFreq_plateau_summary.csv")
dir.create(dirname(out_curve), recursive = TRUE, showWarnings = FALSE)

files <- list.files(res_dir, pattern = "^nFreq_[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No nFreq_*.rds found in: ", res_dir)

rows <- list()
k <- 1
for (f in files) {
  base <- basename(f)
  n_val <- as.integer(sub("^nFreq_([0-9]+)\\.rds$", "\\1", base))
  x <- readRDS(f)
  for (sc in c("ICM", "DCM", "CTL")) {
    rows[[k]] <- data.frame(
      Subclass = sc,
      nFreq = n_val,
      mean_AUROC = x[[sc]]$meanAUC,
      sd_AUROC = x[[sc]]$sdAUC,
      nGenes = x[[sc]]$nGenes,
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }
}

curve <- bind_rows(rows) %>% arrange(Subclass, nFreq)
write.csv(curve, out_curve, row.names = FALSE)

selected <- c("ICM" = 17L, "DCM" = 8L, "CTL" = 2L)
plateau <- curve %>%
  group_by(Subclass) %>%
  summarise(
    best_nFreq = nFreq[which.max(mean_AUROC)][1],
    best_AUROC = max(mean_AUROC),
    selected_nFreq = unname(selected[first(Subclass)]),
    selected_AUROC = mean_AUROC[nFreq == selected_nFreq][1],
    delta_selected_vs_best = selected_AUROC - best_AUROC,
    plateau_min_nFreq = min(nFreq[mean_AUROC >= (best_AUROC - 0.01)]),
    plateau_max_nFreq = max(nFreq[mean_AUROC >= (best_AUROC - 0.01)]),
    selection_basis = "Visual balance-point selection from Fig. B (AUC plateau vs DEG count trade-off): CTL=2, DCM=8, ICM=17",
    delta_method = "direct from scanned nFreq grid",
    .groups = "drop"
  )
write.csv(plateau, out_plateau, row.names = FALSE)

cat("Rebuilt R3D files:\n")
cat(" -", out_curve, "\n")
cat(" -", out_plateau, "\n")
cat("\nCheck selected rows:\n")
print(plateau)
