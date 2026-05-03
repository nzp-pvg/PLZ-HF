#!/usr/bin/env Rscript
# Merged and organized from:
#   - Raw_Data_Processing.R
#   - process_GSE_dataset.R
#   - counts_processing.R
#   - CVD_classification.R

suppressPackageStartupMessages({
  library(dplyr)
  library(GEOquery)
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
  cands <- unique(c(
    if (!is.na(from_arg)) from_arg else NULL,
    if (nzchar(from_env)) from_env else NULL,
    getwd()
  ))
  for (seed in cands) {
    cur <- normalizePath(seed, winslash = "/", mustWork = FALSE)
    for (i in 0:6) {
      probe <- normalizePath(file.path(cur, paste(rep("..", i), collapse = "/")), winslash = "/", mustWork = FALSE)
      if (dir.exists(file.path(probe, "data", "PLZ-HF-main", "Script"))) return(probe)
    }
  }
  stop("Unable to resolve project root. Use --project-root=/path/to/CVD_MS_1")
}

standardize_group <- function(x) {
  x <- tolower(trimws(x))
  ifelse(grepl("non-failing|healthy|control|not failing", x), "CTL",
  ifelse(grepl("ischemic", x), "ICM",
  ifelse(grepl("dilated", x), "DCM",
  ifelse(grepl("hypertrophic", x), "HCM",
  ifelse(grepl("restrictive", x), "RCM", NA_character_)))))
}

infer_batch_from_gse <- function(gse_id) {
  digits <- sub("^GSE", "", gse_id)
  paste0("B", substr(digits, nchar(digits) - 1L, nchar(digits)))
}

inspect_geo_series <- function(gse_id, destdir = ".") {
  gset <- GEOquery::getGEO(gse_id, destdir = destdir, AnnotGPL = FALSE, getGPL = FALSE)
  list(
    expr = Biobase::exprs(gset[[1]]),
    pdata = Biobase::pData(gset[[1]])
  )
}

fetch_geo_sample_info <- function(gse_id) {
  gse_obj <- GEOquery::getGEO(gse_id, GSEMatrix = FALSE)
  gsm_list <- GEOquery::GSMList(gse_obj)
  do.call(rbind, lapply(gsm_list, function(gsm) {
    meta <- Meta(gsm)
    disease_line <- grep("disease", meta$characteristics_ch1, value = TRUE)
    disease <- if (length(disease_line) > 0) sub(".*disease: ?", "", disease_line) else NA_character_
    data.frame(
      sample_id = meta$geo_accession,
      title = meta$title,
      group = standardize_group(disease),
      stringsAsFactors = FALSE
    )
  }))
}

align_counts_with_metadata <- function(expr_matrix, sample_info, batch) {
  sample_info <- sample_info[!is.na(sample_info$group), , drop = FALSE]
  common_ids <- intersect(sample_info$sample_id, colnames(expr_matrix))
  sample_info <- sample_info[sample_info$sample_id %in% common_ids, , drop = FALSE]
  expr_matrix <- expr_matrix[, common_ids, drop = FALSE]
  expr_t <- as.data.frame(t(expr_matrix))
  expr_t$sample_id <- rownames(expr_t)
  expr_t <- expr_t[match(sample_info$sample_id, expr_t$sample_id), , drop = FALSE]
  stopifnot(all(expr_t$sample_id == sample_info$sample_id))
  merged <- cbind(sample_info, expr_t[, -ncol(expr_t), drop = FALSE])
  merged$batch <- batch
  merged[, c("sample_id", "title", "group", "batch", setdiff(colnames(merged), c("sample_id", "title", "group", "batch")))]
}

process_public_gse_dataset <- function(gse_id, tsv_file) {
  expr_matrix <- read.delim(tsv_file, row.names = 1, check.names = FALSE)
  sample_info <- fetch_geo_sample_info(gse_id)
  align_counts_with_metadata(expr_matrix, sample_info, infer_batch_from_gse(gse_id))
}

process_gse141910 <- function(tsv_file, sra_csv) {
  expr_matrix <- read.delim(tsv_file, row.names = 1, check.names = FALSE)
  sra <- read.csv(sra_csv, stringsAsFactors = FALSE)
  sample_info <- data.frame(
    sample_id = sra$GEO_Accession..exp.,
    title = sra$Sample.Name,
    group = standardize_group(sra$disease_state),
    stringsAsFactors = FALSE
  )
  align_counts_with_metadata(expr_matrix, sample_info, "B10")
}

process_gse55296 <- function(tsv_file, sra_csv) {
  expr_matrix <- read.delim(tsv_file, row.names = 1, check.names = FALSE)
  sra <- read.csv(sra_csv, stringsAsFactors = FALSE)
  sample_info <- data.frame(
    sample_id = sra$GEO_Accession..exp.,
    title = sra$Sample.Name,
    group = standardize_group(sra$disease_state),
    stringsAsFactors = FALSE
  )
  align_counts_with_metadata(expr_matrix, sample_info, "B96")
}

merge_processed_matrices <- function(matrix_list) {
  cleaned <- lapply(matrix_list, function(df) df[complete.cases(df), , drop = FALSE])
  merged <- do.call(rbind, cleaned)
  rownames(merged) <- merged$sample_id
  merged
}

build_counts_objects <- function(merged_matrix) {
  sample_info <- merged_matrix[, c("sample_id", "title", "group", "batch"), drop = FALSE]
  counts <- merged_matrix[, setdiff(colnames(merged_matrix), c("sample_id", "title", "group", "batch")), drop = FALSE]
  counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
  list(
    exp_no_0 = counts,
    sample_info = sample_info,
    count = counts,
    samples = sample_info
  )
}

main <- function() {
  project_root <- resolve_project_root()
  raw_dir <- get_arg_value("raw-dir", file.path(project_root, "data", "RNA-Seq", "raw_counts"))
  out_dir <- get_arg_value("out-dir", file.path(project_root, "data", "PLZ-HF-main", "Script", "organized_codes", "outputs"))
  mode <- get_arg_value("mode", "merge")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (mode == "public-gse") {
    gse_id <- get_arg_value("gse-id", NULL)
    tsv_file <- get_arg_value("tsv-file", NULL)
    stopifnot(!is.null(gse_id), !is.null(tsv_file))
    merged <- process_public_gse_dataset(gse_id, tsv_file)
    saveRDS(merged, file.path(out_dir, paste0(gse_id, "_matrix.rds")))
  } else if (mode == "gse141910") {
    merged <- process_gse141910(
      tsv_file = get_arg_value("tsv-file", file.path(raw_dir, "GSE141910.tsv")),
      sra_csv = get_arg_value("sra-file", file.path(raw_dir, "SraRunTable_GSE141910.csv"))
    )
    saveRDS(merged, file.path(out_dir, "GSE141910_matrix.rds"))
  } else if (mode == "gse55296") {
    merged <- process_gse55296(
      tsv_file = get_arg_value("tsv-file", file.path(raw_dir, "GSE55296.tsv")),
      sra_csv = get_arg_value("sra-file", file.path(raw_dir, "SraRunTable_GSE55296.csv"))
    )
    saveRDS(merged, file.path(out_dir, "GSE55296_matrix.rds"))
  } else if (mode == "merge") {
    rds_dir <- get_arg_value("rds-dir", out_dir)
    rds_files <- list.files(rds_dir, pattern = "_matrix\\.rds$", full.names = TRUE)
    matrix_list <- lapply(rds_files, readRDS)
    merged <- merge_processed_matrices(matrix_list)
    objs <- build_counts_objects(merged)
    saveRDS(merged, file.path(out_dir, "All_CVD_samples_merged_matrix.rds"))
    save(objs$exp_no_0, objs$sample_info, file = file.path(out_dir, "GSE_processed_counts.RData"))
    save(objs$count, objs$samples, file = file.path(out_dir, "HF_processed_counts.RData"))
  } else {
    stop("Unknown mode: ", mode)
  }
}

if (sys.nframe() == 0L) {
  main()
}
