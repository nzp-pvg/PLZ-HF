rm(list=ls())
gc()

library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}

resolve_project_root <- function(start = get_script_dir()) {
  cur <- normalizePath(start)
  repeat {
    if (file.exists(file.path(cur, "Revision")) && file.exists(file.path(cur, "data"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Could not locate project root from script path.")
    cur <- parent
  }
}

project_root <- resolve_project_root()
hub_core_dir <- file.path(project_root, "data", "GO", "cytoscape", "hub_core")
setwd(hub_core_dir)

icm_tbl <- readr::read_csv("ICM_hub_top60_edges_nodes.csv", show_col_types = FALSE)
dcm_tbl <- readr::read_csv("DCM_hub_top60_edges_nodes.csv", show_col_types = FALSE)



load("HF_discovery_cohort_logCPM_for_PCA.RData")  
classes_df <- classes_df %>% filter(batch != "B66")

## ====================== 2. ENSEMBL → SYMBOL ======================
map_df <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = rownames(logCPM),
                                keytype = "ENSEMBL",
                                columns = c("ENSEMBL","SYMBOL"))

map_df <- map_df %>% filter(!is.na(SYMBOL)) %>% distinct(ENSEMBL, .keep_all = TRUE)

logCPM_symbol <- logCPM[map_df$ENSEMBL, ]
rownames(logCPM_symbol) <- map_df$SYMBOL

calc_logfc <- function(subclass){
  batches <- unique(classes_df$batch)
  results <- list()
  
  for (b in batches){
    case_ids <- classes_df %>% filter(Classes == subclass, batch == b) %>% pull(ID)
    ctrl_ids <- classes_df %>% filter(Classes == "CTL", batch == b) %>% pull(ID)
    
    if(length(case_ids) > 1 && length(ctrl_ids) > 1){
      case_mean <- rowMeans(logCPM_symbol[, case_ids, drop = FALSE])
      ctrl_mean <- rowMeans(logCPM_symbol[, ctrl_ids, drop = FALSE])
      logFC <- case_mean - ctrl_mean
      results[[b]] <- logFC
    }
  }
  
  results_df <- do.call(cbind, results)
  colnames(results_df) <- paste0(names(results), "|", subclass, "_vs_CTL")
  return(results_df)
}

lfc_icm <- calc_logfc("ICM")
lfc_dcm <- calc_logfc("DCM")

icm_hubs <- icm_tbl %>%
  dplyr::arrange(desc(Degree)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(Gene, Cluster, Degree)

dcm_hubs <- dcm_tbl %>%
  dplyr::arrange(desc(Degree)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(Gene, Cluster, Degree)

icm_hubs$Disease <- "ICM"
dcm_hubs$Disease <- "DCM"

hub_nodes <- bind_rows(icm_hubs, dcm_hubs)

lfc_icm_hubs <- data.frame(Gene = hub_nodes$Gene[hub_nodes$Disease=="ICM"]) %>%
  left_join(
    lfc_icm %>% as.data.frame() %>% tibble::rownames_to_column("Gene"),
    by = "Gene"
  ) %>%
  left_join(
    lfc_dcm %>% as.data.frame() %>% tibble::rownames_to_column("Gene"),
    by = "Gene",
    suffix = c("_ICM","_DCM")
  ) %>%
  mutate(Disease = "ICM")

lfc_dcm_hubs <- data.frame(Gene = hub_nodes$Gene[hub_nodes$Disease=="DCM"]) %>%
  left_join(
    lfc_dcm %>% as.data.frame() %>% tibble::rownames_to_column("Gene"),
    by = "Gene"
  ) %>%
  left_join(
    lfc_icm %>% as.data.frame() %>% tibble::rownames_to_column("Gene"),
    by = "Gene",
    suffix = c("_DCM","_ICM")
  ) %>%
  mutate(Disease = "DCM")

lfc_all <- bind_rows(lfc_icm_hubs, lfc_dcm_hubs) %>%
  left_join(hub_nodes, by = c("Gene","Disease"))

lfc_all_expanded <- lfc_all %>%
  dplyr::mutate(Cluster_raw = Cluster) %>%
  tidyr::separate_rows(Cluster_raw, sep = "\\|") %>%
  dplyr::mutate(value = 1) %>%
  tidyr::pivot_wider(
    id_cols = c(Gene, Disease, Degree, dplyr::starts_with("B")),
    names_from = Cluster_raw,
    values_from = value,
    values_fill = list(value = 0)
  ) %>%
  dplyr::mutate(
    A = ifelse(is.na(A), 0, A),
    B = ifelse(is.na(B), 0, B),
    C = ifelse(is.na(C), 0, C),
    D = ifelse(is.na(D), 0, D),
    E = ifelse(is.na(E), 0, E),
    F = ifelse(is.na(F), 0, F)
  ) %>%
  dplyr::select(Gene, Disease, Degree, dplyr::starts_with("B"), A, B, C, D, E, F)


logfc_to_status_num <- function(x, cutoff){
  ifelse(is.na(x), NA,
         ifelse(x > cutoff, 1,
                ifelse(x < -cutoff, -1, 0)))
}

calc_consistency <- function(status_mat){
  apply(status_mat, 1, function(x){
    vals <- na.omit(x)
    if(length(vals) == 0) return(NA)
    max(table(vals)) / length(vals)
  })
}

scan_cutoffs <- function(lfc_all, cutoffs = seq(0.2, 1.0, by = 0.05)){
  results <- list()
  
  for(cut in cutoffs){
    lfc_status <- lfc_all %>%
      dplyr::mutate(across(matches("ICM|DCM"), ~logfc_to_status_num(.x, cut)))
    
    icm_mat <- lfc_status %>% dplyr::select(Gene, matches("ICM")) %>%
      tibble::column_to_rownames("Gene") %>% as.matrix()
    dcm_mat <- lfc_status %>% dplyr::select(Gene, matches("DCM")) %>%
      tibble::column_to_rownames("Gene") %>% as.matrix()
    
    icm_score <- calc_consistency(icm_mat)
    dcm_score <- calc_consistency(dcm_mat)
    
    diff_score <- abs(icm_score - dcm_score)
    
    results[[as.character(cut)]] <- data.frame(
      Gene = rownames(icm_mat),
      ICM_score = icm_score,
      DCM_score = dcm_score,
      Diff = diff_score,
      Cutoff = cut
    )
  }
  
  final_df <- bind_rows(results)
  
  summary_df <- final_df %>%
    group_by(Cutoff) %>%
    summarise(mean_diff = mean(Diff, na.rm=TRUE),
              max_diff = max(Diff, na.rm=TRUE)) %>%
    arrange(desc(mean_diff))
  
  list(detail = final_df, summary = summary_df, best = summary_df$Cutoff[1])
}

res <- scan_cutoffs(lfc_all, cutoffs = seq(0.2, 1.0, by=0.05))

print(res$best)


library(openxlsx)
library(dplyr)

logfc_to_status_num <- function(x, cutoff = 0.25){
  ifelse(is.na(x), NA,
         ifelse(x > cutoff, 1,
                ifelse(x < -cutoff, -1, 0)))
}

lfc_status_num <- lfc_all_expanded

logfc_cols <- grep("vs_CTL", colnames(lfc_status_num), value = TRUE)

for (col in logfc_cols){
  lfc_status_num[[col]] <- logfc_to_status_num(lfc_status_num[[col]], cutoff = 0.25)
}




library(dplyr)
cluster_map <- c(A=1, B=2, C=3, D=4, E=5, F=6)

df_clusters <- lfc_status_num %>%
  dplyr::select(Gene, A, B, C, D, E, F)

df_numbered <- df_clusters %>%
  dplyr::mutate(across(A:F, ~ ifelse(. == 1, cluster_map[cur_column()], 0)))

compact_row <- function(x){
  vals <- x[x != 0]
  c(vals, rep(0, length(x) - length(vals)))
}

df_compact <- df_numbered
df_compact[,2:7] <- t(apply(df_numbered[,2:7], 1, compact_row))

colnames(df_compact)[2:7] <- paste0("Cluster_", 1:6)

print(df_compact)

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "logFC_values")
writeData(wb, "logFC_values", lfc_all_expanded)

addWorksheet(wb, "UpDownNum")
writeData(wb, "UpDownNum", lfc_status_num)

addWorksheet(wb, "ClusterCompact")
writeData(wb, "ClusterCompact", df_compact)

saveWorkbook(wb, "HubGenes_logFC_and_StatusNum.xlsx", overwrite = TRUE)

message("Exported: HubGenes_logFC_and_StatusNum.xlsx (3 worksheets)")
