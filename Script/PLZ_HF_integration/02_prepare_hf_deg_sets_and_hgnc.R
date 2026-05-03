# Task: Prepare HF DEG sets, apply nFreq thresholds, and convert Ensembl IDs to HGNC symbols. Corresponds to the HF DEG preparation stage used before PLZ-HF integration.
# Manuscript mapping: PLZ-HF integration module; see 00_CODE_DIRECTORY.md for task scope.

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
go_s4_dir <- file.path(project_root, "data", "GO", "GO_s4")
hub_core_dir <- file.path(project_root, "data", "GO", "cytoscape", "hub_core")

## Option 1
setwd(go_s4_dir)

# load HF subtyping classification result
res <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")
str(res, max.level = 2)
names(res)
res[[1]]

# extrat DEGs
get_model_features <- function(res, subclass) {
  unlist(lapply(res, function(rep) {
    lapply(rep[[subclass]], function(fold) fold$DEGs)
  }), use.names = FALSE)
}

DEGs_ICM <- get_model_features(res, "ICM")
DEGs_DCM <- get_model_features(res, "DCM")
DEGs_CTL <- get_model_features(res, "CTL")

DEG_ICM <- names(which(table(DEGs_ICM) >= 17))
DEG_DCM <- names(which(table(DEGs_DCM) >= 8))
DEG_CTL <- names(which(table(DEGs_CTL) >= 2))

HF_spec <- setdiff(union(DEG_ICM, DEG_DCM), DEG_CTL)
HF_core <- intersect(DEG_ICM, DEG_DCM)
ICM <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))
ICM_DCM_diff <- setdiff(union(DEG_ICM, DEG_DCM), intersect(DEG_ICM, DEG_DCM))

length(intersect(DCM_hgnc, ICM_hgnc)) 






HF_gene_sets <- list(
  HF_spec = HF_spec,
  ICM = ICM,
  DCM = DCM,
  ICM_DCM_diff = ICM_DCM_diff,
  HF_core = HF_core,
  CTL = CTL
)

sapply(HF_gene_sets, length)

# save
save(CTL, DCM, ICM, HF_spec, HF_core, ICM_DCM_diff, HF_gene_sets, file="HF_DEGs_top560.RData")


## Gene ID link
library(biomaRt)
#load("HF_DEGs_top560.RData")

#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  host = "https://useast.ensembl.org"
)
deg_ids <- unique(c(CTL, DCM, ICM, HF_spec, HF_core, ICM_DCM_diff))

id_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = deg_ids,
  mart = ensembl
)

head(id_mapping)

convert_deg <- function(deg_list, mapping_table) {
  mapping_table %>%
    filter(ensembl_gene_id %in% deg_list, hgnc_symbol != "") %>%
    pull(hgnc_symbol) %>%
    unique()
}

ICM_hgnc <- convert_deg(ICM, id_mapping)
