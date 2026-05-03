# Task: Generate UpSet-style overlap summaries for HF DEG sets and plasticizer cluster gene sets.
# Manuscript mapping: PLZ-HF integration module; see 00_CODE_DIRECTORY.md for task scope.

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
setwd(go_s4_dir)

##### draw Interscetion map using UpSetR

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
library(UpSetR)
library(ComplexHeatmap)
library(reshape2)
library(patchwork) 

setwd(go_s4_dir)


load("HF_DEGs_hgnc_top560.RData")
load("HF_Plasticizer_Intersection_7S_Data.RData")

deg_sets <- list(
  DCM_only   = DCM_hgnc,
  ICM_only   = ICM_hgnc,
  CTL_only   = CTL_hgnc,
  HF_specific    = HF_spec_hgnc,
  HF_core    = HF_core_hgnc,
  ICM_DCM_diff. = ICM_DCM_diff_hgnc
)

input <- fromList(deg_sets)

venn.plot <- venn.diagram(
  x = list(
    ICM = ICM_hgnc,
    DCM = DCM_hgnc,
    CTL = CTL_hgnc
  ),
  filename = NULL,
  fill = c("#C34062", "#F5BD4D", "#005493"),
  alpha = 0.5,
  cex = 2.5,
  fontface = "bold",
  cat.cex = 2.5,
  cat.fontface = "bold",
  cat.col = c("#C34062", "#F5BD4D", "#005493"),
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)

UpSetR::upset(
  input,
  sets = c("DCM_only", "ICM_only", "CTL_only", "HF_specific", "HF_core", "ICM_DCM_diff."),
  order.by = "freq",
  keep.order = TRUE,
  main.bar.color = "grey30",              
  matrix.color = "black", 
  sets.bar.color = "black",
  point.size = 15,
  line.size = 1.5,
  mb.ratio = c(0.7, 0.3),
  text.scale = c(7, 7, 5, 6, 4, 6),
)

load("HF_Plasticizer_Intersection_7S_Data.RData")


plz_input <- fromList(cluster_union_lists)
unique_counts <- colSums(plz_input == 1 & rowSums(plz_input) == 1)

print(unique_counts)
cluster_colors <- c(
  "Cluster A" = "#984EA3",
  "Cluster B" = "#377EB8",
  "Cluster C" = "#4DAF4A",
  "Cluster D" = "#778899",
  "Cluster E" = "#B92202",
  "Cluster F" = "#483D8B"
)

names(cluster_union_lists) <- c(
  "Cluster A", "Cluster B", "Cluster C",
  "Cluster D", "Cluster E", "Cluster F"
)

plz_input <- fromList(cluster_union_lists)
plz_input

UpSetR::upset(
  plz_input,
  sets = rev(names(cluster_union_lists)),
  keep.order = TRUE,
  order.by = "freq",
  main.bar.color = "grey30",              
  matrix.color = "black",                 
  sets.bar.color =  rev(cluster_colors),
  point.size = 6,
  line.size = 0.9,
  mb.ratio = c(0.6, 0.4),
  text.scale = c(0, 6, 4, 2.5, 0, 0)
)


