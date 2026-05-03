# Task: Build plasticizer cluster gene unions and compute HF-versus-plasticizer cluster intersections.
# Manuscript mapping: PLZ-HF integration module; see 00_CODE_DIRECTORY.md for task scope.

DCM_hgnc <- convert_deg(DCM, id_mapping)
CTL_hgnc <- convert_deg(CTL, id_mapping)
ICM_DCM_diff_hgnc  <- convert_deg(ICM_DCM_diff, id_mapping)
HF_spec_hgnc  <- convert_deg(HF_spec, id_mapping)
HF_core_hgnc  <- convert_deg(HF_core, id_mapping)


HF_gene_sets <- list(
  HF_spec_hgnc = HF_spec_hgnc, 
  ICM_hgnc = ICM_hgnc, 
  DCM_hgnc = DCM_hgnc, 
  ICM_DCM_diff_hgnc = ICM_DCM_diff_hgnc, 
  HF_core_hgnc = HF_core_hgnc, 
  CTL_hgnc = CTL_hgnc)
sapply(HF_gene_sets, length)

# save
save(HF_spec_hgnc, ICM_hgnc, DCM_hgnc, ICM_DCM_diff_hgnc, HF_core_hgnc, CTL_hgnc, file="HF_DEGs_hgnc_top560.RData")

library(readxl)
library(dplyr)

functional_clusters <- list(
  A = c("ATBC", "ATE", "TEC", "TBC", "TOC", "THC"),
  B = c("DCHP", "BBP", "DBP", "DEP", "DIBP"),
  C = c("IPPP", "TCP", "TPP"),
  D = c("TBP", "TOP", "DEHA", "DOS", "DOZ", "DIDA", "DBS", "DOA"),
  E = c("EPS", "DOTP", "NODTM", "TOTM", "DEHP", "DIHP", "DOP", "DINCH", "DIOP", "DIDP", "DINP"),
  F = c("AMG", "DMP", "ESO", "DPOP")
)

# plasticizer_df <- read_excel("Unique_Gene_list_4s.xlsx")
plasticizer_df <- read_excel("Sorted_Plasticizer_Gene_7s.xlsx")

gene_counts <- data.frame(
  Plasticizer = names(plasticizer_df),
  Gene_Count = sapply(plasticizer_df, function(x) sum(!is.na(x)))
)

cluster_map <- unlist(lapply(names(functional_clusters), function(cluster) {
  setNames(rep(cluster, length(functional_clusters[[cluster]])), functional_clusters[[cluster]])
}))
gene_counts$Cluster <- cluster_map[gene_counts$Plasticizer]

cluster_order <- c("A", "B", "C", "D", "E", "F")
gene_counts$ClusterOrder <- match(gene_counts$Cluster, cluster_order)
gene_counts <- gene_counts %>%
  arrange(ClusterOrder, Cluster, match(Plasticizer, unlist(functional_clusters)))

print(gene_counts[, c("Cluster", "Plasticizer", "Gene_Count")])


functional_clusters <- list(
  A = c("ATBC", "ATE", "TEC", "TBC", "TOC", "THC"),
  B = c("DCHP", "BBP", "DBP", "DEP", "DIBP"),
  C = c("IPPP", "TCP", "TPP"),
  D = c("TBP", "TOP", "DEHA", "DOS", "DOZ", "DIDA", "DBS", "DOA"),
  E = c("EPS", "DOTP", "NODTM", "TOTM", "DEHP", "DIHP", "DOP", "DINCH", "DIOP", "DIDP", "DINP"),
  F = c("AMG", "DMP", "ESO", "DPOP")
)

plasticizer_gene_lists_raw <- as.list(plasticizer_df)
names(plasticizer_gene_lists_raw) <- names(plasticizer_df)

plasticizer_gene_lists_raw <- lapply(plasticizer_gene_lists_raw, function(x) x[!is.na(x)])

cluster_union_lists <- lapply(functional_clusters, function(plasticizers) {
  genes <- unlist(plasticizer_gene_lists_raw[plasticizers], use.names = FALSE)
  unique(genes)
})

sapply(cluster_union_lists, length)


cluster_A_genes <- cluster_union_lists[["A"]]
cluster_B_genes <- cluster_union_lists[["B"]]
cluster_C_genes <- cluster_union_lists[["C"]]
cluster_D_genes <- cluster_union_lists[["D"]]
cluster_E_genes <- cluster_union_lists[["E"]]
cluster_F_genes <- cluster_union_lists[["F"]]


## HF vs PLZ
HF_gene_sets <- list(
  HF_spec_hgnc = HF_spec_hgnc, 
  ICM_hgnc = ICM_hgnc, 
  DCM_hgnc = DCM_hgnc, 
  ICM_DCM_diff_hgnc = ICM_DCM_diff_hgnc, 
  HF_core_hgnc = HF_core_hgnc, 
  CTL_hgnc = CTL_hgnc)



hf_cluster_intersections <- lapply(HF_gene_sets, function(deg_genes) {
  lapply(cluster_union_lists, function(cluster_genes) {
    intersect(deg_genes, cluster_genes)
  })
})

intersection_count_matrix <- sapply(HF_gene_sets, function(deg_genes) {
  sapply(cluster_union_lists, function(cluster_genes) {
    length(intersect(deg_genes, cluster_genes))
  })
})

intersection_count_matrix <- t(intersection_count_matrix)
print(intersection_count_matrix)
