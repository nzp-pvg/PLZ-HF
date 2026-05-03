# Task: Construct hub-gene PPI, pathway, and Cytoscape export objects for DCM and ICM cluster-centered analyses.
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

setwd(go_s4_dir)



library(STRINGdb)
library(igraph)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)

# load gene sets
load("HF_DEGs_hgnc_top560.RData") ## HF_gene_set
load("HF_Plasticizer_Intersection_7S_Data.RData") ## cluster_union_lists

HF_gene_sets <- list(
  HF_spec_hgnc = HF_spec_hgnc, 
  ICM_hgnc = ICM_hgnc, 
  DCM_hgnc = DCM_hgnc, 
  ICM_DCM_diff_hgnc = ICM_DCM_diff_hgnc, 
  HF_core_hgnc = HF_core_hgnc, 
  CTL_hgnc = CTL_hgnc)

sapply(HF_gene_sets, length)

DCM_only <- HF_gene_sets$DCM_hgnc

DCM_gene_sets <- lapply(cluster_union_lists, function(cluster_genes) {
  intersect(DCM_only, cluster_genes)
})
names(DCM_gene_sets) <- names(cluster_union_lists)

## HF_gene_sets、DCM_only、cluster_union_lists

top_n <- 60

go_list <- lapply(DCM_gene_sets, function(genes) {
  enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 5, 
    maxGSSize     = 500, 
    readable      = TRUE
  )
})

top_terms <- lapply(go_list, function(er) {
  if (is.null(er) || nrow(er@result) == 0) return(NULL)
  er@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n)
})

extract_genes_from_terms <- function(df) {
  if (is.null(df)) return(NULL)
  unique(unlist(strsplit(df$geneID, "/")))
}

cluster_genes_for_ppi <- lapply(top_terms, extract_genes_from_terms)

ppi_input <- do.call(rbind, lapply(names(cluster_genes_for_ppi), function(cl) {
  genes <- cluster_genes_for_ppi[[cl]]
  if (is.null(genes)) return(NULL)
  data.frame(Cluster = cl, Gene = genes, stringsAsFactors = FALSE)
}))

ppi_merged <- ppi_input %>%
  group_by(Gene) %>%
  summarise(Clusters = paste(unique(Cluster), collapse = ";"),
            .groups = "drop")

string_db <- STRINGdb$new(
  version = "12", 
  species = 9606,
  score_threshold = 400,
  input_directory = ""
)

mapped_global <- string_db$map(data.frame(Gene = DCM_only, stringsAsFactors = FALSE), 
                               "Gene", removeUnmappedRows = TRUE)
ppi_edges_global <- string_db$get_interactions(mapped_global$STRING_id)
ppi_nodes_global <- mapped_global[, c("STRING_id","Gene")]

write.csv(ppi_edges_global, "DCM_global_STRING_edges.csv", row.names = FALSE, quote = FALSE)
write.csv(ppi_nodes_global, "DCM_global_STRING_nodes.csv", row.names = FALSE)

cat("Global PPI:", nrow(ppi_nodes_global), "nodes,", nrow(ppi_edges_global), "edges\n")

ppi_merged_df <- as.data.frame(ppi_merged, stringsAsFactors = FALSE)

mapped_hub <- string_db$map(ppi_merged_df, "Gene", removeUnmappedRows = TRUE)
ppi_edges_hub <- string_db$get_interactions(mapped_hub$STRING_id)
ppi_nodes_hub <- mapped_hub[, c("STRING_id","Gene","Clusters")]

library(dplyr)
ppi_edges_hub_unique <- ppi_edges_hub %>%
  dplyr::mutate(pair = ifelse(from < to,
                              paste(from, to, sep = "_"),
                              paste(to, from, sep = "_"))) %>%
  dplyr::group_by(pair) %>%
  dplyr::summarise(
    from = dplyr::first(from),
    to = dplyr::first(to),
    combined_score = max(combined_score),
    .groups = "drop"
  ) %>%
  dplyr::select(from, to, combined_score)

cat("After deduplication: ", nrow(ppi_edges_hub_unique), " edges\n", sep = "")

write.csv(ppi_edges_hub_unique, "DCM_hub_STRING_top60_edges.csv", row.names = FALSE, quote = FALSE)
write.csv(ppi_nodes_hub, "DCM_hub_STRING_top60_nodes.csv", row.names = FALSE)

cat("Hub PPI:", nrow(ppi_nodes_hub), "nodes,", nrow(ppi_edges_hub), "edges\n")


clusters <- c("A","B","C","D","E","F")

ppi_nodes_expanded <- mapped_hub %>%
  mutate(Clusters_split = strsplit(Clusters, ";")) %>%
  unnest(Clusters_split) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols = c(STRING_id, Gene),
    names_from = Clusters_split,
    values_from = value,
    values_fill = list(value = 0)
  )

missing_clusters <- setdiff(clusters, colnames(ppi_nodes_expanded))
if(length(missing_clusters) > 0){
  for(cl in missing_clusters){
    ppi_nodes_expanded[[cl]] <- 0
  }
}

ppi_nodes_for_cyto <- ppi_nodes_expanded %>%
  mutate(across(all_of(clusters), ~ as.integer(replace_na(., 0)))) %>%
  dplyr::rename(`shared name` = STRING_id) %>%
  rowwise() %>%
  mutate(Cluster = paste(
    clusters[as.logical(c_across(all_of(clusters)))],
    collapse = ";"
  )) %>%
  ungroup()

cluster_colors <- c(
  A = "#984EA3", B = "#377EB8", C = "#4DAF4A",
  D = "#778899", E = "#B92202", F = "#483D8B"
)

ppi_nodes_for_cyto <- ppi_nodes_for_cyto %>%
  rowwise() %>%
  mutate(Color = ifelse(Cluster == "",
                        "#CCCCCC",
                        paste(cluster_colors[strsplit(Cluster, ";")[[1]]],
                              collapse = ";"))) %>%
  ungroup()

write.csv(ppi_nodes_for_cyto, "DCM_hub_STRING_nodes_top60_expanded.csv", row.names = FALSE)

save(HF_gene_sets, cluster_union_lists, 
     DCM_only, DCM_gene_sets, 
     mapped_global, mapped_hub,
     ppi_edges_global, ppi_nodes_global,
     ppi_edges_hub, ppi_nodes_hub, ppi_nodes_for_cyto,
     file ="DCM_PPI_all_res.RData")


## DCM  cluster - gene - KEGG
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

setwd("~/Science/CVD/CVD_MS_1/data/GO/cytoscape/hub_core")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggalluvial)
library(dplyr)
# ==============================
# ==============================
input_file <- "DCM_hub_top60_edges_nodes.csv"
output_prefix <- "DCM"

# ==============================
# ==============================
hub_tbl <- read.csv(input_file, stringsAsFactors = FALSE)

hub_tbl <- hub_tbl %>%
  dplyr::arrange(desc(Degree)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  dplyr::select(Gene, Cluster)

cat("Loaded and retained the top 15 hub genes:\n")
print(hub_tbl$Gene)
hub_tbl$Cluster
# ==============================
# ==============================
hub_tbl <- hub_tbl %>%
  dplyr::mutate(Cluster = strsplit(Cluster, "\\|")) %>%
  tidyr::unnest(Cluster)

# ==============================
# ==============================
gene_entrez <- bitr(unique(hub_tbl$Gene),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

kegg_res <- clusterProfiler::enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.03,
  qvalueCutoff = 0.2
)

if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
  stop("No KEGG enrichment results were returned; check the input gene set.")
}

kegg_df <- as.data.frame(kegg_res)
cat("KEGG enrichment result rows: ", nrow(kegg_df), "\n", sep = "")

# ==============================
# ==============================
kegg_gene_map <- kegg_df %>%
  dplyr::select(Description, geneID) %>%
  dplyr::mutate(Gene = strsplit(geneID, "/")) %>%
  tidyr::unnest(Gene) %>%
  dplyr::select(Gene, Pathway = Description)

kegg_gene_map <- kegg_gene_map %>%
  dplyr::left_join(gene_entrez, by = c("Gene" = "ENTREZID")) %>%
  dplyr::mutate(Gene = SYMBOL) %>%
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::select(Gene, Pathway)

# ==============================
# ==============================
sankey_data <- hub_tbl %>%
  dplyr::left_join(kegg_gene_map, by = "Gene") %>%
  dplyr::filter(!is.na(Pathway)) %>%
  dplyr::distinct(Cluster, Gene, Pathway)

cat("Sankey input rows: ", nrow(sankey_data), "\n", sep = "")

# ==============================
# ==============================
library(dplyr)
library(ggalluvial)
library(ggplot2)
library(tidyr)
sankey_data_equal <- sankey_data %>%
  group_by(Cluster) %>%
  mutate(weight = 1 / n()) %>%
  ungroup()

ggplot(sankey_data_equal,
       aes(axis1 = Cluster, axis2 = Gene, axis3 = Pathway, y = weight)) +
  geom_alluvium(aes(fill = Cluster), width = 1/12, alpha = 0.5, knot.pos = 0.2) +
  geom_stratum(width = 1/10, fill = "grey95", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 0, color = "black") +
  scale_x_discrete(
    limits = c("Plasticizer Cluster", "Hub Genes", "KEGG Pathways"),
    expand = c(0.05, 0.05)
  ) +
  scale_fill_manual(values = c(
    A = "#984EA3", B = "#377EB8", C = "#4DAF4A",
    D = "#778899", E = "#B92202", F = "#483D8B"
  )) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )



colnames(sankey_data)

## ICM cluster - gene - KEGG
# ==============================
# ==============================

input_file <- "ICM_hub_top60_edges_nodes.csv"
output_prefix <- "ICM"

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(clusterProfiler)
library(org.Hs.eg.db)

# ==============================
# ==============================
hub_tbl <- read.csv(input_file, stringsAsFactors = FALSE)

hub_tbl <- hub_tbl %>%
  dplyr::arrange(desc(Degree)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  dplyr::select(Gene, Cluster)

cat("Loaded and retained the top 15 hub genes:\n")
print(hub_tbl$Gene)
hub_tbl$Cluster

# ==============================
# ==============================
hub_tbl <- hub_tbl %>%
  dplyr::mutate(Cluster = strsplit(Cluster, "\\|")) %>%
  tidyr::unnest(Cluster)

# ==============================
# ==============================
gene_entrez <- bitr(unique(hub_tbl$Gene),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

kegg_res <- clusterProfiler::enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
  stop("No KEGG enrichment results were returned; check the input gene set.")
}

kegg_df <- as.data.frame(kegg_res)
# ==============================
# ==============================
kegg_df <- kegg_df %>%
  dplyr::filter(p.adjust < 0.02) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 11)

cat("Retained top 10 KEGG pathways: ", nrow(kegg_df), "\n", sep = "")
cat("KEGG enrichment result rows: ", nrow(kegg_df), "\n", sep = "")

# ==============================
# ==============================
kegg_gene_map <- kegg_df %>%
  dplyr::select(Description, geneID) %>%
  dplyr::mutate(Gene = strsplit(geneID, "/")) %>%
  tidyr::unnest(Gene) %>%
  dplyr::select(Gene, Pathway = Description)

kegg_gene_map <- kegg_gene_map %>%
  dplyr::left_join(gene_entrez, by = c("Gene" = "ENTREZID")) %>%
  dplyr::mutate(Gene = SYMBOL) %>%
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::select(Gene, Pathway)

# ==============================
# ==============================
sankey_data <- hub_tbl %>%
  dplyr::left_join(kegg_gene_map, by = "Gene") %>%
  dplyr::filter(!is.na(Pathway)) %>%
  dplyr::distinct(Cluster, Gene, Pathway)

top_hubs <- kegg_df %>%
  dplyr::select(Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::count(geneID, sort = TRUE) %>%
  dplyr::slice_head(n = 14) %>%
  dplyr::pull(geneID)

sankey_data <- sankey_data %>%
  dplyr::filter(Gene %in% gene_entrez$SYMBOL[gene_entrez$ENTREZID %in% top_hubs])


cat("Sankey input rows: ", nrow(sankey_data), "\n", sep = "")

# ==============================
# ==============================
sankey_data_equal <- sankey_data %>%
  group_by(Cluster) %>%
  mutate(weight = 1 / n()) %>%
  ungroup()

ggplot(sankey_data_equal,
            aes(axis1 = Cluster, axis2 = Gene, axis3 = Pathway, y = weight)) +
  geom_alluvium(aes(fill = Cluster), width = 1/12, alpha = 0.6, knot.pos = 0.2) +
  geom_stratum(width = 1/10, fill = "grey95", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 0, color = "black") +
  scale_x_discrete(
    limits = c("Plasticizer Cluster", "Hub Genes", "KEGG Pathways"),
    expand = c(0.05, 0.05)
  ) +
  scale_fill_manual(values = c(
    A = "#984EA3", B = "#377EB8", C = "#4DAF4A",
    D = "#778899", E = "#B92202", F = "#483D8B"
  )) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )


colnames(sankey_data)

sankey_data$Gene

sankey_data$Pathway


cat(head(unique(sankey_data$Pathway), 20), sep = "\n")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

cluster_cols <- c(
  "A"="#984EA3","B"="#377EB8","C"="#4DAF4A",
  "D"="#778899","E"="#B92202","F"="#483D8B"
)

links_C2A <- tibble::tribble(
  ~Cluster, ~Axis,
  "A", "Stress/DNA damage (p53–Skp2–p27/p21)",
  "B", "Cytokine / NF-κB / IL-17",
  "C", "AGE–RAGE (diabetic complications)",
  "C", "Chemokine receptors (CCR/CXCR)",
  "C", "Viral / infection-like (KSHV/HPV)",
  "D", "MMP activation",
  "D", "Matrix remodeling",
  "E", "MAPK signaling cascade",
  "E", "VEGF receptor activation",
  "F", "Endocrine resistance",
  "F", "Hormonal regulation"
)

# A → DCM leaning；B → ICM dominant；
# C：AGE–RAGE → DCM>ICM；Chemokine/Viral → ICM>DCM；
# D（MMP/Matrix）→ Shared；E（MAPK/VEGF）→ DCM dominant；
links_A2D <- tibble::tribble(
  ~Axis, ~Disease, ~w,
  # A：DCM-leaning
  "Stress/DNA damage (p53–Skp2–p27/p21)", "DCM", 2,
  "Stress/DNA damage (p53–Skp2–p27/p21)", "ICM", 1,
  
  # B：ICM-dominant
  "Cytokine / NF-κB / IL-17", "ICM", 2,
  "Cytokine / NF-κB / IL-17", "DCM", 1,
  
  "AGE–RAGE (diabetic complications)", "DCM", 3,
  "AGE–RAGE (diabetic complications)", "ICM", 1,
  
  "Chemokine receptors (CCR/CXCR)", "ICM", 4,
  "Chemokine receptors (CCR/CXCR)", "DCM", 1,
  
  "Viral / infection-like (KSHV/HPV)", "ICM", 5,
  "Viral / infection-like (KSHV/HPV)", "DCM", 0,
  
  "MMP activation",    "ICM", 1,
  "MMP activation",    "DCM", 1,
  "Matrix remodeling", "ICM", 1,
  "Matrix remodeling", "DCM", 1,
  
  # E：DCM-dominant
  "MAPK signaling cascade", "DCM", 2,
  "MAPK signaling cascade", "ICM", 1,
  "VEGF receptor activation","DCM", 2,
  "VEGF receptor activation","ICM", 1,
  
  "Endocrine resistance", "DCM", 2,
  "Endocrine resistance", "ICM", 1,
  "Hormonal regulation",  "DCM", 2,
  "Hormonal regulation",  "ICM", 1
)

cluster_bias <- tibble::tribble(
  ~Cluster, ~Disease, ~bias,
  "A","DCM",1.6,  "A","ICM",0.6,
  "B","ICM",1.6,  "B","DCM",0.6,
  "C","DCM",1.0,  "C","ICM",1.0,
  "D","DCM",1.0,  "D","ICM",1.0,
  "E","DCM",1.5,  "E","ICM",0.7,
  "F","DCM",1.3,  "F","ICM",0.9
)

df <- links_C2A %>%
  left_join(links_A2D,   by = "Axis") %>%
  left_join(cluster_bias, by = c("Cluster","Disease")) %>%
  mutate(bias = ifelse(is.na(bias), 1, bias),
         w_adj = w * bias) %>%
  group_by(Cluster) %>%
  mutate(y = w_adj / sum(w_adj)) %>%
  ungroup()

axis_levels <- c(
  "Stress/DNA damage (p53–Skp2–p27/p21)",
  axis_B_title,
  "AGE–RAGE (diabetic complications)",
  "Chemokine receptors (CCR/CXCR)",
  "Viral / infection-like (KSHV/HPV)",
  "MMP activation",
  "Matrix remodeling",
  "MAPK signaling cascade",
  "VEGF receptor activation",
  "Endocrine resistance",
  "Hormonal regulation"
)
df$Axis    <- factor(df$Axis,    levels = axis_levels)
df$Cluster <- factor(df$Cluster, levels = c("A","B","C","D","E","F"))
df$Disease <- factor(df$Disease, levels = c("DCM","ICM"))
#df$Disease <- factor(df$Disease, levels = c("ICM", "DCM"))
ggplot(
  df,
  aes(axis1 = Cluster, axis2 = Axis, axis3 = Disease, y = y)
) +
  geom_alluvium(aes(fill = Cluster), width = 1/12, alpha = 0.7, knot.pos = 0.30) +
  geom_stratum(width = 1/10, fill = "grey95", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Plasticizer Clusters", "Mechanistic Axes", "HF subtype"),
    expand = c(0.04, 0.04)
  ) +
  scale_fill_manual(values = cluster_cols) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold")
  )









ICM_only <- HF_gene_sets$ICM_hgnc

ICM_gene_sets <- lapply(cluster_union_lists, function(cluster_genes) {
  intersect(ICM_only, cluster_genes)
})
names(ICM_gene_sets) <- names(cluster_union_lists)

## HF_gene_sets、ICM_only、cluster_union_lists

top_n <- 60

go_list <- lapply(ICM_gene_sets, function(genes) {
  enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 5, 
    maxGSSize     = 500, 
    readable      = TRUE
  )
})

top_terms <- lapply(go_list, function(er) {
  if (is.null(er) || nrow(er@result) == 0) return(NULL)
  er@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n)
})

extract_genes_from_terms <- function(df) {
  if (is.null(df)) return(NULL)
  unique(unlist(strsplit(df$geneID, "/")))
}

cluster_genes_for_ppi <- lapply(top_terms, extract_genes_from_terms)

ppi_input <- do.call(rbind, lapply(names(cluster_genes_for_ppi), function(cl) {
  genes <- cluster_genes_for_ppi[[cl]]
  if (is.null(genes)) return(NULL)
  data.frame(Cluster = cl, Gene = genes, stringsAsFactors = FALSE)
}))

ppi_merged <- ppi_input %>%
  group_by(Gene) %>%
  summarise(Clusters = paste(unique(Cluster), collapse = ";"),
            .groups = "drop")

string_db <- STRINGdb$new(
  version = "12", 
  species = 9606,
  score_threshold = 700,
  input_directory = ""
)

mapped_global <- string_db$map(
  data.frame(Gene = ICM_only, stringsAsFactors = FALSE), 
  "Gene", removeUnmappedRows = TRUE
)

ppi_edges_global <- string_db$get_interactions(mapped_global$STRING_id)
ppi_nodes_global <- mapped_global[, c("STRING_id","Gene")]

library(dplyr)

degree_table_global <- ppi_edges_global %>%
  dplyr::select(from, to) %>%
  tidyr::pivot_longer(cols = everything(), values_to = "STRING_id") %>%
  count(STRING_id, name = "Degree")

ppi_nodes_global <- ppi_nodes_global %>%
  left_join(degree_table_global, by = "STRING_id") %>%
  mutate(Degree = replace_na(Degree, 0))

write.csv(ppi_edges_global, "ICM_global_STRING_edges.csv", row.names = FALSE, quote = FALSE)
write.csv(ppi_nodes_global, "ICM_global_STRING_nodes.csv", row.names = FALSE)

cat("Global PPI:", nrow(ppi_nodes_global), "nodes,", 
    nrow(ppi_edges_global), "edges\n")

ppi_merged_df <- as.data.frame(ppi_merged, stringsAsFactors = FALSE)

mapped_hub <- string_db$map(ppi_merged_df, "Gene", removeUnmappedRows = TRUE)
ppi_edges_hub <- string_db$get_interactions(mapped_hub$STRING_id)
ppi_nodes_hub <- mapped_hub[, c("STRING_id","Gene","Clusters")]

library(dplyr)
ppi_edges_hub_unique <- ppi_edges_hub %>%
  dplyr::mutate(pair = ifelse(from < to,
                              paste(from, to, sep = "_"),
                              paste(to, from, sep = "_"))) %>%
  dplyr::group_by(pair) %>%
  dplyr::summarise(
    from = dplyr::first(from),
    to = dplyr::first(to),
    combined_score = max(combined_score),
    .groups = "drop"
  ) %>%
  dplyr::select(from, to, combined_score)

cat("After deduplication: ", nrow(ppi_edges_hub_unique), " edges\n", sep = "")

write.csv(ppi_edges_hub_unique, "ICM_hub_STRING_top60_edges.csv", row.names = FALSE, quote = FALSE)
write.csv(ppi_nodes_hub, "ICM_hub_STRING_top60_nodes.csv", row.names = FALSE)

cat("Hub PPI:", nrow(ppi_nodes_hub), "nodes,", nrow(ppi_edges_hub), "edges\n")


clusters <- c("A","B","C","D","E","F")

ppi_nodes_expanded <- mapped_hub %>%
  mutate(Clusters_split = strsplit(Clusters, ";")) %>%
  unnest(Clusters_split) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols = c(STRING_id, Gene),
    names_from = Clusters_split,
    values_from = value,
    values_fill = list(value = 0)
  )

missing_clusters <- setdiff(clusters, colnames(ppi_nodes_expanded))
if(length(missing_clusters) > 0){
  for(cl in missing_clusters){
    ppi_nodes_expanded[[cl]] <- 0
  }
}

ppi_nodes_for_cyto <- ppi_nodes_expanded %>%
  mutate(across(all_of(clusters), ~ as.integer(replace_na(., 0)))) %>%
  dplyr::rename(`shared name` = STRING_id) %>%
  rowwise() %>%
  mutate(Cluster = paste(
    clusters[as.logical(c_across(all_of(clusters)))],
    collapse = ";"
  )) %>%
  ungroup()

cluster_colors <- c(
  A = "#984EA3", B = "#377EB8", C = "#4DAF4A",
  D = "#778899", E = "#B92202", F = "#483D8B"
)

ppi_nodes_for_cyto <- ppi_nodes_for_cyto %>%
  rowwise() %>%
  mutate(Color = ifelse(Cluster == "",
                        "#CCCCCC",
                        paste(cluster_colors[strsplit(Cluster, ";")[[1]]],
                              collapse = ";"))) %>%
  ungroup()

write.csv(ppi_nodes_for_cyto, "ICM_hub_STRING_nodes_top60_expanded.csv", row.names = FALSE)

save(HF_gene_sets, cluster_union_lists, 
     ICM_only, ICM_gene_sets, 
     mapped_global, mapped_hub,
     ppi_edges_global, ppi_nodes_global,
     ppi_edges_hub, ppi_nodes_hub, ppi_nodes_for_cyto,
     file ="ICM_PPI_all_res.RData")



