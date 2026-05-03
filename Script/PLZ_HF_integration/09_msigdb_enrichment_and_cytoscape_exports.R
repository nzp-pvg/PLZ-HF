# Task: Run MSigDB enrichment and prepare Cytoscape-friendly exports for PLZ-HF integration results.
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

library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)


msigdb_cgp <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CGP"
) %>%
  dplyr::select(gs_name, gene_symbol) 

head(msigdb_cgp)

library(msigdbr)
library(clusterProfiler)
library(dplyr)

msigdb_cgp <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CGP"
) %>%
  dplyr::select(gs_name, gene_symbol)

head(msigdb_cgp)

HF_spec_gene_sets <- lapply(cluster_union_lists, function(cluster_genes) {
  intersect(HF_gene_sets$HF_spec_hgnc, cluster_genes)
})
names(HF_spec_gene_sets) <- names(cluster_union_lists)

sapply(HF_spec_gene_sets, length)

## ====================== GSEA / ORA ======================
ego_cgp_list <- lapply(HF_spec_gene_sets, function(genes) {
  enricher(
    gene          = genes,
    TERM2GENE     = msigdb_cgp,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
})

names(ego_cgp_list) <- names(HF_spec_gene_sets)

top_n <- 10
cgp_cluster_df <- do.call(rbind, lapply(names(ego_cgp_list), function(cl) {
  res <- ego_cgp_list[[cl]]@result
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res <- res[order(res$p.adjust), ][1:min(nrow(res), top_n), ]
  data.frame(
    Cluster     = cl,
    ID          = res$ID,
    Description = res$Description,
    p.adjust    = res$p.adjust,
    Count       = res$Count,
    stringsAsFactors = FALSE
  )
}))

head(cgp_cluster_df)

cluster_colors <- c(
  A = "#984EA3", B = "#377EB8", C = "#4DAF4A",
  D = "#778899", E = "#B92202", F = "#483D8B"
)

cgp_cluster_df <- cgp_cluster_df %>%
  mutate(
    log10p = -log10(p.adjust),
    Description = factor(Description, levels = rev(unique(Description))),
    Cluster = factor(Cluster, levels = names(cluster_colors))
  )

ggplot(cgp_cluster_df, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, fill = log10p),
             shape = 21, stroke = 0.8, color = "black") +
  scale_fill_gradientn(
    colors = c("#4F4F4F", "#008080", "#10C2AA"),
    limits = c(3, 10),
    breaks = c(0, 4, 6, 8, 10),   
    oob = scales::squish,
    name = "-log10(p.adjust)"
  ) +
  scale_size(
    range = c(4, 16),
    breaks = seq(2, 18, 2), 
    name = "Gene Count") +
  labs(title = "HF-specific × Cluster – C2:CGP Enrichment",
       x = "Cluster", y = "Perturbation Gene Set") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    plot.title  = element_text(size = 2, face = "bold"),
    panel.grid.major = element_line(color = "grey8", size = 0.2),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 15)
  )
range(cgp_cluster_df$Count)
range(cgp_cluster_df$log10p)

cat(paste(seq_along(unique(as.character(cgp_cluster_df$Description))),
          unique(as.character(cgp_cluster_df$Description))),
    sep = "\n")

library(dplyr)
top_terms_global_unique <- cgp_cluster_df %>%
  arrange(p.adjust) %>%
  distinct(Description, .keep_all = TRUE) %>%
  slice_head(n = 12)

cat("Top global unique terms:\n", 
    paste(top_terms_global_unique$Description, collapse = "\n"))



top_terms_global_unique$Cluster


library(msigdbr)
library(clusterProfiler)
library(dplyr)

msigdb_tf <- msigdbr(
  species = "Homo sapiens",
  collection = "C3",
  subcollection = "TFT:GTRD"   # Transcription Factor Targets
) %>%
  dplyr::select(gs_name, gene_symbol)

ego_tf_list <- lapply(HF_spec_gene_sets, function(genes) {
  enricher(
    gene          = genes,
    TERM2GENE     = msigdb_tf,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
})

names(ego_tf_list) <- names(HF_spec_gene_sets)

top_n <- 10
tf_cluster_df <- do.call(rbind, lapply(names(ego_tf_list), function(cl) {
  res <- ego_tf_list[[cl]]@result
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res <- res[order(res$p.adjust), ][1:min(nrow(res), top_n), ]
  data.frame(
    Cluster     = cl,
    ID          = res$ID,
    Description = res$Description,
    p.adjust    = res$p.adjust,
    Count       = res$Count,
    stringsAsFactors = FALSE
  )
}))

head(tf_cluster_df)
range(tf_cluster_df$Count)
range(-log10(tf_cluster_df$p.adjust))

ggplot(tf_cluster_df, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)),
             shape = 21, stroke = 0.8, color = "black") +
  scale_fill_gradientn(
    colors = c("darkblue", "#A52A2A","#FC8F8F"),
    limits = c(0.2, 1),
    #breaks = c(0, 4, 6, 8, 10),   
    oob = scales::squish,
    name = "-log10(p.adjust)"
  ) +
  scale_size(
    range = c(2, 13),
    breaks = seq(1, 6, 2), 
    name = "Gene Count") +
  labs(title = "HF-specific × Cluster – C3:TFT Enrichment",
       x = "Cluster", y = "Transcription Factor Target Set") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    plot.title  = element_text(size = 2, face = "bold"),
    panel.grid.major = element_line(color = "grey8", size = 0.2),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 15)
  )


cat(paste(seq_along(unique(as.character(cgp_cluster_df$Description))),
          unique(as.character(cgp_cluster_df$Description))),
    sep = "\n")
tf_cluster_df$p.adjust


tf_cluster_df$ID


## =======================
library(dplyr)


edges_tf <- c3_cluster_df %>%
  transmute(source = Cluster,
            target = Description,
            interaction = "Cluster-TF",
            weight = Count,
            score = -log10(p.adjust))

edges_kegg <- top3_kegg_per_cluster %>%
  transmute(source = Cluster,
            target = Description,
            interaction = "Cluster-KEGG",
            weight = Count,
            score = -log10(p.adjust))

edges_reactome <- reactome_cluster_df %>%
  transmute(source = Cluster,
            target = Pathway,
            interaction = "Cluster-Reactome",
            weight = Count,
            score = -log10(p.adjust))

edges_all <- bind_rows(edges_tf, edges_kegg, edges_reactome)

clusters <- data.frame(
  name = c("A","B","C","D","E","F"),
  Type = "Cluster",
  Cluster = c("A","B","C","D","E","F"),
  ClusterColor = c("#984EA3","#377EB8","#4DAF4A","#778899","#B92202","#483D8B"),
  Label = c("A","B","C","D","E","F")
)

other_nodes <- data.frame(
  name = unique(c(edges_all$target)),
  Type = gsub(".*-", "", edges_all$interaction)[match(unique(c(edges_all$target)), edges_all$target)],
  Cluster = NA,
  ClusterColor = NA,
  Label = unique(c(edges_all$target))
)

nodes_all <- bind_rows(clusters, other_nodes)

write.csv(edges_all, "HF_specific_Cytoscape_edges.csv", row.names = FALSE, quote = FALSE)
write.csv(nodes_all, "HF_specific_Cytoscape_nodes.csv", row.names = FALSE, quote = FALSE)

cat("✅ Exported:", nrow(nodes_all), "nodes and", nrow(edges_all), "edges\n")
