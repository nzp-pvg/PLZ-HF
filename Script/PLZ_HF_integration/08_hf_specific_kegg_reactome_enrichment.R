# Task: Run HF-specific enrichment analyses, including KEGG and Reactome summaries.
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
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(openxlsx)

setwd(go_s4_dir)

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

library(eulerr)
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



length(DEG_ICM)
length(DEG_DCM)
length(DEG_CTL)

length(intersect(DEG_ICM, DEG_DCM))   # ICM ∩ DCM
length(intersect(DEG_ICM, DEG_CTL))   # ICM ∩ CTL
length(intersect(DEG_DCM, DEG_CTL))   # DCM ∩ CTL

length(Reduce(intersect, list(DEG_ICM, DEG_DCM, DEG_CTL)))  # ICM ∩ DCM ∩ CTL

length(union(DEG_ICM, DEG_DCM))       # ICM ∪ DCM
length(union(DEG_ICM, DEG_CTL))       # ICM ∪ CTL
length(union(DEG_DCM, DEG_CTL))       # DCM ∪ CTL

length(Reduce(union, list(DEG_ICM, DEG_DCM, DEG_CTL)))  # ICM ∪ DCM ∪ CTL

icm_dcm_intersect <- intersect(DEG_ICM, DEG_DCM)
length(icm_dcm_intersect)
icm_dcm_not_ctl <- setdiff(icm_dcm_intersect, DEG_CTL)
length(icm_dcm_not_ctl)

icm_ctl_intersect <- intersect(DEG_ICM, DEG_CTL)
length(icm_ctl_intersect)
icm_ctl_not_dcm <- setdiff(icm_ctl_intersect, DEG_DCM)
length(icm_ctl_not_dcm)


dcm_ctl_intersect <- intersect(DEG_DCM, DEG_CTL)
length(dcm_ctl_intersect)
dcm_ctl_not_icm <- setdiff(dcm_ctl_intersect, DEG_ICM)
length(dcm_ctl_not_icm)

HF_spec <- setdiff(union(DEG_ICM, DEG_DCM), DEG_CTL)
HF_core <- intersect(DEG_ICM, DEG_DCM)
ICM_only <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL)) # ICM-only
DCM_only <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL)) # DCM-only
CTL_only <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM)) # CTL-only
ICM_DCM_diff <- setdiff(union(DEG_ICM, DEG_DCM), intersect(DEG_ICM, DEG_DCM))

sapply(list(
  HF_spec = HF_spec,
  HF_core = HF_core,
  ICM_only = ICM_only,
  DCM_only = DCM_only,
  CTL_only = CTL_only,
  ICM_DCM_diff = ICM_DCM_diff
), length)





HF_spec <- setdiff(union(DEG_ICM, DEG_DCM), DEG_CTL)
HF_core <- intersect(DEG_ICM, DEG_DCM)
ICM <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))
ICM_DCM_diff <- setdiff(union(DEG_ICM, DEG_DCM), intersect(DEG_ICM, DEG_DCM))

sapply(list(
  HF_spec = HF_spec,
  HF_core = HF_core,
  ICM_only = ICM,
  DCM_only = DCM,
  CTL_only = CTL,
  ICM_DCM_diff = ICM_DCM_diff
), length)






hf_spec_genes <- HF_gene_sets$HF_spec_hgnc

clusters <- cluster_union_lists

hf_spec_intersections <- lapply(clusters, function(cl) {
  intersect(hf_spec_genes, cl)
})

hf_spec_intersections_df <- data.frame(
  Cluster = names(hf_spec_intersections),
  Genes   = sapply(hf_spec_intersections, function(x) paste(x, collapse = ";")),
  Count   = sapply(hf_spec_intersections, length)
)

hf_spec_genes <- unique(unlist(hf_cluster_intersections$HF_spec_hgnc))




rm(list=ls())
gc()

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)



## ====================== enrichGO (BP) ======================
sapply(cluster_union_lists, function(cl) {
  length(intersect(HF_gene_sets$HF_spec_hgnc, cl))
})


## ====================== enrichGO (HF_spec) ======================
HF_spec_gene_sets <- lapply(cluster_union_lists, function(cluster_genes) {
  intersect(HF_gene_sets$HF_spec_hgnc, cluster_genes)
})
names(HF_spec_gene_sets) <- names(cluster_union_lists)
sapply(HF_spec_gene_sets, length)

ego_list <- lapply(cluster_union_lists, function(cluster_genes) {
  enrichGO(
    gene          = intersect(HF_gene_sets$HF_spec_hgnc, cluster_genes),
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
names(ego_list) <- names(cluster_union_lists)

top_n <- 10
go_cluster_df <- do.call(rbind, lapply(names(ego_list), function(cl) {
  go <- ego_list[[cl]]@result
  if (nrow(go) == 0) return(NULL)
  go <- go[order(go$p.adjust), ][1:min(nrow(go), top_n), ]
  data.frame(
    Cluster     = cl,
    ID          = go$ID,
    Description = paste0(go$ID, "  ", go$Description),
    p.adjust    = go$p.adjust,
    Count       = go$Count,
    stringsAsFactors = FALSE
  )
}))

term_freq <- go_cluster_df %>%
  count(Description, name = "Freq") %>%
  arrange(desc(Freq))

selected_terms <- term_freq %>%
  filter(Freq >= 2) %>%
  pull(Description)

dotplot_df <- go_cluster_df %>%
  filter(Description %in% selected_terms) %>%
  mutate(
    log10p = -log10(p.adjust),
    Description = factor(Description, levels = rev(selected_terms)),
    Cluster = factor(Cluster, levels = c("A","B","C","D","E","F"))
  )


range(dotplot_df$Count)
range(dotplot_df$log10p)
## Step 3: dotplot
ggplot(dotplot_df, aes(x = Cluster, y = Description)) +
  geom_point(
    aes(fill = log10p, size = Count),
    shape = 21, stroke = 0.8, color = "black"
  ) +
  scale_fill_gradientn(
    colors = c("#005493", "#F5BD4D", "#C34062"),
    limits = c(2, 4),
    oob = scales::squish,
    name = "-log10(p.adjust)"
  ) +
  scale_size(
    range = c(4, 20),
    breaks = seq(2, 10, 2),
    name = "Gene Count") +
  scale_x_discrete(labels = names(cluster_colors)) +
  labs(title = "HF_spec – GO Dotplot (BP)", x = "Cluster", y = "GO Term") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey8", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = ggtext::element_markdown(size = 16, face = "bold", hjust = 0.35),
    axis.text.y = element_text(size = 13, face = "bold")
  )



cat("Dotplot display order:\n")
cat(levels(dotplot_df$Description), sep = "\n")
cat(rev(levels(dotplot_df$Description)), sep = "\n")



## KEGG
## ====================== enrichKEGG (HF-specific × Cluster) ======================
library(clusterProfiler)
library(org.Hs.eg.db)

HF_spec_gene_sets_entrez <- lapply(HF_spec_gene_sets, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
})
names(HF_spec_gene_sets_entrez) <- names(HF_spec_gene_sets)

ego_kegg_list <- lapply(HF_spec_gene_sets_entrez, function(genes) {
  enrichKEGG(
    gene          = genes,
    organism      = "hsa",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    minGSSize     = 5,
    maxGSSize     = 500
  )
})
names(ego_kegg_list) <- names(HF_spec_gene_sets_entrez)
head(ego_kegg_list$A@result, 10)

dotplot(ego_kegg_list$A, showCategory = 10)
barplot(ego_kegg_list$A, showCategory = 10)

kegg_df <- ego_kegg_list$A@result
View(kegg_df)

library(ggplot2)
library(dplyr)

kegg_df <- do.call(rbind, lapply(names(ego_kegg_list), function(cl) {
  df <- ego_kegg_list[[cl]]@result
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$Cluster <- cl
  df
}))

top_n <- 10
kegg_df_top <- kegg_df %>%
  group_by(Cluster) %>%
  slice_min(order_by = p.adjust, n = top_n) %>%
  ungroup()

ggplot(kegg_df_top, aes(x = Cluster, y = reorder(Description, -log10(p.adjust)))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low="blue", high="darkred") +
  labs(title = "KEGG Enrichment (Top 10 per Cluster)", 
       x = "Cluster", y = "Pathway") +
  theme_minimal(base_size = 14)



barplot(ego_kegg_list[["A"]], showCategory = 10, title = "Cluster A – KEGG")
barplot(ego_reactome_list[["A"]], showCategory = 10, title = "Cluster A – Reactome")


library(enrichplot)
ego_kegg_list_symbol <- lapply(ego_kegg_list, function(x) {
  if (!is.null(x) && nrow(as.data.frame(x)) > 0) {
    setReadable(x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  } else {
    NULL
  }
})

cnetplot(ego_kegg_list_symbol[["A"]], categorySize="pvalue", foldChange=NULL,
         showCategory = 5, circular = FALSE, colorEdge = TRUE)



library(dplyr)
library(ggplot2)

top3_kegg_per_cluster <- do.call(rbind, lapply(names(ego_kegg_list), function(cl) {
  go <- ego_kegg_list[[cl]]@result
  if (is.null(go) || nrow(go) == 0) return(NULL)
  go <- go[order(go$p.adjust), ][1:min(3, nrow(go)), ]
  go$Cluster <- cl
  go
}))

top3_kegg_per_cluster <- top3_kegg_per_cluster %>%
  mutate(
    Cluster = factor(Cluster, levels = c("A","B","C","D","E","F")),
    Description = factor(Description, levels = unique(Description[order(p.adjust)]))
  )

p_kegg_merge <- ggplot(top3_kegg_per_cluster, aes(x = -log10(p.adjust), 
                                  y = Description, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = seq(1.5, length(levels(top3_kegg_per_cluster$Description)) - 0.5, 1),
             linetype = "solid", color = "black", linewidth = 0.4) +
  labs(title = "Top 3 KEGG Pathways per Cluster", 
       x = "-log10(p.adjust)", y = "Pathway") +
  scale_fill_manual(values = c(
    A="#984EA3", B="#377EB8", C="#4DAF4A",
    D="#778899", E="#B92202", F="#483D8B"
  )) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", lineheight = 2),
    axis.text.x = element_text(size = 40),
    plot.title  = element_text(size = 24, face = "bold")
  )


top3_kegg_per_cluster$Description

## Reactome
library(ReactomePA)

ego_reactome_list <- lapply(HF_spec_gene_sets_entrez, function(genes) {
  enrichPathway(
    gene          = genes,
    organism      = "human",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE,
    minGSSize = 5, maxGSSize = 500
  )
})
names(ego_reactome_list) <- names(HF_spec_gene_sets)


library(dplyr)
library(ggplot2)

top_n <- 3
reactome_cluster_df <- do.call(rbind, lapply(names(ego_reactome_list), function(cl) {
  go <- ego_reactome_list[[cl]]@result
  if (is.null(go) || nrow(go) == 0) return(NULL)
  go <- go[order(go$p.adjust), ][1:min(nrow(go), top_n), ]
  data.frame(
    Cluster  = cl,
    Pathway  = go$Description,
    p.adjust = go$p.adjust,
    Count    = go$Count,
    stringsAsFactors = FALSE
  )
}))

reactome_cluster_df <- reactome_cluster_df %>%
  mutate(
    log10p = -log10(p.adjust),
    Cluster = factor(Cluster, levels = c("A","B","C","D","E","F"))
  )

cluster_colors <- c(
  A = "#984EA3", B = "#377EB8", C = "#4DAF4A",
  D = "#778899", E = "#B92202", F = "#483D8B"
)

p_reactome_merge <- ggplot(reactome_cluster_df, aes(x = log10p, y = Pathway, fill = Cluster)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = seq(1.5, length(levels(top3_kegg_per_cluster$Description)) - 0.5, 1),
             linetype = "solid", color = "black", linewidth = 0.4) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Top 3 Reactome Pathways per Cluster",
       x = "-log10(p.adjust)", y = "Pathway") +
  theme_minimal(base_size = 14) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", lineheight = 1.3),
    axis.text.x = element_text(size = 40),
    plot.title  = element_text(size = 24, face = "bold")
  )




library(patchwork)

p_kegg_merge <- p_kegg_merge +
  labs(title = "Top 3 KEGG Pathways per Cluster") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

p_reactome_merge <- p_reactome_merge +
  labs(title = "Top 3 Reactome Pathways per Cluster") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

combined_plot <- p_kegg_merge / p_reactome_merge +
  plot_annotation(
    title = "Pathway Enrichment (KEGG + Reactome)",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

combined_plot

ggsave("HF_spec_Pathway_KEGG_Reactome.pdf", combined_plot, width = 9, height = 12)
ggsave("HF_spec_Pathway_KEGG_Reactome.png", combined_plot, width = 9, height = 12, dpi = 300)

p_kegg_merge
p_reactome_merge

cat("Reactome pathways:\n", paste(reactome_cluster_df$Pathway, collapse = "\n"))
cat("Top 3 KEGG pathways per cluster:\n", 
    paste(unique(as.character(top3_kegg_per_cluster$Description)), collapse = "\n"))





