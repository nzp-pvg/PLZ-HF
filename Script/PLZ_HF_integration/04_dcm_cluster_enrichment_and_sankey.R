# Task: Run DCM-focused cluster enrichment analyses and generate DCM cluster/pathway Sankey-ready objects and related exports.
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
library(ggplot2)
library(enrichplot)
library(dplyr)
library(ggnewscale)
library(stringr)
library(org.Hs.eg.db)

cluster_union_lists <- lapply(functional_clusters, function(plasticizers) {
  genes <- unlist(plasticizer_gene_lists_raw[plasticizers], use.names = FALSE)
  unique(genes)
})

all_clusters_genes <- unique(unlist(cluster_union_lists))
DCM_all <- intersect(HF_gene_sets$DCM_hgnc, all_clusters_genes)
ICM_all <- intersect(HF_gene_sets$ICM_hgnc, all_clusters_genes)
HF_core_all <- intersect(HF_gene_sets$HF_core_hgnc, all_clusters_genes)
length(HF_gene_sets$DCM_hgnc)
length(HF_gene_sets$ICM_hgnc)
length(all_clusters_genes)
length(DCM_all)
length(ICM_all)
length(HF_core_all)


# enrichGO for DCM_only
go_res_DCM_all <- enrichGO(
  gene = DCM_all,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.3,
  readable = TRUE
)

go_res_DCM_all@result

dotplot(go_res_DCM_all, showCategory = 15, font.size = 12, 
        title = "GO Enrichment for DCM and All Clusters")
library(ggplot2)

n_genes <- length(DCM_all)

dotplot(go_res_DCM_all, showCategory = 15, font.size = 12) +
  ggtitle("GO Enrichment for DCM and All Clusters",
          subtitle = paste0("n = ", n_genes))

DCM_gene_sets <- lapply(cluster_union_lists, function(cluster_genes) {
  intersect(HF_gene_sets$DCM_hgnc, cluster_genes)
})
names(DCM_gene_sets) <- names(cluster_union_lists)
sapply(DCM_gene_sets, length)


library(clusterProfiler)
library(org.Hs.eg.db)

go_list <- lapply(DCM_gene_sets, function(genes) {
  enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize = 5, 
    maxGSSize = 500, 
    readable      = TRUE
  )
})


names(go_list) <- names(DCM_gene_sets)

cluster_ids <- names(go_list)

top_n <- 6

go_cluster_df <- do.call(rbind, lapply(cluster_ids, function(cl) {
  go <- go_list[[cl]]@result
  if (nrow(go) == 0) return(NULL)
  go <- go[order(go$p.adjust), ][1:min(nrow(go), top_n), ]
  data.frame(Cluster = cl, ID = go$ID, Description = go$Description, stringsAsFactors = FALSE)
}))


library(dplyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(ggtext)

shorten_term <- function(term, n_words = 5) {
  sapply(strsplit(term, " "), function(x) paste(head(x, n_words), collapse = " "))
}

p_thresh    <- 0.05
count_thresh <- 3
cluster_colors <- c(
  A = "#984EA3", B = "#377EB8", C = "#4DAF4A",
  D = "#778899", E = "#B92202", F = "#483D8B"
)

all_top_go <- do.call(rbind, lapply(names(go_list), function(cl) {
  df <- go_list[[cl]]@result 
  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$p.adjust), ][1:min(nrow(df), top_n), ]
  df$Cluster    <- cl
  df$short_desc <- shorten_term(df$Description)
  df$label      <- paste0(df$ID, "  ", df$short_desc)
  df
}))

#colnames(go_list[[1]]@result)
all_top_go <- as_tibble(all_top_go)


used_go_terms <- unique(all_top_go$label)
used_go_terms <- used_go_terms[order(sub("^[^ ]+  ", "", used_go_terms))]

sankey_df <- all_top_go %>%
  dplyr::select(Cluster, label) %>%
  dplyr::count(Cluster, label, name = "Freq") %>%
  dplyr::rename(Description = label)

sankey_df$Cluster     <- factor(sankey_df$Cluster, levels = names(cluster_colors))
sankey_df$Description <- factor(sankey_df$Description, levels = used_go_terms)

dotplot_df <- do.call(rbind, lapply(names(go_list), function(cl) {
  df <- go_list[[cl]]@result
  if (nrow(df) == 0) return(NULL)
  df$Cluster <- cl
  df$short_desc <- shorten_term(df$Description)
  df$label <- paste0(df$ID, "  ", df$short_desc)
  df <- df[df$label %in% used_go_terms, ]
  df[, c("Cluster", "label", "GeneRatio", "p.adjust", "Count")]
})) %>%
  as_tibble() %>%
  mutate(
    sig      = ifelse(p.adjust < p_thresh & Count >= count_thresh, TRUE, FALSE),
    log10p   = ifelse(sig, -log10(p.adjust), NA),
    GeneRatio = sapply(GeneRatio, function(x) {
      nums <- strsplit(x, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    })
  ) %>%
  dplyr::rename(Description = all_of("label"))

dotplot_df$Cluster     <- factor(dotplot_df$Cluster, levels = names(cluster_colors))
dotplot_df$Description <- factor(dotplot_df$Description, levels = used_go_terms)

p_sankey_base <- ggplot(sankey_df,
                        aes(axis1 = Cluster, axis2 = Description, y = Freq)) +
  geom_alluvium(aes(fill = Cluster), width = 1/12, alpha = 0.85) +
  geom_stratum(width = 1/5, fill = "white", color = "black", size = 0.4) +
  scale_fill_manual(values = cluster_colors) +
  scale_x_discrete(limits = c("Cluster", "GO term"), expand = c(.1, .05)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = ggplot2::margin(5, 30, 5, 5))

ggdat <- ggplot_build(p_sankey_base)$data[[2]]

go_lab <- ggdat[ggdat$x == 2, ] |>
  dplyr::group_by(stratum) |>
  dplyr::summarise(x = mean(x), y = mean(y), .groups="drop") |>
  dplyr::mutate(xlab = x + 0.001)

clu_lab <- ggdat[ggdat$x == 1, ] %>%
  dplyr::group_by(stratum) %>%
  dplyr::summarise(x = mean(x), y = mean(y), .groups="drop")


p_sankey <- p_sankey_base +
  geom_text(data = clu_lab, aes(x = x, y = y, label = stratum),
            size = 3, inherit.aes = FALSE) +
  geom_text(data = go_lab, aes(x = xlab, y = y, label = stratum),
            size = 3, hjust = 0, lineheight = .95, inherit.aes = FALSE)

## ==== dotplot ====
dotplot_df_clean <- dotplot_df %>%
  group_by(Cluster, Description) %>%
  slice_min(order_by = p.adjust, with_ties = FALSE) %>%
  ungroup()

colored_labels <- sapply(levels(dotplot_df$Cluster), function(cl) {
  paste0("<span style='color:", cluster_colors[cl], "'><b>", cl, "</b></span>")
})

range(dotplot_df_clean$Count)
range(-log10(dotplot_df_clean$p.adjust))


dotplot <- ggplot(dotplot_df_clean, aes(x = Cluster, y = Description)) +
  geom_point(
    aes(fill = -log10(p.adjust), size = Count),
    shape = 21, stroke = 0.8, color = "black"
  ) +
  scale_fill_gradientn(
    colors = c("#005493", "#F5BD4D", "#C34062"),
    limits = c(1, 4.5),
    oob = scales::squish,
    name = "-log10(p.adjust)") +
  scale_size(range = c(1, 14), breaks = seq(1, 7, 2), name = "Count") +
  scale_x_discrete(labels = colored_labels) +
  scale_y_discrete(limits = rev(used_go_terms)) +
  labs(title = "GO Dotplot", x = "Cluster", y = "GO term") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey8", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = ggtext::element_markdown(size = 16, face = "bold", hjust = 0.35),
    axis.text.y = element_text(size = 13, face = "bold")
  )

#p_sankey + dotplot + plot_layout(widths = c(2, 1))
dotplot
p_sankey

cat("Sankey order:\n")
cat(levels(sankey_df$Description), sep = "\n")

cat("Dotplot display order:\n")
cat(levels(dotplot_df_clean$Description), sep = "\n")


DCM_only_top5_GO_CC_snk <- p_sankey
DCM_only_top5_GO_CC_dot <- dotplot
DCM_only_top5_GO_CC_snk_sq <- levels(sankey_df$Description)
DCM_only_top5_GO_CC_dot_sq <- levels(sankey_df$Description)
save(DCM_only_top5_GO_CC_snk, DCM_only_top5_GO_CC_dot, DCM_only_top5_GO_CC_snk_sq, DCM_only_top5_GO_CC_dot_sq, file ="DCM_only_top5_GO_CC_res.RData")

DCM_only_top5_GO_MF_snk <- p_sankey
DCM_only_top5_GO_MF_dot <- dotplot
DCM_only_top5_GO_MF_snk_sq <- levels(sankey_df$Description)
DCM_only_top5_GO_MF_dot_sq <- levels(sankey_df$Description)
save(DCM_only_top5_GO_MF_snk, DCM_only_top5_GO_MF_dot, DCM_only_top5_GO_MF_snk_sq, DCM_only_top5_GO_MF_dot_sq, file ="DCM_only_top5_GO_MF_res.RData")



DCM_only_top5_GO_BP_snk <- p_sankey
DCM_only_top5_GO_BP_dot <- dotplot
DCM_only_top5_GO_BP_snk_sq <- levels(sankey_df$Description)
DCM_only_top5_GO_BP_dot_sq <- levels(sankey_df$Description)
save(DCM_only_top5_GO_BP_snk, DCM_only_top5_GO_BP_dot, DCM_only_top5_GO_BP_snk_sq, DCM_only_top5_GO_BP_dot_sq, file ="DCM_only_top5_GO_BP_res.RData")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggalluvial)
library(dplyr)

# =========================
# =========================
hub_tbl <- read.csv("DCM_hub_top60_edges_nodes.csv")
hub_tbl <- hub_tbl %>% 
  arrange(desc(Degree)) %>%
  slice_head(n = 15) %>%
  select(Gene, Cluster)

# =========================
# =========================
kegg_res <- enrichKEGG(
  gene = hub_tbl$Gene,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05
)
kegg_df <- as.data.frame(kegg_res)

kegg_gene_map <- kegg_df %>%
  select(Description, geneID) %>%
  mutate(Gene = strsplit(geneID, "/")) %>%
  tidyr::unnest(Gene) %>%
  select(Gene, Pathway = Description)

# =========================
# =========================
sankey_data <- hub_tbl %>%
  left_join(kegg_gene_map, by = "Gene") %>%
  distinct(Cluster, Gene, Pathway)

# =========================
# =========================
ggplot(sankey_data,
       aes(axis1 = Cluster, axis2 = Gene, axis3 = Pathway,
           y = 1)) +
  geom_alluvium(aes(fill = Cluster), width = 1/8) +
  geom_stratum(width = 1/8, fill = "grey90", color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.2, color = "black") +
  scale_x_discrete(limits = c("Cluster", "Hub Genes", "KEGG Pathways"),
                   expand = c(0.05, 0.05)) +
  theme_minimal(base_size = 13) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        panel.grid = element_blank())




