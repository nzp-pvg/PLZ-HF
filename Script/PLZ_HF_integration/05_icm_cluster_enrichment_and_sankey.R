# Task: Run ICM-focused cluster enrichment analyses and generate ICM cluster/pathway Sankey-ready objects and related exports.
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

ICM_gene_sets <- lapply(cluster_union_lists, function(cluster_genes) {
  intersect(HF_gene_sets$ICM_hgnc, cluster_genes)
})
names(ICM_gene_sets) <- names(cluster_union_lists)
sapply(ICM_gene_sets, length)


library(clusterProfiler)
library(org.Hs.eg.db)

go_list <- lapply(ICM_gene_sets, function(genes) {
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

names(go_list) <- names(ICM_gene_sets)

cluster_ids <- names(go_list)

top_n <- 10

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
    limits = c(0.5, 6),
    oob = scales::squish,
    name = "-log10(p.adjust)"
  ) +
  scale_size(
    range = c(3, 10),
    breaks = seq(1, 8, 2),
    name = "Count"
  ) +
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

dim(dotplot_df_clean)

ICM_only_top5_GO_BP_snk <- p_sankey
ICM_only_top5_GO_BP_dot <- dotplot
ICM_only_top5_GO_BP_snk_sq <- levels(sankey_df$Description)
ICM_only_top5_GO_BP_dot_sq <- levels(sankey_df$Description)
save(ICM_only_top5_GO_BP_snk, ICM_only_top5_GO_BP_dot, ICM_only_top5_GO_BP_snk_sq, ICM_only_top5_GO_BP_dot_sq, file ="ICM_only_top5_GO_BP_res.RData")

ICM_only_top5_GO_CC_snk <- p_sankey
ICM_only_top5_GO_CC_dot <- dotplot
ICM_only_top5_GO_CC_snk_sq <- levels(sankey_df$Description)
ICM_only_top5_GO_CC_dot_sq <- levels(sankey_df$Description)
save(ICM_only_top5_GO_CC_snk, ICM_only_top5_GO_CC_dot, ICM_only_top5_GO_CC_snk_sq, ICM_only_top5_GO_CC_dot_sq, file ="ICM_only_top5_GO_CC_res.RData")


ICM_only_top5_GO_MF_snk <- p_sankey
ICM_only_top5_GO_MF_dot <- dotplot
ICM_only_top5_GO_MF_snk_sq <- levels(sankey_df$Description)
ICM_only_top5_GO_MF_dot_sq <- levels(sankey_df$Description)
save(ICM_only_top5_GO_MF_snk, ICM_only_top5_GO_MF_dot, ICM_only_top5_GO_MF_snk_sq, ICM_only_top5_GO_MF_dot_sq, file ="ICM_only_top5_GO_MF_res.RData")


