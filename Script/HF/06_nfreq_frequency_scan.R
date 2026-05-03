# === HF_DEG_nFreq_Modeling.R ===
rm(list = ls())
gc()

library(dplyr)
library(edgeR)
library(limma)
library(caret)
library(randomForest)
library(doParallel)


# Load data
load("Discovery_start_data.RData")  # expr_mat, classes_df
Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")  # ICM/DCM/CTL list

# Parameters
nFreq_range <- c(
  1:10,
  seq(12, 30, by = 2),
  seq(35, 70, by = 5),
  80, 90, 100
)
subclasses <- c("ICM", "DCM", "CTL")
repeat_n <- 20
k <- 10
set.seed(2025)

cv_ctrl <- trainControl(
  method = "cv", number = k, classProbs = TRUE, savePredictions = "final"
)

# Setup parallel
cl <- makeCluster(12)
registerDoParallel(cl)

# Start loop
All_Results <- list()

for (sub in subclasses) {
  cat("== Subclass:", sub, "==\n")
  sub_results <- list()

  feature_pool <- names(Feature_list[[sub]])  # sorted by freq desc
  freq_table <- Feature_list[[sub]]

  for (nFreq in nFreq_range) {
    cat("-- nFreq =", nFreq, "--\n")

    selected <- names(freq_table[freq_table >= nFreq])
    if (length(selected) < 10) {
      warning(paste("Too few DEGs at nFreq", nFreq, "for", sub))
      next
    }

    aucs <- c()
    for (r in 1:repeat_n) {
      idx <- createFolds(classes_df$Classes, k = k, list = TRUE)
      pred_all <- data.frame()

      for (i in 1:k) {
        test_idx <- idx[[i]]
        train_idx <- setdiff(1:nrow(classes_df), test_idx)

        TrainData <- expr_mat[selected, train_idx]
        TestData <- expr_mat[selected, test_idx]

        TrainPheno <- classes_df[train_idx, ]
        TestPheno <- classes_df[test_idx, ]

        # Recode One-vs-Rest
        TrainPheno$Label <- ifelse(TrainPheno$Classes == sub, "One", "Others")
        TestPheno$Label <- ifelse(TestPheno$Classes == sub, "One", "Others")

        model <- train(
          x = t(TrainData),
          y = factor(TrainPheno$Label),
          method = "rf",
          trControl = cv_ctrl,
          metric = "Kappa",
          tuneGrid = data.frame(mtry = floor(length(selected) / 3))
        )

        prob <- predict(model, newdata = t(TestData), type = "prob")
        actual <- TestPheno$Label

        library(pROC)
        roc_obj <- roc(response = factor(actual), predictor = prob$One, quiet = TRUE)
        aucs <- c(aucs, as.numeric(auc(roc_obj)))
      }
    }
    sub_results[[paste0("nFreq_", nFreq)]] <- list(
      nGenes = length(selected),
      AUCs = aucs,
      meanAUC = mean(aucs),
      sdAUC = sd(aucs)
    )
  }
  All_Results[[sub]] <- sub_results
}

stopCluster(cl)
saveRDS(All_Results, "HF_nFreq_Scan_25xCV_Results.rds")
cat("✅ All subclass modeling complete. Result saved.\n")


results <- readRDS("HF_nFreq_Scan_25xCV_top560_Results.rds")


names(results$DCM)
str(results$DCM[[1]])


extract_auc_curve <- function(subclass) {
  res_list <- results[[subclass]]
  nFreq_values <- as.numeric(gsub("nFreq_", "", names(res_list)))
  mean_auc <- sapply(res_list, function(x) x$meanAUC)
  sd_auc <- sapply(res_list, function(x) x$sdAUC)
  
  df <- data.frame(
    subclass = subclass,
    nFreq = nFreq_values,
    meanAUC = mean_auc,
    sdAUC = sd_auc
  )
  df[order(df$nFreq), ]
}

ICM_auc <- extract_auc_curve("ICM")
DCM_auc <- extract_auc_curve("DCM")
CTL_auc <- extract_auc_curve("CTL")

#results <- readRDS("HF_nFreq_Scan_25xCV_Results.rds")
all_auc <- bind_rows(ICM_auc, DCM_auc, CTL_auc)


library(dplyr)
library(tidyr)
library(purrr)

raw_auc_list <- list()

for (sub in names(results)) {
  sublist <- results[[sub]]
  for (nFreq_name in names(sublist)) {
    auc_vec <- sublist[[nFreq_name]]$AUCs
    df <- data.frame(
      subclass = sub,
      nFreq = as.numeric(gsub("nFreq_", "", nFreq_name)),
      repeat_id = 1:length(auc_vec),
      AUC = auc_vec
    )
    raw_auc_list[[paste(sub, nFreq_name, sep = "_")]] <- df
  }
}

AUC_raw_matrix <- bind_rows(raw_auc_list)

library(dplyr)
library(tidyr)

AUC_raw_matrix_wide <- AUC_raw_matrix %>%
  group_by(subclass, nFreq) %>%
  mutate(AUC_index = paste0("AUC_", row_number())) %>%
  ungroup() %>%
  select(subclass, nFreq, AUC_index, AUC) %>%
  pivot_wider(
    names_from = AUC_index,
    values_from = AUC
  )

glimpse(AUC_raw_matrix_wide)

library(writexl)
write_xlsx(
  list(
    Summary_AUC = all_auc,
    Raw_AUC = AUC_raw_matrix_wide
  ),
  path = "HF_nFreq_25xCV_top560_AUC_RawMatrix.xlsx"
)


library(ggplot2)
ggplot(all_auc, aes(x = nFreq, y = meanAUC, color = subclass)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = meanAUC - sdAUC, ymax = meanAUC + sdAUC), width = 1) +
  theme_minimal() +
  labs(
    title = "AUC vs nFreq across subclasses",
    x = "nFreq threshold (minimum DEG frequency)",
    y = "Mean AUC ± SD"
  )

ICM_auc[which.max(ICM_auc$meanAUC), ]
DCM_auc[which.max(DCM_auc$meanAUC), ]
CTL_auc[which.max(CTL_auc$meanAUC), ]




Feature_list <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")

deg_file <- readRDS("HF_DEG_Frequency_top560_dynamic_degTable.rds")

DEG_ICM_n22 <- names(deg_file$ICM[deg_file$ICM >= 22])
DEG_DCM_n7  <- names(deg_file$DCM[deg_file$DCM >= 7])
DEG_CTL_n70 <- names(deg_file$CTL[deg_file$CTL >= 70])

final_features_opt <- unique(c(DEG_ICM_n22, DEG_DCM_n7, DEG_CTL_n70))
length(final_features_opt)


ICM_specific <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM_specific <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))

length(DEG_CTL)
length(DEG_ICM)
length(DEG_DCM)
length(ICM_specific)
length(DCM_specific)

library(ComplexUpset)
library(ggplot2)
all_genes <- union(union(DEG_ICM, DEG_DCM), DEG_CTL)

upset_df <- data.frame(
  gene = all_genes,
  ICM = all_genes %in% DEG_ICM,
  DCM = all_genes %in% DEG_DCM,
  CTL = all_genes %in% DEG_CTL
)

ComplexUpset::upset(
  upset_df,
  intersect = c("ICM", "DCM", "CTL"),
  name = "DEG Overlap",
  base_annotations = list(
    'Intersection size' = intersection_size()
  )
)



load("HF_discovery_cohort_logCPM_for_PCA.RData")
logCPM <- readRDS(file="")

library(FactoMineR)
library(factoextra)
library(umap)
library(ggplot2)
library(dplyr)

final_features <- final_features_opt
expr_pca <- logCPM[rownames(logCPM) %in% final_features, ]
expr_pca_t <- t(expr_pca)

pca_res <- PCA(expr_pca_t, graph = FALSE)

aligned_classes <- classes_df[rownames(expr_pca_t), , drop = FALSE]
group_vector <- as.factor(aligned_classes$Classes)

fviz_pca_ind(
  pca_res,
  geom = "point",
  col.ind = group_vector,
  addEllipses = TRUE,
  palette = c("ICM" = "#d8222b", "DCM" = "#00aacb", "CTL" = "#5d6a85"),
  pointsize = 3,
  legend.title = "Group"
) +
  theme_minimal(base_size = 16) +
  theme(
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black")
  )

fviz_pca_ind(pca_res,
             geom.ind = "point",
             col.ind = group_vector,
             palette = c("ICM" = "#d8222b", "DCM" = "#00aacb", "CTL" = "#5d6a85"),
             addEllipses = TRUE,
             ellipse.type = "t",
             ellipse.level = 0.95,
             ellipse.alpha = 0,
             ellipse.border.size = 1.2,
             legend.title = "Group",
             pointsize = 3) +
  theme_minimal(base_size = 16) +
  theme(
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black")
  )


library(ggplot2)
library(FactoMineR)
library(factoextra)
library(dplyr)

coord_df <- as.data.frame(pca_res$ind$coord)
coord_df$SampleID <- rownames(coord_df)
coord_df <- left_join(coord_df, classes_df, by = c("SampleID" = "ID"))

get_hull <- function(df) df[chull(df$Dim.1, df$Dim.2), ]
hull_df <- coord_df %>% group_by(Classes) %>% do(get_hull(.))

centers_df <- coord_df %>%
  group_by(Classes) %>%
  summarise(
    CenterX = mean(Dim.1),
    CenterY = mean(Dim.2)
  )

ggplot(coord_df, aes(x = Dim.1, y = Dim.2, color = Classes)) +
  stat_ellipse(aes(fill = Classes), geom = "polygon", alpha = 0.2, level = 0.65, type = "norm", color = NA) +
  
  geom_point(size = 5, alpha = 0.9) +
  
  geom_point(data = centers_df, aes(x = CenterX, y = CenterY, fill = Classes),
             shape = 21, size = 12, stroke = 1.2, color = "black", inherit.aes = FALSE) +
  
  geom_text(data = centers_df, aes(x = CenterX, y = CenterY, label = Classes),
            color = "black", fontface = "bold", vjust = -1.8, size = 6, inherit.aes = FALSE) +
  
  scale_color_manual(values = c("ICM" = "#d8222b", "DCM" = "#00aacb", "CTL" = "#5d6a85")) +
  scale_fill_manual(values = c("ICM" = "#d8222b", "DCM" = "#00aacb", "CTL" = "#5d6a85")) +
  
  labs(title = "PCA with Confidence Ellipses and Subtype Centroids") +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 40, color = "black"),
    axis.title = element_text(size = 40),
    legend.title = element_blank(),
    legend.text = element_text(size = 40)
  )

library(ggplot2)
library(dplyr)


coord_df <- coord_df %>%
  group_by(Classes) %>%
  mutate(
    CenterX = mean(Dim.1),
    CenterY = mean(Dim.2),
    Dist2Center = sqrt((Dim.1 - CenterX)^2 + (Dim.2 - CenterY)^2)
  ) %>%
  ungroup()

filtered_df <- coord_df %>%
  group_by(Classes) %>%
  filter(Dist2Center <= quantile(Dist2Center, 0.9)) %>%
  ungroup()

hull_df <- filtered_df %>%
  group_by(Classes) %>%
  group_split() %>%
  lapply(function(df) df[chull(df$Dim.1, df$Dim.2), ]) %>%
  bind_rows()

centers_df <- coord_df %>%
  group_by(Classes) %>%
  summarize(CenterX = mean(Dim.1), CenterY = mean(Dim.2))

colors <- c("ICM" = "#d8222b", "DCM" = "#000080", "CTL" = "#5d6a85")

ggplot(coord_df, aes(x = Dim.1, y = Dim.2, color = Classes)) +
  geom_polygon(data = hull_df, aes(fill = Classes, group = Classes), alpha = 0.25, color = NA) +
  geom_path(data = hull_df, aes(group = Classes), color = "black", size = 0.5) +
  geom_point(size = 4, alpha = 0.85) +
  geom_point(data = centers_df, aes(x = CenterX, y = CenterY, fill = Classes),
             shape = 21, size = 8, stroke = 1.2, color = "black", show.legend = FALSE) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Dim1 (22.5%)",
    y = "Dim2 (10.4%)",
    title = "PCA with Compact Convex Hull"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    axis.text = element_text(size = 40, color = "black"),
    axis.title = element_text(size = 40),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 40)
  )
