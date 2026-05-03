# Manuscript outputs:
#   - Discovery-stage one-vs-rest RF training
#   - Inputs for Fig. S1C and downstream Fig. S1D-I analyses
# Source basis:
#   - data/PLZ-HF-main/Script/HF_subtyping_v5.R
# === HF_subtyping_v3.R ===
rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)
library(limma)
library(caret)
library(randomForest)
library(glmnet)
library(doParallel)
get_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

resolve_project_root <- function() {
  script_dir <- get_script_dir()
  cands <- unique(c(Sys.getenv("CVD_MS1_PROJECT_ROOT", unset = ""), script_dir, getwd()))
  cands <- cands[nzchar(cands)]
  for (seed in cands) {
    cur <- normalizePath(seed, winslash = "/", mustWork = FALSE)
    for (i in 0:8) {
      probe <- normalizePath(file.path(cur, paste(rep("..", i), collapse = "/")), winslash = "/", mustWork = FALSE)
      if (dir.exists(file.path(probe, "data", "RNA-Seq", "subtype", "top560"))) return(probe)
    }
  }
  stop("Unable to resolve project root. Set CVD_MS1_PROJECT_ROOT or run from the project tree.")
}

project_root <- resolve_project_root()
script_dir <- get_script_dir()
hf_dir <- file.path(project_root, "data", "RNA-Seq", "subtype", "top560")
Cluster <-makeCluster(12)
Cluster <- registerDoParallel(Cluster)

setwd(hf_dir)


# Load discovery matrix (count_annotated, sample_info_discovery)
load("Discovery_start_data.RData") 
source(file.path(script_dir, "03_one_vs_each_training_function.R"))

# Parameters
set.seed(2025)
k <- 10
repeat_n <- 100
subclasses <- c("ICM", "DCM", "CTL")
cv_ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                        classProbs = TRUE, verboseIter = TRUE, returnData = FALSE,
                        savePredictions = "final")

# Save results container
All.Repeats <- list()

for (r in 1:repeat_n) {
  cat("==== Repeat", r, "====\n")
  folds <- SplitkFold(expr_mat, classes_df$Classes, K = k)
  result_one_repeat <- list()

  for (sub in subclasses) {
    cat("  ---- Subclass:", sub, "----\n")
    sub_results <- list()

    for (i in 1:k) {
      cat("    Fold", i, "...\n")
        nDEG_value <- 560
        fold_result <- OnevsEach.HF(
        Mat = expr_mat,
        classes.df = classes_df,
        Indices = folds$samples[[i]],
        nDEG = nDEG_value,
        subclass = sub,
        cv_control = cv_ctrl
      )
      sub_results[[paste0("Fold", i)]] <- fold_result
    }
    result_one_repeat[[sub]] <- sub_results
  }
  All.Repeats[[paste0("Repeat", r)]] <- result_one_repeat
  saveRDS(result_one_repeat, file = paste0("HF_repeat_", r, ".rds"))
  cat("Saved: HF_repeat_", r, ".rds\n")
}

saveRDS(All.Repeats, file = "HF_OnevsEach_CV_top560_degTable_results.rds")

# Extract and save DEG frequency tables
Features.list <- list()
for (sub in subclasses) {
  all_features <- unlist(
    lapply(All.Repeats, function(rep) {
      lapply(rep[[sub]], function(fold) fold$DEGs)
    })
  )
  Features.list[[sub]] <- sort(table(all_features), decreasing = TRUE)
}
saveRDS(Features.list, "HF_DEG_Frequency_top560_dynamic_degTable.rds")


library(pROC)
library(ggplot2)
library(dplyr)
library(dplyr)

results <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

extract_predictions <- function(results_list, subclass) {
  preds <- list()
  for (r in names(results_list)) {
    for (f in names(results_list[[r]][[subclass]])) {
      df <- results_list[[r]][[subclass]][[f]]$TestPred
      df$Repeat <- r
      df$Fold <- f
      df$Subclass <- subclass
      preds[[paste(r, f, subclass, sep = "_")]] <- df
    }
  }
  bind_rows(preds)
}

subclasses <- names(results[[1]])
all_predictions <- bind_rows(lapply(subclasses, function(sc) extract_predictions(results, sc)))
roc_list <- list()
auc_table <- data.frame()

for (sc in subclasses) {
  df <- all_predictions %>% filter(Subclass == sc)
  
  if (!"One" %in% colnames(df)) {
    warning(paste("No 'One' column in prediction for subclass:", sc))
    next
  }
  
  df$Label <- ifelse(df$ActualClass == sc, 1, 0)
  
  roc_obj <- roc(response = df$Label, predictor = df$One, quiet = TRUE)
  roc_list[[sc]] <- roc_obj
  
  ci_obj <- ci.auc(roc_obj)
  auc_table <- rbind(auc_table, data.frame(
    Subclass = sc,
    AUC = auc(roc_obj),
    CI_lower = ci_obj[1],
    CI_upper = ci_obj[3]
  ))
}

roc_df <- do.call(rbind, lapply(names(roc_list), function(class) {
  r <- roc_list[[class]]
  data.frame(
    Spec = 1 - r$specificities,
    Sens = r$sensitivities,
    Class = class
  )
}))

roc_df$Class <- factor(roc_df$Class, levels = subclasses)

subclass_colors <- c("ICM" = "#E6AB2A", "DCM" = "#941651", "CTL" = "#4582B0")

ggplot(roc_df, aes(x = Spec, y = Sens, color = Class)) +
  geom_line(aes(color = Class), size = 3) +
  scale_color_manual(values = subclass_colors) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80", size = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 48),
    axis.title = element_text(color = "black", size = 48)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves (RF Discovery Set)")

print(auc_table)


roc_list_smooth <- list()
auc_ci_list <- list()

for (sc in subclasses) {
  df <- all_predictions %>% filter(Subclass == sc)
  roc_obj <- pROC::roc(
    response = factor(ifelse(df$ActualClass == sc, sc, "Other")),
    predictor = df$One,
    ci = TRUE,
    smooth = TRUE,
    quiet = TRUE
  )
  roc_list_smooth[[sc]] <- roc_obj
  auc_ci_list[[sc]] <- ci.auc(roc_obj)
}

roc_df <- purrr::map2_df(roc_list_smooth, names(roc_list_smooth), function(roc, label) {
  coords <- coords(roc, ret = c("specificity", "sensitivity"))
  data.frame(Spec = 1 - coords$specificity, Sens = coords$sensitivity, Class = label)
})

subclass_colors <- c("ICM" = "#E6AB2A", "DCM" = "#941651", "CTL" = "#4582B0")

ci_ICM <- ci.auc(roc_list_smooth[["ICM"]])
ci_DCM <- ci.auc(roc_list_smooth[["DCM"]])
ci_CTL <- ci.auc(roc_list_smooth[["CTL"]])

auc_label <- paste0(
  "ICM:", round(ci_ICM[2], 3), " (95%CI ", round(ci_ICM[1], 3), "-", round(ci_ICM[3], 3), ")\n",
  "DCM:", round(ci_DCM[2], 3), " (95%CI ", round(ci_DCM[1], 3), "-", round(ci_DCM[3], 3), ")\n",
  "CTL:", round(ci_CTL[2], 3), " (95%CI ", round(ci_CTL[1], 3), "-", round(ci_CTL[3], 3), ")"
)

roc_df$Class <- factor(roc_df$Class, levels = c("ICM", "DCM", "CTL"))
ggplot(roc_df, aes(x = Spec, y = Sens, color = Class)) +
  geom_line(size = 3) +
  scale_color_manual(values = c("ICM" = "#E6AB2A", "DCM" = "#941651", "CTL" = "#4582B0")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80" , size = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 50),
    axis.title = element_text(color = "black", size = 50)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Smoothed ROC Curves") +
  annotate("text", 
           x = 0.2, y = 0.07, 
           label = auc_label,fontface = "italic",
           size = 12, hjust = 0)


rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)
library(limma)
library(caret)


setwd(hf_dir)
load("Discovery_start_data.RData")
All.Repeats <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

library(dplyr)
library(caret)

test_preds_all <- list()

for (r in names(All.Repeats)) {
  rep_data <- All.Repeats[[r]]
  
  for (fold_id in names(rep_data[[1]])) {
    subclass_preds <- list()
    
    for (sub in c("ICM", "DCM", "CTL")) {
      fold_pred <- rep_data[[sub]][[fold_id]]$TestPred
      if (!is.null(fold_pred)) {
        subclass_preds[[sub]] <- fold_pred %>%
          select(ActualClass, matches("One|Others")) %>%
          rename(!!sub := One) %>%
          mutate(SampleID = rownames(fold_pred))
      }
    }
    
    merged_pred <- Reduce(function(x, y) full_join(x, y, by = "SampleID"), subclass_preds)
    
    merged_pred$TrueClass <- merged_pred$ActualClass.x
    merged_pred <- merged_pred %>%
      select(SampleID, TrueClass, ICM, DCM, CTL)
    
    merged_pred$PredictedClass <- apply(merged_pred[, c("ICM", "DCM", "CTL")], 1, function(x) {
      colnames(merged_pred[, c("ICM", "DCM", "CTL")])[which.max(x)]
    })
    
    test_preds_all[[paste0(r, "_", fold_id)]] <- merged_pred
  }
}

all_predictions <- bind_rows(test_preds_all)

all_predictions$TrueClass <- factor(all_predictions$TrueClass, levels = c("ICM", "DCM", "CTL"))
all_predictions$PredictedClass <- factor(all_predictions$PredictedClass, levels = c("ICM", "DCM", "CTL"))

conf_mat <- confusionMatrix(all_predictions$PredictedClass, all_predictions$TrueClass)
print(conf_mat)
library(openxlsx)

conf_table <- as.data.frame.matrix(conf_mat$table)
conf_overall <- as.data.frame(t(conf_mat$overall))
conf_by_class <- as.data.frame(conf_mat$byClass)

wb <- createWorkbook()

addWorksheet(wb, "Confusion Matrix")
writeData(wb, "Confusion Matrix", conf_table, rowNames = TRUE)

addWorksheet(wb, "Overall Statistics")
writeData(wb, "Overall Statistics", conf_overall, rowNames = FALSE)

addWorksheet(wb, "By Class Statistics")
writeData(wb, "By Class Statistics", conf_by_class, rowNames = TRUE)

saveWorkbook(wb, file = "HF_Confusion_top560_Matrix_Result.xlsx", overwrite = TRUE)

table(classes_df$batch)

## PCA
# === HF_subtyping_v3.R ===
rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)
library(limma)
library(caret)


setwd(hf_dir)
library(edgeR)
library(sva)
library(dplyr)
library(ggplot2)

load("Discovery_start_data.RData")
deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

batch <- classes_df$batch
group <- classes_df$Classes
expr_combat <- ComBat_seq(counts = as.matrix(expr_mat), batch = batch)

dge <- DGEList(counts = expr_combat)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)


get_filtered_deg <- function(res, subclass) {
  unlist(lapply(res, function(rep) {
    lapply(rep[[subclass]], function(fold) {
      degs <- fold$DEG_table
      rownames(degs[which(degs$adj.P.Val < 0.01 & abs(degs$logFC) > 1), ])
    })
  }), use.names = FALSE)
}

DEGs_ICM <- get_filtered_deg(deg_result, "ICM")
DEGs_DCM <- get_filtered_deg(deg_result, "DCM")
DEGs_CTL <- get_filtered_deg(deg_result, "CTL")

DEG_ICM <- names(which(table(DEGs_ICM) >= 22))
DEG_DCM <- names(which(table(DEGs_DCM) >= 7))
DEG_CTL <- names(which(table(DEGs_CTL) >= 70))

DEG_ICM <- names(which(table(DEGs_ICM) >= 100))
DEG_DCM <- names(which(table(DEGs_DCM) >= 100))
DEG_CTL <- names(which(table(DEGs_CTL) >= 100))


ICM_only <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM_only <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL_only <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))
ICM_DCM_shared <- intersect(DEG_ICM, DEG_DCM)
final_features <- unique(c(ICM_only, DCM_only))

save(DEG_CTL,DEG_DCM, DEG_ICM, ICM_only,DCM_only,CTL_only,ICM_DCM_shared,final_features, file= "HF_DEGs.RData")

get_top_features <- function(res, subclass) {
  unlist(lapply(res, function(rep) {
    lapply(rep[[subclass]], function(fold) {
      fold$DEGs
    })
  }), use.names = FALSE)
}

DEGs_ICM <- get_top_features(deg_result, "ICM")
DEGs_DCM <- get_top_features(deg_result, "DCM")
DEGs_CTL <- get_top_features(deg_result, "CTL")

DEG_ICM <- unique(DEGs_ICM)
DEG_DCM <- unique(DEGs_DCM)
DEG_CTL <- unique(DEGs_CTL)

ICM_only <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM_only <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL_only <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))
ICM_DCM_shared <- intersect(DEG_ICM, DEG_DCM)

ICM_DCM_only <- setdiff(intersect(DEG_ICM, DEG_DCM), DEG_CTL)
length(ICM_DCM_only)


length(DEG_ICM)
length(ICM_only)
length(DEG_DCM)
length(DCM_only)
length(DEG_CTL)
length(CTL_only)
length(ICM_DCM_shared)


deg_lists <- list(
  ICM = DEG_ICM,
  DCM = DEG_DCM,
  CTL = DEG_CTL
)
ICM_only <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM_only <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL_only <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))
ICM_DCM_shared <- intersect(DEG_ICM, DEG_DCM)


library(UpSetR)

deg_list <- list(
  ICM = DEG_ICM,
  DCM = DEG_DCM,
  CTL = DEG_CTL
)

input_upset <- fromList(deg_list)

upset(
  input_upset,
  sets = c("ICM", "DCM", "CTL"),
  sets.bar.color = c("#d8222b", "#000080", "#5d5d5d"),
  order.by = "freq",
  mainbar.y.label = "Intersection Size",
  sets.x.label = "DEGs per Group",
  text.scale = 1.5
)
ICM_DCM_CTL_shared <- Reduce(intersect, list(DEG_ICM, DEG_DCM, DEG_CTL))
length(ICM_DCM_CTL_shared)




deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

get_model_features <- function(res, subclass) {
  unlist(lapply(res, function(rep) {
    lapply(rep[[subclass]], function(fold) fold$DEGs)
  }), use.names = FALSE)
}

DEGs_ICM <- get_model_features(deg_result, "ICM")
DEGs_DCM <- get_model_features(deg_result, "DCM")
DEGs_CTL <- get_model_features(deg_result, "CTL")

DEG_ICM <- names(which(table(DEGs_ICM) >= 22))
DEG_DCM <- names(which(table(DEGs_DCM) >= 7))
DEG_CTL <- names(which(table(DEGs_CTL) >= 70))

DEG_ICM <- names(which(table(DEGs_ICM) >= 100))
DEG_DCM <- names(which(table(DEGs_DCM) >= 100))
DEG_CTL <- names(which(table(DEGs_CTL) >= 100))



ICM_only <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM_only <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL_only <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))

three_class_DEGs <- unique(c(ICM_only, DCM_only, CTL_only))


cat("ICM_only:", length(ICM_only), "\n")
cat("DCM_only:", length(DCM_only), "\n")
cat("ICM_DCM_shared:", length(ICM_DCM_shared), "\n")
length(DEG_ICM)
length(DEG_DCM)
length(DEG_CTL)
length(DCM_only)
length(ICM_only)
length(CTL_only)
length(ICM_DCM_shared)
length(HF_DEGs)
length(HF_specific)

load("HF_discovery_cohort_logCPM_for_PCA.RData")

library(dplyr)
library(FactoMineR)
library(factoextra)

ctl_ids <- classes_df %>% filter(Classes == "CTL") %>% pull(ID)
expr_ctl <- logCPM[, ctl_ids]  # ENSURE: colnames(expr_ctl) = Sample IDs

pca_res <- PCA(t(expr_ctl), graph = FALSE)

pca_coords <- as.data.frame(pca_res$ind$coord[, 1:2])
pca_coords$SampleID <- rownames(pca_coords)

ctl_center <- colMeans(pca_coords[, 1:2])
pca_coords$Dist <- sqrt((pca_coords$Dim.1 - ctl_center[1])^2 + (pca_coords$Dim.2 - ctl_center[2])^2)

library(ggplot2)
library(ggforce)
ggplot(pca_coords, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Dist), size = 4) +
  geom_text(aes(label = SampleID), vjust = -0.5, size = 4) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw()

library(dplyr)
library(ggplot2)
library(umap)

# ================================
# ================================
expr_pca <- logCPM[rownames(logCPM) %in% three_class_DEGs, ]
expr_pca_t <- t(expr_pca)

umap_df <- data.frame(SampleID = rownames(expr_pca_t)) %>%
  left_join(classes_df, by = c("SampleID" = "ID"))

# ================================
# ================================
set.seed(2025)
umap_result <- umap(expr_pca_t)
umap_df$UMAP1 <- umap_result$layout[, 1]
umap_df$UMAP2 <- umap_result$layout[, 2]

# ================================
# ================================
centers <- umap_df %>%
  group_by(Classes) %>%
  summarise(Center_UMAP1 = mean(UMAP1), Center_UMAP2 = mean(UMAP2))

umap_df <- umap_df %>%
  left_join(centers, by = "Classes") %>%
  mutate(Dist = sqrt((UMAP1 - Center_UMAP1)^2 + (UMAP2 - Center_UMAP2)^2))

# ================================
# ================================
remove_top_n_list <- list(
  CTL = 19,
  DCM = 1,
  ICM = 1
)

to_remove <- c()

for (cls in names(remove_top_n_list)) {
  n_remove <- remove_top_n_list[[cls]]
  if (n_remove > 0) {
    remove_ids <- umap_df %>%
      filter(Classes == cls) %>%
      slice_max(Dist, n = n_remove) %>%
      pull(SampleID)
    to_remove <- c(to_remove, remove_ids)
  }
}

# ================================
# ================================
logCPM_filtered <- logCPM[, !colnames(logCPM) %in% to_remove]
classes_df_filtered <- classes_df %>% filter(!(ID %in% to_remove))

expr_pca_t_filtered <- t(logCPM_filtered[rownames(logCPM_filtered) %in% three_class_DEGs, ])

pca_res <- prcomp(expr_pca_t_filtered, center = TRUE, scale. = TRUE)

umap_input <- pca_res$x[, 1:20]

set.seed(2025)
umap_result_filtered <- umap(umap_input)

umap_plot_df <- data.frame(
  SampleID = rownames(umap_input),
  UMAP1 = umap_result_filtered$layout[, 1],
  UMAP2 = umap_result_filtered$layout[, 2]
) %>%
  left_join(classes_df_filtered, by = c("SampleID" = "ID"))
# ================================
# ================================
centroids <- umap_plot_df %>%
  group_by(Classes) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

subclass_colors <- c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")

umap_var <- apply(umap_result_filtered$layout, 2, var)
umap_var_exp <- round(100 * umap_var / sum(umap_var), 1)
x_lab <- paste0("UMAP1 (", umap_var_exp[1], "%)")
y_lab <- paste0("UMAP2 (", umap_var_exp[2], "%)")

ggplot(umap_plot_df, aes(x = UMAP1, y = UMAP2, color = Classes)) +
  stat_density_2d(aes(fill = Classes), geom = "polygon", contour = TRUE, alpha = 0.15, color = NA) +
  
  geom_point(size = 7, stroke = 0.5, alpha = 0.8, shape = 21, aes(fill = Classes), color = "black") +
  
  geom_polygon(data = umap_plot_df %>% group_by(Classes) %>% slice(chull(UMAP1, UMAP2)),
               aes(fill = Classes), alpha = 0.15, color = "black", linewidth = 0, show.legend = FALSE) +
  
  geom_point(data = centroids, aes(x = UMAP1, y = UMAP2, fill = Classes), 
             shape = 21, size = 18, alpha = 0.5, stroke = 1.2, color = "black") +
  
  scale_fill_manual(values = subclass_colors) +
  
  labs(
    title = "UMAP with Subtype Density and Hulls",
    x = "UMAP1 (33%)", y = "UMAP2 (67%)"
  ) +
  
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
    axis.text = element_text(color = "black", size = 46),
    axis.title = element_text(color = "black", size = 46),
    legend.title = element_blank(),
    legend.text = element_text(size = 46)
  )


# ================================
# PCA : ICM vs DCM 
# ================================
load("Discovery_start_data.RData")
deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

library(sva)
library(edgeR)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)

# ================================
# ================================
batch <- classes_df$batch
expr_combat <- ComBat_seq(counts = as.matrix(expr_mat), batch = batch)

dge <- DGEList(counts = expr_combat)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)

# ================================
# ================================
get_model_features <- function(res, subclass) {
  unlist(lapply(res, function(rep) {
    lapply(rep[[subclass]], function(fold) fold$DEGs)
  }), use.names = FALSE)
}

DEGs_ICM <- get_model_features(deg_result, "ICM")
DEGs_DCM <- get_model_features(deg_result, "DCM")
DEGs_CTL <- get_model_features(deg_result, "CTL")

DEG_ICM <- names(which(table(DEGs_ICM) >= 100))
DEG_DCM <- names(which(table(DEGs_DCM) >= 100))
DEG_CTL <- names(which(table(DEGs_CTL) >= 100))

ICM_only <- setdiff(DEG_ICM, union(DEG_DCM, DEG_CTL))
DCM_only <- setdiff(DEG_DCM, union(DEG_ICM, DEG_CTL))
CTL_only <- setdiff(DEG_CTL, union(DEG_ICM, DEG_DCM))
selected_genes <- unique(c(ICM_only, DCM_only, CTL_only))

# ================================
# ================================
sel_samples <- classes_df$Classes %in% c("ICM", "DCM")
expr_pca <- logCPM[rownames(logCPM) %in% selected_genes, classes_df$ID[sel_samples]]
expr_pca_t <- t(expr_pca)

pca_input_df <- data.frame(SampleID = rownames(expr_pca_t)) %>%
  left_join(classes_df[, c("ID", "Classes")], by = c("SampleID" = "ID")) %>%
  rename(Class = Classes)

# ================================
# ================================
res.pca <- PCA(expr_pca_t, graph = FALSE)
pca_scores <- as.data.frame(res.pca$ind$coord)
pca_scores$Class <- pca_input_df$Class
pca_scores$SampleID <- pca_input_df$SampleID

centers <- pca_scores %>%
  group_by(Class) %>%
  summarise(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2), .groups = "drop")

pc1_var <- round(res.pca$eig[1, 2], 1)
pc2_var <- round(res.pca$eig[2, 2], 1)

# ================================
# ================================
hulls <- pca_scores %>%
  group_by(Class) %>%
  slice(chull(Dim.1, Dim.2))

geom_polygon(
  data = hulls,
  aes(x = Dim.1, y = Dim.2, fill = Class),
  alpha = 0.15,
  color = "black",
  linewidth = 1.2,
  show.legend = FALSE
)
ggplot(pca_scores, aes(x = Dim.1, y = Dim.2)) +
  stat_density_2d(aes(fill = Class), geom = "polygon", alpha = 0.15, color = NA) +

  geom_polygon(data = hulls, aes(x = Dim.1, y = Dim.2, fill = Class),
               alpha = 0.15, color = "black", linewidth = 0.3, show.legend = FALSE) +

  geom_point(aes(fill = Class), shape = 21, size = 8, stroke = 1, alpha = 0.8, color = "black") +

  geom_point(data = centers, aes(x = Dim.1, y = Dim.2, fill = Class),
             shape = 21, size = 18, stroke = 1.2, alpha = 0.8, color = "black") +
  geom_segment(
    data = data.frame(
      x = centers$Dim.1[centers$Class == "ICM"],
      y = centers$Dim.2[centers$Class == "ICM"],
      xend = centers$Dim.1[centers$Class == "DCM"],
      yend = centers$Dim.2[centers$Class == "DCM"]
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black", linetype = "dashed", linewidth = 1.2
  )+
  scale_fill_manual(values = c("DCM" = "#941651", "ICM" = "#E6AB2A")) +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal(base_size = 48) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.9),
    axis.text = element_text(size = 46, color = "black"),
    axis.title = element_text(size = 48, color = "black"),
    axis.title.x = element_text(margin = unit(c(0, 0, 0, 0), "pt")),
    axis.title.y = element_text(margin = unit(c(0, -10, 0, 0), "pt")),
    legend.title = element_blank(),
    legend.text = element_text(size = 46),
    plot.title = element_text(size = 54, face = "bold", hjust = 0.5)
  )

library(cluster)
pca_input <- pca_scores %>%
  select(Dim.1, Dim.2) %>%
  as.matrix()

class_labels <- pca_scores$Class  # "ICM", "DCM"

sil <- silhouette(as.numeric(factor(class_labels)), dist(pca_input))

mean_sil <- mean(sil[, "sil_width"])
print(paste("Average silhouette score:", round(mean_sil, 3)))



rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(randomForest)
library(viridis)

setwd(hf_dir)

load("HF_discovery_cohort_logCPM_for_PCA.RData")
deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

# ================================
# ================================
get_model_features <- function(res, subclass) {
  unlist(lapply(res, function(rep) {
    lapply(rep[[subclass]], function(fold) fold$DEGs)
  }), use.names = FALSE)
}

get_model_tables <- function(res, subclass) {
  do.call(rbind, lapply(res, function(rep) {
    do.call(rbind, lapply(rep[[subclass]], function(fold) fold$DEG_table))
  }))
}

DEGs_ICM <- get_model_features(deg_result, "ICM")
DEGs_DCM <- get_model_features(deg_result, "DCM")
DEGs_CTL <- get_model_features(deg_result, "CTL")

DEG_ICM <- names(which(table(DEGs_ICM) >= 1))
DEG_DCM <- names(which(table(DEGs_DCM) >= 1))
DEG_CTL <- names(which(table(DEGs_CTL) >= 1))

# ================================
# ================================

merge_and_filter_deg_tables <- function(deg_result, subclass, logfc_threshold = 1, padj_threshold = 0.01) {
  all_degs <- do.call(rbind, lapply(deg_result, function(rep) {
    do.call(rbind, lapply(rep[[subclass]], function(fold) {
      degs <- fold$DEG_table
      degs$raw_id <- rownames(degs)
      return(degs)
    }))
  }))
  
  all_degs$gene_id <- sub(".*\\.", "", all_degs$raw_id)
  all_degs <- all_degs[all_degs$adj.P.Val < padj_threshold & abs(all_degs$logFC) > logfc_threshold, ]
  all_degs <- all_degs[order(abs(all_degs$logFC), decreasing = TRUE), ]
  all_degs <- all_degs[!duplicated(all_degs$gene_id), ]
  
  rownames(all_degs) <- all_degs$gene_id
  all_degs$gene_id <- NULL
  all_degs$raw_id <- NULL
  return(all_degs)
}

deg_tbl_ICM <- merge_and_filter_deg_tables(deg_result, "ICM")
deg_tbl_DCM <- merge_and_filter_deg_tables(deg_result, "DCM")
deg_tbl_CTL <- merge_and_filter_deg_tables(deg_result, "CTL")

top_n <- 100
top_ICM <- rownames(deg_tbl_ICM[order(-abs(deg_tbl_ICM$logFC)), ])[1:min(top_n, nrow(deg_tbl_ICM))]
top_DCM <- rownames(deg_tbl_DCM[order(-abs(deg_tbl_DCM$logFC)), ])[1:min(top_n, nrow(deg_tbl_DCM))]
top_CTL <- rownames(deg_tbl_CTL[order(-abs(deg_tbl_CTL$logFC)), ])[1:min(top_n, nrow(deg_tbl_CTL))]

# ================================
# ================================
selected_genes <- unique(c(top_ICM, top_DCM, top_CTL))
selected_genes <- selected_genes[selected_genes %in% rownames(logCPM)]

expr_mat <- logCPM[selected_genes, ]
cv_values <- apply(expr_mat, 1, function(x) sd(x) / (mean(x) + 1e-5))
selected_genes_cv_filtered <- names(cv_values[cv_values > 0.4])
final_genes <- intersect(selected_genes, selected_genes_cv_filtered)

# ================================
# ================================
heat_expr <- logCPM[final_genes, classes_df$ID]
zscore <- t(scale(t(heat_expr)))

annotation_col <- data.frame(Group = classes_df$Classes)
rownames(annotation_col) <- classes_df$ID

gene_category <- rep("Other", length(final_genes))
gene_category[final_genes %in% top_ICM] <- "ICM_top"
gene_category[final_genes %in% top_DCM] <- "DCM_top"
gene_category[final_genes %in% top_CTL] <- "CTL_top"
gene_category <- factor(gene_category, levels = c("ICM_top", "DCM_top", "CTL_top"))

annotation_row <- data.frame(DEG_Category = gene_category)
rownames(annotation_row) <- final_genes

# ================================
# ================================
ann_colors <- list(
  Group = c("ICM" = "#E6AB2A", "DCM" = "#941651", "CTL" = "#1F78B4"),
  DEG_Category = c("ICM_top" = "#E6AB2A", "DCM_top" = "#941651", "CTL_top" = "#4A90E2")
)

sorted_ids <- classes_df %>%
  arrange(factor(Classes, levels = c("ICM", "DCM", "CTL"))) %>%
  pull(ID)

zscore_sorted <- zscore[, sorted_ids]
annotation_col_sorted <- annotation_col[sorted_ids, , drop = FALSE]

annotation_row_sorted <- annotation_row[rownames(zscore_sorted), , drop = FALSE]

# ================================
# ===============================
color_palette <- colorRampPalette(c("#01C4C4", "black", "gold"))(255)
breaks_seq <- seq(-1.5, 1.5, length.out = 255)
n_ICM <- sum(annotation_col_sorted$Group == "ICM")
n_DCM <- sum(annotation_col_sorted$Group == "DCM")

library(ComplexHeatmap)
library(circlize)
library(grid)

column_split = annotation_col_sorted$Group

ha_col <- HeatmapAnnotation(
  Group = annotation_col_sorted$Group,
  col = list(Group = ann_colors$Group)
)

top_ICM_filtered <- intersect(DEG_ICM, top_genes_group$ICM)
top_DCM_filtered <- intersect(DEG_DCM, top_genes_group$DCM)
top_CTL_filtered <- intersect(DEG_CTL, top_genes_group$CTL)

deg_category <- sapply(rownames(zscore_sorted), function(gene) {
  if (gene %in% top_ICM_filtered) {
    return("ICM_top")
  } else if (gene %in% top_DCM_filtered) {
    return("DCM_top")
  } else if (gene %in% top_CTL_filtered) {
    return("CTL_top")
  } else {
    return("Other")
  }
})

deg_colors <- c("ICM_top" = "#E6AB2A", "DCM_top" = "#941651", "CTL_top" = "#1F78B4")

row_anno <- rowAnnotation(
  DEG_Category = annotation_row$DEG_Category,
  col = list(DEG_Category = deg_colors),
  width = unit(5, "mm")
)



Heatmap(
  zscore_sorted,
  name = "Z-score",
  col = colorRamp2(c(-1, 0, 1), c("#00b9b9", "black",  "#DE4967")),
  #col = colorRamp2(c(-1, 0, 1), c( "black", "#CCC19C", "#DE4967")),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = as.dist(1 - cor(t(zscore_sorted), method = "spearman")),
  #clustering_distance_rows = as.dist(1 - cor(t(zscore_sorted), method = "pearson")),
  #clustering_distance_rows = as.dist(1 - cor(t(zscore_sorted), method = "kendall")),
  clustering_method_rows = "ward.D2",
  #clustering_method_rows = "complete",
  #clustering_method_rows = "average",
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = annotation_col_sorted$Group,
  top_annotation = HeatmapAnnotation(
    Group = annotation_col_sorted$Group,
    col = list(Group = ann_colors$Group),
    annotation_name_side = "right",
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 10)
    ),
    gp = gpar(col = NA),
    simple_anno_size = unit(20, "mm")
  ),
  row_dend_width = unit(20, "mm"),
  row_dend_gp = gpar(lwd = 3),
  heatmap_legend_param = list(title = "Z-score"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "black", lwd = 0.3, fill = NA))
  }
)



load("HF_discovery_cohort_logCPM_for_PCA.RData")  # logCPM, classes_df

group_list <- split(classes_df$ID, classes_df$Classes)
logCPM_grouped <- lapply(group_list, function(ids) logCPM[, ids])

topN <- 50
top_genes_group <- lapply(logCPM_grouped, function(mat) {
  gene_sd <- apply(mat, 1, sd)
  names(sort(gene_sd, decreasing = TRUE))[1:topN]
})

selected_genes <- unique(unlist(top_genes_group))
selected_genes <- selected_genes[selected_genes %in% rownames(logCPM)]

heat_expr <- logCPM[selected_genes, classes_df$ID]
zscore <- t(scale(t(heat_expr)))

library(dplyr)
sorted_ids <- classes_df %>%
  arrange(factor(Classes, levels = c("ICM", "DCM", "CTL"))) %>%
  pull(ID)
zscore_sorted <- zscore[, sorted_ids]

annotation_col <- data.frame(Group = classes_df$Classes)
rownames(annotation_col) <- classes_df$ID
annotation_col_sorted <- annotation_col[sorted_ids, , drop = FALSE]

deg_category <- sapply(rownames(zscore_sorted), function(g) {
  if (g %in% top_genes_group$ICM) {
    "ICM_top"
  } else if (g %in% top_genes_group$DCM) {
    "DCM_top"
  } else if (g %in% top_genes_group$CTL) {
    "CTL_top"
  } else {
    "Other"
  }
})
annotation_row <- data.frame(DEG_Category = deg_category)
rownames(annotation_row) <- rownames(zscore_sorted)

ann_colors <- list(
  Group = c("ICM" = "#E6AB2A", "DCM" = "#941651", "CTL" = "#1F78B4"),
  DEG_Category = c("ICM_top" = "#E6AB2A", "DCM_top" = "#941651", "CTL_top" = "#4A90E2")
)
color_palette <- colorRampPalette(c("blue", "black", "yellow"))(255)
breaks_seq <- seq(-2.5, 2.5, length.out = 255)

n_ICM <- sum(annotation_col_sorted$Group == "ICM")
n_DCM <- sum(annotation_col_sorted$Group == "DCM")

library(pheatmap)
pheatmap(
  zscore_sorted,
  annotation_col = annotation_col_sorted,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "correlation",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 10,
  scale = "none",
  color = color_palette,
  breaks = breaks_seq,
  gaps_col = c(n_ICM, n_ICM + n_DCM),
  main = "Top Variable Genes by Group (ICM / DCM / CTL)"
)

# Notes retained from the original exploratory plotting block:
# clustering_distance_rows = "correlation"
# clustering_method = "complete"



## HF Disease Index
library(writexl)
library(dplyr)


extract_disease_index <- function(deg_result, subtype = c("ICM", "DCM")) {
  subtype <- match.arg(subtype)
  prob_list <- lapply(deg_result, function(rep) {
    lapply(rep[[subtype]], function(fold) {
      df <- fold$TestPred
      df$Sample_ID <- rownames(df)
      df[, c("Sample_ID", "One", "ActualClass")]
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  prob_summary <- prob_list %>%
    group_by(Sample_ID) %>%
    summarise(
      Index = mean(One),
      True_Label = dplyr::first(ActualClass)
    ) %>%
    ungroup()
  
  colnames(prob_summary)[2] <- paste0(subtype, "_index")
  return(prob_summary)
}

ICM_index_df <- extract_disease_index(deg_result, "ICM")
DCM_index_df <- extract_disease_index(deg_result, "DCM")

merged_df <- full_join(ICM_index_df, DCM_index_df, by = c("Sample_ID", "True_Label"))

write_xlsx(hf_index_df, path = "HF_index_from_binary_models.xlsx")

library(ggplot2)

ICM_index_df <- extract_disease_index_2step(deg_result, subclass = "ICM")
DCM_index_df <- extract_disease_index_2step(deg_result, subclass = "DCM")

hf_index_all_df <- ICM_index_df %>%
  rename(ICM_index = Mean) %>%
  inner_join(DCM_index_df %>% rename(DCM_index = Mean), by = c("Sample_ID", "ActualClass"))

hf_index_all_df$HF_index <- rowMeans(hf_index_all_df[, c("ICM_index", "DCM_index")])


library(ggplot2)

hf_index_all_df$True_Label <- hf_index_all_df$ActualClass

ggplot(hf_index_all_df, aes(x = True_Label, y = HF_index, fill = True_Label)) +
  geom_boxplot(alpha = 0.3, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")) +
  scale_color_manual(values = c("DCM" = "#941651", "CTL" = "#4582B0", "ICM" = "#E6AB2A")) +
  labs(
    title = "HF Disease Index",
    x = NULL,
    y = "HF_index"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )

writexl::write_xlsx(hf_index_all_df[, c("Sample_ID", "ActualClass", "ICM_index", "DCM_index", "HF_index")], "HF_Disease_Index.xlsx")

# batch effect reducing

# === HF_subtyping_v3.R ===
rm(list=ls())
gc()
library(dplyr)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(tidyr)
library(edgeR)
library(limma)
library(sva)
library(ggplot2)

setwd(hf_dir)


# Load discovery matrix (count_annotated, sample_info_discovery)
load("Discovery_start_data.RData") 
expr_mat_raw <- expr_mat

batch <- classes_df$batch
group <- classes_df$Classes
expr_combat <- ComBat_seq(counts = as.matrix(expr_mat), batch = batch)

dge_raw <- DGEList(counts = expr_mat_raw)
dge_raw <- calcNormFactors(dge_raw)
logCPM_raw <- cpm(dge_raw, log = TRUE)

dge_combat <- DGEList(counts = expr_combat)
dge_combat <- calcNormFactors(dge_combat)
logCPM_combat <- cpm(dge_combat, log = TRUE)

deg_result <- readRDS("HF_OnevsEach_CV_top560_degTable_results.rds")

get_model_features <- function(res, subclass) {
  unlist(lapply(res, function(rep) {
    lapply(rep[[subclass]], function(fold) fold$DEGs)
  }), use.names = FALSE)
}

top_var_genes <- names(sort(apply(logCPM_raw, 1, mad), decreasing = TRUE))[1:5000]
expr_raw <- logCPM_raw[top_var_genes, ]
expr_combat <- logCPM_combat[top_var_genes, ]

library(ggplot2)
library(dplyr)

pca_raw <- prcomp(t(expr_raw), scale. = TRUE)
pca_combat <- prcomp(t(expr_combat), scale. = TRUE)

pca_df_raw <- data.frame(
  ID = colnames(expr_raw),
  PC1 = pca_raw$x[, 1],
  PC2 = pca_raw$x[, 2],
  BatchStatus = "before BER"
)

pca_df_combat <- data.frame(
  ID = colnames(expr_combat),
  PC1 = pca_combat$x[, 1],
  PC2 = pca_combat$x[, 2],
  BatchStatus = "after BER"
)

pca_df <- rbind(pca_df_raw, pca_df_combat)

pca_df <- left_join(pca_df, classes_df, by = "ID")

pca_df$BatchStatus <- factor(pca_df$BatchStatus, levels = c("before BER", "after BER"))

batch_colors <- c(
  "B69" = "#001f5b",
  "B50" = "#0096a0",
  "B66" = "#6a2ca0",
  "B55" = "#009246",
  "B52" = "#941651",
  "B24" = "#C55A11"
)

expl_var <- 100 * (pca_raw$sdev)^2 / sum((pca_raw$sdev)^2)
x_lab <- paste0("PC1 (", round(expl_var[1], 1), "%)")
y_lab <- paste0("PC2 (", round(expl_var[2], 1), "%)")

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = batch), shape = 21, size = 6, color = "black", stroke = 0.6, alpha = 0.6) +
  facet_wrap(~BatchStatus, dir = "h") +
  scale_fill_manual(values = batch_colors) +
  theme_bw(base_size = 16) +
  labs(
    x = x_lab,
    y = y_lab,
    title = "PCA Before and After Batch Effect Correction",
    fill = "Batch"
  ) +
  theme(
    strip.text = element_text(size = 40),
    axis.text = element_text(size = 40),
    axis.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 40),
    panel.border = element_rect(color = "black", size = 1)
  )
