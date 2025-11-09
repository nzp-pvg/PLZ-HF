# === one_vs_each_HF.R ===

SplitkFold <- function(Mat, Classes, K){
  require(dplyr)
  require(caret)
  df <- data.frame(ID = colnames(Mat), Classes = Classes)
  samples <- createFolds(df$Classes, k = K, list = TRUE)
  return(list(df = df, samples = samples))
}

OnevsEach.HF <- function(Mat, classes.df, Indices, nDEG, subclass, cv_control) {
  require(limma)
  require(edgeR)
  require(caret)
  require(randomForest)
  require(dplyr)

  # Split into train/test
  TestData <- Mat[, Indices]
  TestPheno <- classes.df[Indices, ]
  TrainData <- Mat[, !colnames(Mat) %in% TestPheno$ID]
  TrainPheno <- classes.df %>% filter(!ID %in% TestPheno$ID)

  # Create DGEList and remove low-expressed genes
  dge <- DGEList(counts = TrainData)
  dge <- calcNormFactors(dge)

  # Construct group: subclass vs Others
  group <- factor(ifelse(TrainPheno$Classes == subclass, subclass, "Others"))

  # Design matrix with group and batch
  design <- model.matrix(~ 0 + group + batch, data = TrainPheno)
  colnames(design) <- make.names(colnames(design))

  # voom transformation
  v <- voom(dge, design, plot = FALSE)

  # Fit model with contrast
  fit <- lmFit(v, design)
  contrast_str <- paste0("group", subclass, " - groupOthers")
  cont <- makeContrasts(contrasts = contrast_str, levels = design)

  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2, trend = TRUE)

  # Get DEGs
  degs <- topTable(fit2, number = Inf)
  degs <- degs[order(degs$t), ]
  half <- floor(nDEG / 2)
  top_features <- rownames(rbind(degs[1:half, ], degs[(nrow(degs)-half+1):nrow(degs), ]))

  # Train classifier
  NewAnn <- ifelse(TrainPheno$Classes == subclass, "One", "Others")
  mtry.val <- length(top_features) / 3
  rf_model <- train(
    x = t(TrainData[top_features, ]),
    y = factor(NewAnn),
    method = "rf",
    trControl = cv_control,
    tuneGrid = expand.grid(.mtry = mtry.val),
    metric = "Kappa"
  )

  # Predict on test data
  probs <- predict(rf_model, newdata = t(TestData), type = "prob") %>% data.frame
  probs$ActualClass <- TestPheno$Classes
  probs$PredictedClass <- predict(rf_model, newdata = t(TestData), type = "raw")

  return(list(Model = rf_model, TestPred = probs, DEGs = top_features))
}