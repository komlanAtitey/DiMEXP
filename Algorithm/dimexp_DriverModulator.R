library(Matrix)
library(MASS)  # for ginv
library(expm)  # only for matrix powers / expm, not logm
library(ggplot2)
library(reshape2)
library(vars)
library(parallel)
library(tidyverse)
library(zoo)
library(dplyr)
library(igraph)
library(ggraph)
library(visNetwork)
library(scales)

# ------------------------------------------------------
# TIME-VARYING GENE INTERACTION NETWORKS (EDGE DEFINITION)
# ------------------------------------------------------
# Granger causality score function
compute_granger_score_safe <- function(x, y, max_lag = 2) {
  df <- data.frame(x = x, y = y)
  model <- try(VAR(df, p = max_lag, type = "const"), silent = TRUE)
  if (inherits(model, "try-error")) return(0)
  test <- try(causality(model, cause = "x"), silent = TRUE)
  if (inherits(test, "try-error")) return(0)
  -log10(test$Granger$p.value + 1e-10)
}

# ------------------------------------------------------
# CAUSAL INFLUENCE SCORE (CIS)
# ------------------------------------------------------
# Compute CIS matrix between non-drivers and drivers
compute_CIS_matrix_parallel <- function(expr_mat, non_drivers, drivers, max_lag = 2, ncores = 4) {
  # Make rownames unique
  rownames(expr_mat) <- make.unique(rownames(expr_mat))
  
  # Clean gene lists
  non_drivers <- intersect(non_drivers, rownames(expr_mat))
  drivers <- intersect(drivers, rownames(expr_mat))
  
  # Setup cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("expr_mat", "drivers", "max_lag", "compute_granger_score_safe"), envir = environment())
  clusterEvalQ(cl, library(vars))
  
  # Function to compute one row
  compute_row <- function(gene_i) {
    tryCatch({
      x <- as.numeric(expr_mat[gene_i, ])
      sapply(drivers, function(gene_j) {
        y <- as.numeric(expr_mat[gene_j, ])
        compute_granger_score_safe(x, y, max_lag)
      })
    }, error = function(e) rep(NA, length(drivers)))
  }
  
  # Parallel computation
  CIS_list <- parLapply(cl, non_drivers, compute_row)
  stopCluster(cl)
  
  CIS_matrix <- do.call(rbind, CIS_list)
  rownames(CIS_matrix) <- non_drivers
  colnames(CIS_matrix) <- drivers
  return(CIS_matrix)
}

# Identify top-k non-driver influencers per driver gene
get_top_influencers_per_driver <- function(CIS_matrix, top_n = 5) {
  # Remove rows with NA rownames or NA influence scores
  valid_rows <- !is.na(rownames(CIS_matrix)) & rowSums(is.na(CIS_matrix)) < ncol(CIS_matrix)
  CIS_matrix <- CIS_matrix[valid_rows, , drop = FALSE]
  
  results <- list()
  
  for (driver in colnames(CIS_matrix)) {
    top_df <- data.frame(
      Driver_Gene = driver,
      Non_Driver_Gene = rownames(CIS_matrix),
      Influence_Score = CIS_matrix[, driver],
      stringsAsFactors = FALSE
    ) %>%
      arrange(desc(Influence_Score)) %>%
      slice_head(n = top_n)
    
    results[[driver]] <- top_df
  }
  
  return(results)
}


# ------------------------------------------------------
# TRANSFORMATION LIKELIHOOD FUNCTION
# ------------------------------------------------------
# Fit logistic regression to predict transformation likelihood
compute_transformation_likelihood <- function(CIS_matrix, driver_genes, all_genes) {
  # Ensure unique rownames
  rownames(CIS_matrix) <- make.unique(rownames(CIS_matrix))
  
  # Initialize influence scores
  influence_score <- rep(0, length(all_genes))
  names(influence_score) <- all_genes
  
  # Map CIS_matrix rows to all_genes
  idx <- match(rownames(CIS_matrix), all_genes)
  valid <- !is.na(idx)
  if (any(valid)) {
    influence_score[idx[valid]] <- rowSums(CIS_matrix[valid, , drop = FALSE])
  }
  influence_score[is.na(influence_score)] <- 0
  
  # Label driver genes
  label <- ifelse(all_genes %in% driver_genes, 1, 0)
  
  # Build model dataframe
  model_df <- data.frame(
    Gene = all_genes,
    Label = label,
    Influence = influence_score,
    stringsAsFactors = FALSE
  )
  
  # Fit logistic regression
  fit <- glm(Label ~ Influence, data = model_df, family = "binomial")
  
  # Predict safely
  model_df$Transformation_Likelihood <- predict(fit, newdata = model_df, type = "response")
  
  return(list(model_df = model_df, model = fit))
}


# ------------------------------------------------------
# OPTIMIZATION FRAMEWORK
# ------------------------------------------------------

# Perform k-fold cross-validation and compute AUC
evaluate_model_performance <- function(model_df, k = 5) {
  library(pROC)    # For computing ROC curves and AUC
  library(caret)   # For creating stratified folds
  
  folds <- createFolds(model_df$Label, k = k, list = TRUE)  # Create k stratified folds based on the 'Label' column
  aucs <- c()  # Initialize an empty vector to store AUCs for each fold
  
  for (i in seq_along(folds)) {  # Loop over each fold
    test_idx <- folds[[i]]  # Get indices for the test set
    train_df <- model_df[-test_idx, ]  # Training data (excluding the test set)
    test_df <- model_df[test_idx, ]    # Testing data
    
    # Skip fold if test set has only one class (ROC requires both positive and negative cases)
    if (length(unique(test_df$Label)) < 2) {
      warning(paste("Skipping fold", i, "due to single-class test set"))
      next
    }
    
    fit <- glm(Label ~ Influence, data = train_df, family = "binomial")  # Fit logistic regression using Influence to predict Label
    preds <- predict(fit, newdata = test_df, type = "response")  # Predict probabilities on the test set
    
    roc_obj <- roc(test_df$Label, preds)  # Compute ROC curve
    aucs <- c(aucs, auc(roc_obj))  # Store AUC for the current fold
  }
  
  # If no valid folds were found (e.g., every test set had only one class), stop with an error
  if (length(aucs) == 0) {
    stop("No valid folds with both classes. Try increasing sample size or stratifying folds.")
  }
  
  return(mean(aucs))  # Return mean AUC across all valid folds
}
