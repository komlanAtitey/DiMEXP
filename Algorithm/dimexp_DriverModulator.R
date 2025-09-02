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

# ------------------------------------------------------
# MAIN PIPELINE EXECUTION
# ------------------------------------------------------
load("data_GSE21422_Research_7.rdata")
load("DiMEXP_data.rdata")
samples_dataset <- data_GSE21422_Research_7
genes <- rownames(samples_dataset) 
rownames(DiMEXP_data) <- genes
driver_genes <- c("LTBP4", "BMP4","GREM1") 
non_driver_genes <- setdiff(rownames(DiMEXP_data), driver_genes)

# Step 1-3: CIS Computation and Top Influencers
CIS_matrix <- compute_CIS_matrix_parallel(
  expr_mat = DiMEXP_data,
  non_drivers = non_driver_genes,
  drivers = driver_genes,
  max_lag = 2,
  ncores = 4
)
top_influencers <- get_top_influencers_per_driver(CIS_matrix, top_n = 5)

# Step 4: Likelihood estimation
likelihood_output <- compute_transformation_likelihood(
  CIS_matrix, 
  driver_genes, 
  rownames(DiMEXP_data)
)
likelihood_df <- likelihood_output$model_df

# Step 5: Optimization evaluation
auc_mean <- evaluate_model_performance(likelihood_df)
cat("Mean AUC from cross-validation:", auc_mean, "\n")

# ------------------------------------------------------
# OUTPUT & VISUALIZATION
# ------------------------------------------------------

#-------------------- Network plot of top influencers --------------------

# A. Prepare edge list
all_edges <- do.call(rbind, top_influencers)
edge_list <- dplyr::select(
  all_edges,
  from = Non_Driver_Gene,
  to = Driver_Gene,
  weight = Influence_Score
)

# B. Create igraph object
g <- graph_from_data_frame(edge_list, directed = TRUE)

# C. Identify node types
non_driver_nodes <- unique(edge_list$from)
driver_nodes <- unique(edge_list$to)

# D. Assign colors
non_driver_color <- "#1f77b4"  # blue for all non-drivers
driver_palette <- hue_pal()(length(driver_nodes))
driver_colors <- setNames(driver_palette, driver_nodes)

# E. Assign group and color to each node
node_df <- data.frame(
  id = V(g)$name,
  label = V(g)$name,
  color = sapply(V(g)$name, function(node) {
    if (node %in% driver_nodes) driver_colors[[node]] else non_driver_color
  }),
  group = sapply(V(g)$name, function(node) {
    if (node %in% driver_nodes) node else "Non-Driver Genes"
  }),
  font.size = 20,#20 
  stringsAsFactors = FALSE
)

# F. Prepare edge data with width scaled
edge_df <- edge_list %>%
  dplyr::mutate(width = scales::rescale(weight, to = c(1, 10))) %>%  # scale edge width
  dplyr::rename(from = from, to = to) %>%
  dplyr::select(from, to, width)

# G. Build custom node legend: non-drivers and drivers
node_legend <- data.frame(
  label = c("Non-Driver Genes", driver_nodes),
  color = c(non_driver_color, driver_colors),
  shape = "dot",
  font.size = 20, #20
  stringsAsFactors = FALSE
)

# H. Build custom edge legend for influence
edge_legend <- data.frame(
  #label = c("Low influence →", "Highe →"),
  label = c(" →", " →"),
  arrows = "to",
  width = c(1, 10),
  color = "black",
  stringsAsFactors = FALSE
)

# I. Plot the interactive network with labeled legend
visNetwork(node_df, edge_df) %>%
  visEdges(arrows = "to", scaling = list(min = 1, max = 10)) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend(
    addNodes = node_legend,
    addEdges = edge_legend,
    useGroups = FALSE,  # We manually specify nodes
    position = "right"
  ) %>%
  visLayout(randomSeed = 123)

#-------------------- Knockout in-silico execution --------------------

# ----------------------------
# Define function to compute Granger causality
# ----------------------------
compute_gc <- function(geneA, geneB, expr_data, max_lag = 2) {
  # Create a data frame with time-series for geneA and geneB
  df <- data.frame(t(expr_data[c(geneA, geneB), ]))
  colnames(df) <- c("A", "B")
  
  # Fit VAR model (A and B are time series)
  model <- tryCatch(
    VAR(df, p = max_lag, type = "const"),
    error = function(e) return(NULL)
  )
  if (is.null(model)) return(c(F = NA, p = NA))
  
  # Perform Granger causality test: Does A Granger-cause B?
  test <- tryCatch(
    causality(model, cause = "A"),
    error = function(e) return(NULL)
  )
  if (is.null(test) || is.null(test$Granger)) return(c(F = NA, p = NA))
  
  # Extract test statistics
  stat <- tryCatch(test$Granger$statistic, error = function(e) NA)
  pval <- tryCatch(test$Granger$p.value, error = function(e) NA)
  
  return(c(F = stat, p = pval))
}

#driver_genes <- c("LTBP4", "BMP4","GREM1") 
# ----------------------------
#  Define non-driver and driver genes
# ----------------------------
non_drivers <- top_influencers$LTBP4[,2]
drivers <- c("LTBP4")

# ----------------------------
# Run perturbation consistency analysis
# ----------------------------

# Trim gene names
non_drivers <- trimws(non_drivers)
drivers <- trimws(drivers)

# Check existence
stopifnot(all(non_drivers %in% rownames(DiMEXP_data)))
stopifnot(all(drivers %in% rownames(DiMEXP_data)))

results <- data.frame()

for (nd in non_drivers) {
  for (d in drivers) {
    # Compute Granger causality with original data
    stat_orig <- compute_gc(nd, d, DiMEXP_data)
    
    # Knock out the non-driver by replacing its expression with near-zero noise
    expr_knockout <- DiMEXP_data
    expr_knockout[nd, ] <- rnorm(ncol(expr_knockout), mean = 0, sd = 1e-4)
    
    # Re-compute Granger causality after knockout
    stat_knock <- compute_gc(nd, d, expr_knockout)
    
    # Collect results
    results <- rbind(results, data.frame(
      non_driver = nd,
      driver = d,
      F_orig = stat_orig[1],
      p_orig = stat_orig[2],
      F_knock = stat_knock[1],
      p_knock = stat_knock[2]
    ))
  }
}

# ----------------------------
# Post-process for visualization
# ----------------------------
# Convert p-values to -log10(p-values) for easier visualization
results$logp_orig <- -log10(results$p_orig)
results$logp_knock <- -log10(results$p_knock)

# ----------------------------
# Visualize results with ggplot2
# ----------------------------
# Reshape the data for plotting
df_plot <- reshape2::melt(results, id.vars = c("non_driver", "driver"),
                          measure.vars = c("logp_orig", "logp_knock"),
                          variable.name = "Condition",
                          value.name = "-log10(p)")

# Rename conditions for plot clarity
df_plot$Condition <- factor(df_plot$Condition, levels = c("logp_orig", "logp_knock"),
                            labels = c("Original", "Knockout"))

# Plot bar chart comparing -log10(p-value) before and after knockout
ggplot(df_plot, aes(x = interaction(non_driver, driver), y = `-log10(p)`, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(
      "Original" = "blue",
      "Knockout" = "red"
      # Add more if needed
    )
  ) +
  labs(
    x = "Non-driver → Driver Pair",
    y = "-log10(p-value)",
    title = "Effect of Non-driver Knockout on Granger Causality"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  )

