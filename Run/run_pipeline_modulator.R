
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

