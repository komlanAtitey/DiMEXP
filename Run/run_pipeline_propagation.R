
################################################################
################################################################ Run pipeline
################################################################
# =============================
# Run the integrated sample-level DiMEXP pipeline
# =============================
load("data_GSE21422_Research_7.rdata")
samples_dataset <- data_GSE21422_Research_7

####### preprocess the data
process_genomic_data <- function(df, npc_use = NULL) {
  # 1) Identify gene symbol column (character or factor)
  gene_col <- which(sapply(df, function(x) is.character(x) || is.factor(x)))
  if (length(gene_col) != 1) {
    stop("Data must contain exactly one gene symbol column (character or factor).")
  }
  
  # 2) Identify sample columns (numeric)
  sample_cols <- which(sapply(df, is.numeric))
  if (length(sample_cols) == 0) {
    stop("No numeric sample columns found.")
  }
  
  # 3) Reorder: gene symbol first, then samples
  df_ordered <- df[, c(gene_col, sample_cols), drop = FALSE]
  
  # 4) Convert to list for $ access
  df_list <- as.list(df_ordered)
  
  # 5) Optional: PCA/SVD
  if (!is.null(npc_use) && npc_use > 0) {
    mat <- as.matrix(df_ordered[, -1])   # exclude gene column
    mat <- scale(mat, center = TRUE, scale = TRUE)  # normalize
    
    svd_res <- svd(mat)
    max_pc <- min(nrow(mat), ncol(mat))
    
    # auto-adjust npc_use if too high
    npc_final <- min(npc_use, max_pc)
    if (npc_use > max_pc) {
      message("⚠️ npc_use reduced to ", npc_final, " (max possible PCs).")
    }
    
    df_list$svd_u <- svd_res$u[, 1:npc_final, drop = FALSE]
    df_list$svd_v <- svd_res$v[, 1:npc_final, drop = FALSE]
    df_list$svd_d <- svd_res$d[1:npc_final]
  }
  df_list <- data.frame(df_list)
  return(df_list)
}

if (is.character(samples_dataset[[1]]) || is.factor(samples_dataset[[1]])) {
  # Case A: first column is gene symbols
  data_processing <- process_genomic_data(samples_dataset)
  data_pipeline <- data_processing[, -1]
  data_pipeline <- as.matrix(data_pipeline)
  rownames(data_pipeline) <- NULL
  data_pipeline[is.na(data_pipeline)] <- 0
} else {
  # Case B: gene symbols are in rownames
  data_pipeline <- as.matrix(samples_dataset)
  data_pipeline[is.na(data_pipeline)] <- 0
}

####### run the pipeline
res <- sample_markov_pipeline(data_pipeline,
                              k=10,
                              alpha=0.99,
                              t_discrete=5,
                              t_continuous=5,
                              n_clusters=3) 

####### collect the discrete data
DiMEXP_data <- res$X_discrete
sample_names <- colnames(data_pipeline)
colnames(DiMEXP_data) <- sample_names

################################################################
################################################################ Comparing lines plots
################################################################
# ----------------------------
# Compare original vs propagated for first genes
# ----------------------------
data_orig_simul <- data_pipeline[,1:10]
data_DiMEXP_simul <- DiMEXP_data[,1:10]
gene_vector <- rownames(data_pipeline)

N <- 10 # sample size
data_plot <- data.frame(
  Sample = rep(1:N, 3),
  Expression = c(data_orig_simul[1,],
                 data_orig_simul[2,], 
                 data_orig_simul[3,]),
  Expression_Markov = c(data_DiMEXP_simul[1,],
                        data_DiMEXP_simul[2,], 
                        data_DiMEXP_simul[3,]),
  Gene = factor(rep(gene_vector[1:3], each=N))
)

# ----------------------------
# Plot original vs propagated expression
# ----------------------------
ggplot(data_plot, aes(x=Sample)) +
  geom_line(aes(y = Expression, color = "Original"), 
            size = 1.2) +
  geom_line(aes(y = Expression_Markov, color = "DiMEXP Propagated"), 
            size = 1.2) +
  facet_wrap(~Gene, scales="free_y") +
  labs(title="Sample-level DiMEXP Propagation of Gene Expression",
       y="Expression", color="Legend") +
  scale_color_manual(values = c("Original" = "blue",         # Blue
                                "DiMEXP Propagated" = "red2")) +  # Red
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  )

################################################################
################################################################ Heatmap
################################################################
# ----------------------------
# Prepare data for heatmap
# ----------------------------
# heatmap of the 10×10 gene-expression matrix to visualize the effect of sample-level Markov propagation. 
# show original vs propagated expression side by side.
library(reshape2)


# Original matrix
data_orig <- melt(data_orig_simul[1:6,])
colnames(data_orig) <- c("Gene", "Sample", "Expression")
data_orig$Type <- "Original"

# Propagated matrix
data_DiMEXP <- melt(data_DiMEXP_simul[1:6,])
colnames(data_DiMEXP) <- c("Gene", "Sample", "Expression")
data_DiMEXP$Type <- "Markov Propagated"

# Combine
data_combined <- rbind(data_orig, data_DiMEXP)
data_combined$Type <- factor(data_combined$Type, levels=c("Original","Markov Propagated"))

# ----------------------------
# Plot heatmap with both panels
# ----------------------------
ggplot(data_combined, aes(x=Sample, y=Gene, fill=Expression)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=mean(data_combined$Expression)) +
  facet_wrap(~Type, ncol=1) +
  labs(title="Gene Expression Heatmap: Original vs Markov Propagated",
       x="Sample", y="Gene", fill="Expression") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# This small example demonstrates how gene expressions “diffuse” across similar
# samples according to the Markov transition matrix.







