
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


# =============================
# DiMEXP sample propagation pipeline
# =============================

# ----------------------------
# Build robust sample-level transition matrix P
# ----------------------------
build_sample_markov <- function(X, k = 15, sigma = NULL, local_scaling = TRUE,
                                alpha = 0.995, use_pca = TRUE, npc = 30) {
  if (!is.matrix(X)) stop("X must be a numeric matrix genes x samples")
  N <- ncol(X)
  
  # Samples x genes
  S <- t(X)
  
  # PCA reduction (safe)
  if (use_pca) {
    Scentered <- scale(S, center = TRUE, scale = FALSE)
    svd_res <- svd(Scentered)
    max_pc <- min(nrow(Scentered), ncol(Scentered))
    npc_use <- min(npc, max_pc)
    
    if (npc > max_pc) {
      message("⚠️ npc reduced to ", max_pc, " (max possible PCs).")
    }
    
    Sred <- Scentered %*% svd_res$v[, 1:npc_use, drop = FALSE]
  } else {
    Sred <- S
  }
  
  # Squared Euclidean distances
  d2 <- as.matrix(dist(Sred))^2
  
  if (is.null(sigma)) {
    nonzero <- d2[lower.tri(d2)]
    sigma0 <- sqrt(median(nonzero[nonzero > 0], na.rm = TRUE))
    if (is.na(sigma0) || sigma0 == 0) sigma0 <- 1e-6
    sigma <- sigma0
  }
  
  # Affinity matrix
  if (local_scaling) {
    knn_dists <- apply(d2, 1, function(row) {
      sort_s <- sort(row, decreasing = FALSE)
      kth <- sort_s[min(length(sort_s), k + 1)]
      sqrt(kth)
    })
    knn_dists[knn_dists == 0] <- min(knn_dists[knn_dists > 0])
    W <- matrix(0, nrow = N, ncol = N)
    for (i in 1:N) {
      for (j in 1:N) {
        W[i, j] <- exp(-d2[i, j] / (knn_dists[i] * knn_dists[j] + 1e-12))
      }
    }
  } else {
    W <- exp(-d2 / (2 * sigma^2))
  }
  
  # Keep k-NN
  for (i in 1:N) {
    ord <- order(W[i, ], decreasing = TRUE)
    keep <- ord[1:min(length(ord), k + 1)]
    drop <- setdiff(1:N, keep)
    W[i, drop] <- 0
  }
  
  # Symmetrize
  W <- pmax(W, t(W))
  diag(W) <- 0
  
  # Zero-row handling
  d <- rowSums(W)
  zero_rows <- which(d == 0)
  if (length(zero_rows) > 0) {
    W[zero_rows, ] <- 1 / N
    d[zero_rows] <- 1
  }
  
  # Row-stochastic transition matrix
  P <- W / rowSums(W)
  
  # Teleportation
  if (!is.null(alpha) && alpha < 1) {
    P <- alpha * P + (1 - alpha) / N
  }
  
  # Final renormalization
  P <- P / rowSums(P)
  
  return(list(P = P, W = W, sigma = sigma, degrees = d))
}

# ----------------------------
# Safe continuous-time generator (first-order approximation)
# ----------------------------
estimate_Q_from_P <- function(P, delta_t=1) {
  Q <- (P - diag(nrow(P))) / delta_t
  diag(Q) <- -rowSums(Q) + diag(Q)
  return(Q)
}

# ----------------------------
# Stationary distribution
# ----------------------------
compute_stationary <- function(P) {
  es <- eigen(t(P))
  idx <- which.max(Re(es$values))
  v <- Re(es$vectors[, idx])
  pi <- abs(v / sum(v))
  pi
}

# ----------------------------
# Leading eigenvectors (symmetric)
# ----------------------------
compute_leading_eig <- function(P, k=10) {
  Psym <- (P + t(P))/2
  es <- eigen(Psym)
  list(values=es$values[1:k], vectors=es$vectors[,1:k])
}

# ----------------------------
# MFPT (discrete)
# ----------------------------
mfpt_discrete <- function(P, B_idx) {
  N <- nrow(P); C_idx <- setdiff(seq_len(N), B_idx)
  if(length(C_idx)==0) return(rep(0,N))
  Pcc <- P[C_idx,C_idx,drop=FALSE]; Icc <- diag(length(C_idx))
  Nmat <- tryCatch(solve(Icc-Pcc), error=function(e) MASS::ginv(Icc-Pcc))
  t <- numeric(N); t[B_idx] <- 0; t[C_idx] <- rowSums(Nmat); t
}

# ----------------------------
# Committor (discrete)
# ----------------------------
committor_discrete <- function(P, A_idx, B_idx) {
  N <- nrow(P); T_idx <- setdiff(seq_len(N), c(A_idx,B_idx))
  q <- numeric(N); q[B_idx] <- 1; q[A_idx] <- 0
  if(length(T_idx)==0) return(q)
  Ptt <- P[T_idx,T_idx,drop=FALSE]; Ptb <- P[T_idx,B_idx,drop=FALSE]
  rhs <- rowSums(Ptb)
  q[T_idx] <- as.numeric(tryCatch(solve(diag(length(T_idx))-Ptt, rhs),
                                  error=function(e) MASS::ginv(diag(length(T_idx))-Ptt) %*% rhs))
  q
}

# ----------------------------
# Reactive flux
# ----------------------------
reactive_flux_discrete <- function(P, pi, q) {
  N <- nrow(P); kron <- matrix(0,N,N)
  for(i in 1:N) for(j in 1:N) kron[i,j] <- pi[i]*P[i,j]*q[i]*(1-q[j])
  kron
}

# ----------------------------
# Spectral coarse-graining / clustering
# ----------------------------
coarse_grain_spectral <- function(P, n_clusters=4, n_eigs=NULL) {
  N <- nrow(P)
  if(is.null(n_eigs)) n_eigs <- min(n_clusters+3,N-1)
  eig <- compute_leading_eig(P, k=n_eigs)
  evecs <- Re(eig$vectors)
  embed <- evecs[,1:n_clusters,drop=FALSE]
  emb_norm <- embed / sqrt(rowSums(embed^2)+1e-12)
  km <- kmeans(emb_norm, centers=n_clusters, nstart=20)
  list(membership=km$cluster, centers=km$centers, kmeans=km, embedding=emb_norm)
}

# ----------------------------
# Propagate gene-expression
# ----------------------------

# ---- Discrete-time stochastic diffusion with noticeable stochasticity --------
propagate_expression_discrete <- function(
    X, P, t=1, noise_level_rel=0.05, noise_level_abs=0.05, clamp_negatives=TRUE
) {
  Xprop <- X
  n_genes <- nrow(X)
  n_samples <- ncol(X)
  
  for(tt in 1:t) {
    # deterministic propagation
    Xdet <- t(P %*% t(Xprop))
    
    # stochastic perturbation: relative + absolute
    noise <- matrix(0, nrow=n_genes, ncol=n_samples)
    for(g in 1:n_genes) {
      noise[g, ] <- rnorm(n_samples, mean=0,
                          sd = noise_level_abs + noise_level_rel * sd(Xdet[g, ]))
    }
    
    Xprop <- Xdet + noise
    if(clamp_negatives) Xprop[Xprop < 0] <- 0
  }
  Xprop
}

# ----- Continuous-time stochastic diffusion with noticeable stochasticity -----
propagate_expression_continuous <- function(
    X, Q, t=1, noise_level_rel=0.05, noise_level_abs=0.05, clamp_negatives=TRUE
) {
  if(!requireNamespace("expm", quietly=TRUE)) stop('Install "expm" package')
  
  Xdet <- t(expm::expm(t*Q) %*% t(X))
  
  n_genes <- nrow(Xdet)
  n_samples <- ncol(Xdet)
  noise <- matrix(0, nrow=n_genes, ncol=n_samples)
  
  for(g in 1:n_genes) {
    noise[g, ] <- rnorm(n_samples, mean=0,
                        sd = noise_level_abs + noise_level_rel * sd(Xdet[g, ]))
  }
  
  Xprop <- Xdet + noise
  if(clamp_negatives) Xprop[Xprop < 0] <- 0
  Xprop
}

# ----------------------------
# Full upgraded wrapper
# ----------------------------
sample_markov_pipeline <- function(X, k=15, alpha=0.995, 
                                   t_discrete=5, t_continuous=5, 
                                   n_clusters=4, delta_t=1) {
  res <- list()
  
  # 1) Build transition matrix
  build <- build_sample_markov(X, k=k, alpha=alpha)
  P <- build$P; res$P <- P; res$W <- build$W; res$sigma <- build$sigma
  
  # 2) Continuous-time generator Q
  Q <- estimate_Q_from_P(P, delta_t=delta_t); res$Q <- Q
  
  # 3) Stationary distribution
  res$pi <- compute_stationary(P)
  
  # 4) Spectral embedding
  res$eig <- compute_leading_eig(P, k=min(10, ncol(P)))
  
  # 5) Coarse-graining (clusters)
  res$clusters <- coarse_grain_spectral(P, n_clusters=n_clusters)
  
  # 6) Mean first-passage time (to cluster 1)
  B_idx <- which(res$clusters$membership==1)
  res$mfpt_discrete <- mfpt_discrete(P, B_idx)
  
  # 7) Committor (cluster 2 vs cluster 1)
  A_idx <- which(res$clusters$membership==2)
  res$committor_discrete <- committor_discrete(P, A_idx, B_idx)
  
  # 8) Reactive flux
  res$reactive_flux <- reactive_flux_discrete(P, res$pi, res$committor_discrete)
  
  # 9) Gene-expression diffusion
  res$X_discrete <- propagate_expression_discrete(X, P, t=t_discrete)
  res$X_continuous <- propagate_expression_continuous(X, Q, t=t_continuous)
  
  return(res)
}


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







