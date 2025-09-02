
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

