
ScoreProjection <- function(
  exprs, result, exprs_distance_matrix, projection_collapse_n, projection_neighbourhood, knn, knn.k, knn.algorithm, knn.distance, verbose
) {
  
  if (is.list(exprs))
    exprs <- do.call(rbind, exprs)
  if (is.list(result$Projection))
    result$Projection <- do.call(rbind, result$Projection)
  
  collapsed <- nrow(exprs) > projection_collapse_n
  
  if (collapsed) {
    if (verbose)
      .msg('\t-> computing k-NN graph for dimension-reduction scoring... ')
    systime <- NA
    knn_res <- ComputekNNMatrix(result$Projection, k = knn.k, method = knn.algorithm, dist = knn.distance, out.systime = systime)
    if (verbose)
      .msg_alt_good(' done in ', round(systime['elapsed'], 2), ' seconds\n')
    
    # Collapsed Co-ranking Matrix (based on k-NNGs)
    if (verbose)
      .msg('\t-> computing collapsed co-ranking matrix for dimension-reduction scoring... ')
    systime <- system.time({
      cQ <- CollapsedCorankingMatrix(knn$Indices, knn_res$Indices)
    })
    if (verbose)
      .msg_alt_good(' done in ', round(systime['elapsed'], 2), ' seconds\n')
    
    # Local Continuity Meta-Criterion (MCMC)
    N <- nrow(knn_res$Indices)
    cN <- nrow(cQ)
    K <- cN - 1
    lcmc <- K / (1 - N) + (1 / (N * K)) * sum(cQ[1:K, 1:K])
    
    # B_{NX} (from Lee and Verleysen, 2008) ... "Intrusiveness"
    Un <- 1 / (K * N) * sum(cQ[1:K, 1:K][upper.tri(cQ[1:K, 1:K])])
    Ux <- 1 / (K * N) * sum(cQ[1:K, 1:K][lower.tri(cQ[1:K, 1:K])])
    Bnx <- Un - Ux # Bnx > 0 <=> intrusive, Bnx < 0 <=> extrusive
    
    eT <- NA
    eC <- NA
    Q <- cQ
    rm(cQ)
  } else {
    
    n <- nrow(exprs)
    k <- projection_neighbourhood
    
    if (verbose)
      .msg('\t-> computing projection distance matrix... ')
    systime <- system.time({
      proj_distance_matrix <- coRanking:::euclidean_C(result$Projection)
    })
    if (verbose) {
      .msg_alt_good('done in ', round(systime['elapsed'], 2), ' seconds\n')
    }
    
    exprs_rankmatrix <- coRanking::rankmatrix(exprs_distance_matrix, input = 'dist')
    proj_rankmatrix <- coRanking::rankmatrix(proj_distance_matrix, input = 'dist')
    
    Q <- coRanking:::coranking_C(exprs_rankmatrix, proj_rankmatrix)
    lcmc <- coRanking::LCMC(Q, K = as.integer(k))
    
    Gk <- ifelse(k < n / 2, n*k*(2*n-3*k-1), n*(n-k)*(n-k-1))
    
    LL <- Q[(k+1):nrow(Q), 1:k] # ext (lower-left quadrant)
    UR <- Q[1:k, (k+1):ncol(Q)] # int (upper-right quadrant)
    
    eT <- 1-2/Gk*sum(LL * 1:nrow(LL))
    eC <- 1-2/Gk*sum(t(t(UR) * 1:ncol(UR)))
    
    Un <- 1 / (k * n) * sum(Q[1:k, 1:k][upper.tri(Q[1:k, 1:k])])
    Ux <- 1 / (k * n) * sum(Q[1:k, 1:k][lower.tri(Q[1:k, 1:k])])
    Bnx <- Un - Ux # Bnx > 0 <=> intrusive, Bnx < 0 <=> extrusive
    
    knn_res <- list(Indices = NA, Distances = NA)
    Q <- as.matrix(Q)
    class(Q) <- 'matrix'
  }
  
  list(
    'Layout k-NNG' = knn_res,
    'Collapsed' = collapsed,
    'Co-Ranking Matrix' = Q,
    'Local Continuity Meta-Criterion' = lcmc,
    'Relative Intrusiveness' = Bnx,
    'Trustworthiness' = eT,
    'Continuity' = eC,
    'Projection Neighbourhood' = if (collapsed) knn.k else projection_neighbourhood
  )
}
