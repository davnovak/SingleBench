
DeployClustering_Repeat <- function(
  input, subpipelines, fTrain, fExtract, fMap, seed, idcs_training, knn, exprs, n_iter, n_param, n_cores, parallelise, h5_path = NULL, idx.subpipeline = NULL, idx.n_param = NULL, out.intermediates = NULL
) {
  
  packages <- c(
    unlist(purrr::map(subpipelines, function(x) x$clustering$r_packages)),
    if (any(purrr::map_lgl(subpipelines, function(x) x$clustering$uses_python))) c('reticulate') else c()
  )
  
  if (parallelise) {
    
    if (is.null(n_cores))
      n_cores <- parallel::detectCores()
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    result <- foreach::foreach(
      idx_iter = seq_len(n_iter),
      .combine =
        function(l1, l2) list(ClusteringVector = cbind(l1$ClusteringVector, l2$ClusteringVector), Timing = c(l1$Timing, l2$Timing)),
      .packages = packages
    ) %dopar% {
      
      systime <- system.time({
        set.seed(seed + idx_iter - 1)
        res <- fTrain(
          input   = if (is.null(idcs_training)) input else input[idcs_training],
          n_param = n_param,
          knn     = knn,
          exprs   = exprs,
          save_intermediates = TRUE,
          h5_path = h5_path,
          idx.subpipeline    = idx.subpipeline,
          idx.n_param        = idx.n_param,
          out.intermediates  = if (!is.null(out.intermediates)) intermediates else NULL
        )
        if (!is.null(out.intermediates))
          eval.parent(substitute(out.intermediates <- intermediates))
        
        res <-
          if (is.null(idcs_training))
            fExtract(res)
        else
          fMap(res, input)
      })
      return(
        list(
          ClusteringVector = res,
          Timing = systime['elapsed']
        )
      )
    }
    parallel::stopCluster(cl)
    
    ClusteringVector <- purrr::map(seq_len(n_iter), function(idx_iter) result$ClusteringVector[, idx_iter])
    Timing           <- unlist(purrr::map(seq_len(n_iter), function(idx_iter) result$Timing[idx_iter]))
    
  } else { # if (!parallelise)
    
    ClusteringVector <- vector(mode = 'list', length = n_iter)
    Timing           <- vector(mode = 'list', length = n_iter)
    for (idx_iter in seq_len(n_iter)) {
      systime <- system.time({
        set.seed(seed + idx_iter - 1)
        res <- fTrain(
          input   = if (is.null(idcs_training)) input else input[idcs_training],
          n_param = n_param,
          knn     = knn,
          exprs   = exprs,
          save_intermediates = TRUE,
          h5_path = h5_path,
          idx.subpipeline    = idx.subpipeline,
          idx.n_param        = idx.n_param,
          out.intermediates  = if (!is.null(out.intermediates)) intermediates else NULL
        )
        if (!is.null(out.intermediates))
          eval.parent(substitute(out.intermediates <- intermediates))
        
        res <-
          if (is.null(idcs_training))
            fExtract(res)
        else
          fMap(res, input)
      })
      ClusteringVector[[idx_iter]] <- res
      Timing[[idx_iter]] <- systime['elapsed']
    }
    Timing <- unlist(Timing)
  }
  
  list(
    ClusteringVector = ClusteringVector,
    Timing = Timing
  )
}
