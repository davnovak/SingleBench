
DeployClustering_SingleRun <- function(
  input, subpipelines, fTrain, fExtract, fMap, seed, idcs_training, knn, exprs, n_param, h5_path = NULL, idx.subpipeline = NULL, idx.n_param = NULL, out.intermediates = NULL
) {
  systime <- system.time({
    set.seed(seed)
    intermediates <- NA
    res <- fTrain(
      input              = if (is.null(idcs_training)) input else input[idcs_training],
      n_param            = n_param,
      knn                = knn,
      exprs              = exprs,
      save_intermediates = TRUE,
      h5_path            = h5_path,
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
  
  list(
    ClusteringVector = res, Timing = systime['elapsed']
  )
}
