
EvalProjection <- function(
  benchmark, verbose
) {
  
  seed <- benchmark$seed.projection
  
  purrr::walk(
    seq_len(benchmark$n_subpipelines), function(idx.subpipeline) {
      
      subpipeline_name  <- GetSubpipelineName(benchmark, idx.subpipeline)
      
      proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
      if (!is.null(proj)) {
        
        if (verbose) {
          .msg('Evaluating subpipeline '); .msg_alt(idx.subpipeline); .msg(' of ' ); .msg_alt(benchmark$n_subpipelines); .msg('\n')
        }
        
        if (!IsClone(proj)) {
          
          n_param_values    <- benchmark$n_params[[idx.subpipeline]]$projection
          train             <- fTrain.ModuleChain(proj)
          extract           <- fExtract.ModuleChain(proj)
          map               <- fMap.ModuleChain(proj)
          
          exprs <- GetExpressionMatrix(benchmark)
          knn   <- if (proj$uses_knn_graph) GetkNNMatrix(benchmark) else NULL
          
          nparam_range <- seq_along(n_param_values)
          if (length(nparam_range) == 0)
            nparam_range <- 'NoNParameter'
          
          purrr::walk(
            nparam_range,
            function(idx.n_param) {
              
              if (idx.n_param != 'NoNParameter')
                idx_ref <- GetProjectionReference(benchmark, idx.subpipeline, idx.n_param)
              
              if (idx.n_param == 'NoNParameter' || idx_ref <= 0) {
                
                if (idx.n_param != 'NoNParameter' && is.na(n_param_values[idx.n_param])) {
                    if (verbose)
                      .msg('\t-> omitting projection step (since n-parameter set to NA)\n')
                    .h5writeProjectionSubpipelineReference(benchmark, idx.subpipeline, idx.n_param, idx_ref = 0)
                    .h5writeProjectionScoring(NA, benchmark, idx.subpipeline, idx.n_param)
                } else {
                  if (verbose) {
                    .msg('\t-> evaluating projection step: ')
                    .msg_alt(GetNParameterIterationName_Projection(benchmark, idx.subpipeline, idx.n_param))
                  }
                  n_param <-
                    if (idx.n_param == 'NoNParameter')
                      NULL
                    else
                      n_param_values[idx.n_param]
                  if (idx.n_param == 'NoNParameter')
                    idx.n_param <- NULL
                  
                  this_exprs <-
                    if (proj$uses_original_expression_matrix)
                      exprs
                    else
                      NULL
                  
                  res <- DeployProjection(exprs, train, extract, map, seed, benchmark$projection.training_set, knn, this_exprs, n_param, benchmark$h5_path, idx.subpipeline = idx.subpipeline, idx.n_param = idx.n_param)
                  if (verbose)
                    .msg_alt_good(' done in ', round(res$Timing, 2), ' seconds\n')
                  res$Projection <- SeparateIntoSamples(res$Projection, benchmark)
                  .h5writeProjectionResult(res, benchmark, idx.subpipeline, idx.n_param)
                  if (benchmark$score_projections) {
                    distmat <- if (benchmark$row_count <= benchmark$projection_collapse_n) GetDistanceMatrix(benchmark) else NULL
                    knn <- if (benchmark$row_count > benchmark$projection_collapse_n) GetkNNMatrix(benchmark) else NULL
                    scores <- ScoreProjection(exprs, res, distmat, benchmark$projection_collapse_n, benchmark$projection_neighbourhood, knn, benchmark$knn.k, benchmark$knn.algorithm, benchmark$knn.distance, verbose)
                    .h5writeProjectionScoring(scores, benchmark, idx.subpipeline, idx.n_param)
                  }
                }
                
              } else {
                if (verbose)
                  .msg('\t-> cloning projection result of n-parameter iteration ', idx_ref, '\n')
              }
              
            }
          )
        } else {
          if (verbose)
            .msg('\t-> cloning projection result of subpipeline ', proj$ref, '\n')
        }
      }
    }
  )
  
  invisible(benchmark)
}

DeployProjection <- function(
  input, fTrain, fExtract, fMap, seed, idcs_training, knn, exprs, n_param, h5_path = NULL, idx.subpipeline = NULL, idx.n_param = NULL, out.intermediates = NULL
) {
  systime <- system.time({
    set.seed(seed)
    intermediates <- NA
    res <- 
      fTrain(
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
  
  colnames(res) <- paste0('component_', seq_len(ncol(res)))
  list(Projection = res, Timing = systime['elapsed'])
}

