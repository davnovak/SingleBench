
EvalProjection <- function(
  benchmark, verbose
) {
  
  seed <- benchmark$seed.projection
  
  ## Iterate over subpipelines
  purrr::walk(
    seq_len(benchmark$n_subpipelines), function(idx.subpipeline) {
      
      subpipeline_name  <- GetSubpipelineName(benchmark, idx.subpipeline)
      proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
      
      if (!is.null(proj)) {
        if (verbose) { .msg('Evaluating subpipeline '); .msg_alt(idx.subpipeline); .msg(' of ' ); .msg_alt(benchmark$n_subpipelines); .msg('\n') }
        
        cloned <- IsClone(proj)
        if (cloned)
          proj <- benchmark$subpipelines[[proj$ref]]$projection
        
        ## See if there is an n-parameter and get the parameter value range if there is
        n_param_values <- benchmark$n_params[[idx.subpipeline]]$projection
        n_param_range <- seq_along(n_param_values)
        if (length(n_param_range) == 0)
          n_param_range <- 'NoNParameter'
        
        if (IsClone(proj) && n_param_range == 'NoNParameter') {
          ## If this projection step was used before and there is no n-parameter, write a reference to the result
          if (verbose) { .msg('\t-> cloning projection result from sub-pipeline'); .msg_alt(proj$ref); .msg('\n') }
          .h5writeProjectionReference(benchmark, idx.subpipeline = idx.subpipeline, idx.subpipeline_ref = proj$ref)
          
        } else {
          
          if (IsClone(proj) && n_param_range != 'NoNParameter') {
            ## If this projection step was used before and there is an n-parameter, check if there are n-parameter values for which we already have the result
            proj_result_locations <-
              purrr::map(
                seq_along(n_param_values),
                function(idx)
                  FindProjectionResultIfAlreadyGenerated(
                    benchmark,
                    idx.subpipeline = idx.subpipeline,
                    idx.n_param = idx,
                    n_param = n_param_values[idx]
                  )
              )
            which_proj_results_not_available_already <-
              purrr::map_lgl(proj_result_locations, function(x) is.null(x))
            
            ## Write reference to those n-parameter iteration results that we already have
            for (idx_res in seq_along(proj_result_locations)) {
              loc <- proj_result_locations[[idx_res]]
              if (!is.null(loc)) {
                .h5writeProjectionReference(
                  benchmark,
                  idx.subpipeline = idx.subpipeline,
                  idx.n_param = idx_res,
                  idx.subpipeline_ref = loc$idx.subpipeline,
                  idx.n_param_ref = loc$idx.n_param
                )
              }
            }
            
            ## Only those n-parameter iterations that have not been evaluated yet will be evaluated in this parameter sweep
            n_param_range <- n_param_range[which_proj_results_not_available_already]
          }
          
          
          if (length(n_param_range) > 0) { # if there are any unevaluated n-parameter iterations left...
            
            ## Get model-building functions
            train <- fTrain.ModuleChain(proj)
            extract <- fExtract.ModuleChain(proj)
            map <- fMap.ModuleChain(proj)
            
            ## Retrieve input expression data (and k-NNG if needed)
            exprs <- GetExpressionMatrix(benchmark)
            knn <- if (proj$uses_knn_graph) GetkNNMatrix(benchmark) else NULL
            
            ## Iterate over n-parameter values
            purrr::walk(
              n_param_range,
              function(idx.n_param) {
                
                if (idx.n_param != 'NoNParameter' && is.na(n_param_values[idx.n_param])) {
                  ## Skip the projection step if n-parameter value is set to NA
                  if (verbose) { .msg('\t-> omitting projection step (since n-parameter value set to NA)\n') }
                  .h5writeProjectionReference(
                    benchmark,
                    idx.subpipeline = idx.subpipeline,
                    idx.n_param = idx.n_param,
                    idx.subpipeline_ref = 0,
                    idx.n_param_ref = NA
                  )
                } else {
                  
                  if (verbose) { .msg('\t-> evaluating projection step: '); .msg_alt(GetNParameterIterationName_Projection(benchmark, idx.subpipeline, idx.n_param)) }
                  
                  if (idx.n_param != 'NoNParameter') {
                    
                    ## Check if n-parameter iteration was evaluated already in *this* subpipeline
                    loc <- FindProjectionResultIfAlreadyGenerated(
                      benchmark,
                      idx.subpipeline = idx.subpipeline,
                      idx.n_param = idx.n_param,
                      n_param = n_param_values[idx.n_param]
                    )
                    
                    if (!is.null(loc)) {
                      ## Write reference to previous n-parameter iteration if it has
                      .h5writeProjectionReference(
                        benchmark,
                        idx.subpipeline = idx.subpipeline,
                        idx.n_param = idx.n_param,
                        idx.subpipeline_ref = loc$idx.subpipeline,
                        idx.n_param_ref = loc$idx.n_param
                      )
                      if (verbose) { .msg_alt_good(' cloned from subpipeline ', loc$idx.subpipeline, ' n-param iteration ', loc$idx.n_param, '\n') }
                    }
                  }
                  
                  if (
                    idx.n_param == 'NoNParameter' ||
                    (idx.n_param != 'NoNParameter' && is.null(loc))
                  ) { # if the result has not not been generated before in this nor another subpipeline...
                    
                    # ...get n-parameter value
                    n_param <- if (idx.n_param == 'NoNParameter') NULL else n_param_values[idx.n_param]
                    if (idx.n_param == 'NoNParameter') idx.n_param <- NULL
                    
                    # ...get copy of expression data if needed
                    this_exprs <- if (proj$uses_original_expression_matrix) exprs else NULL
                    
                    # ...deploy the projection tool
                    res <- DeployProjection(exprs, train, extract, map, seed, benchmark$projection.training_set, knn, this_exprs, n_param, benchmark$h5_path, idx.subpipeline = idx.subpipeline, idx.n_param = idx.n_param)
                    if (verbose) .msg_alt_good(' done in ', round(res$Timing, 2), ' seconds\n')
                    
                    # ...if there were separate samples on input, separate out the projection also
                    res$Projection <- SeparateIntoSamples(res$Projection, benchmark)
                    
                    # ...write the projection result
                    .h5writeProjectionResult(res, benchmark, idx.subpipeline, idx.n_param)
                    
                    # ...score the projection result if needed
                    if (benchmark$score_projections) {
                      distmat <- if (benchmark$row_count <= benchmark$projection_collapse_n) GetDistanceMatrix(benchmark) else NULL
                      knn <- if (benchmark$row_count > benchmark$projection_collapse_n) GetkNNMatrix(benchmark) else NULL
                      scores <- ScoreProjection(exprs, res, distmat, benchmark$projection_collapse_n, benchmark$projection_neighbourhood, knn, benchmark$knn.k, benchmark$knn.algorithm, benchmark$knn.distance, verbose)
                      .h5writeProjectionScoring(scores, benchmark, idx.subpipeline, idx.n_param)
                    }
                  }
                }
              }
            )
          } # endif length(n_param_range) > 0
        }
      } # endif !is.null(proj)
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
    
    res <- if (is.null(idcs_training)) fExtract(res) else fMap(res, input)
  })
  
  colnames(res) <- paste0('component_', seq_len(ncol(res)))
  list(Projection = res, Timing = systime['elapsed'])
}

