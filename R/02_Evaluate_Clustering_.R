
EvalClustering <- function(
  benchmark, verbose, no_parallelisation = FALSE, which_subpipelines = NULL
) {
  
  seed <- benchmark$seed.clustering
  
  seq_pipelines <- if (is.null(which_subpipelines)) seq_len(benchmark$n_subpipelines) else which_subpipelines
  
  ## Iterate over subpipelines
  purrr::walk(
    seq_pipelines, function(idx.subpipeline) {
      
      subpipeline_name  <- GetSubpipelineName(benchmark, idx.subpipeline)

      clus <- benchmark$subpipelines[[idx.subpipeline]]$clustering
      if (!is.null(clus)) {
        
        if (verbose) { .msg('Evaluating subpipeline '); .msg_alt(idx.subpipeline); .msg(' of ' ); .msg_alt(benchmark$n_subpipelines); .msg('\n') }
        
        ## Get n-parameter values for this subpipeline (either from the clustering step or the projection step...)
        npar_proj <- FALSE
        n_param_values <- benchmark$n_params[[idx.subpipeline]]$clustering
        if (length(n_param_values) == 0) {
          n_param_values <- benchmark$n_params[[idx.subpipeline]]$projection
          proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
          idx <- idx.subpipeline
          while (IsClone(proj)) {
            idx <- proj$ref
            proj <- benchmark$subpipelines[[idx]]$projection
          }
          n_param_values <- benchmark$n_params[[idx]]$projection
          if (length(n_param_values) > 0)
            npar_proj <- TRUE
        }
        
        ## Get model-building functions
        train <- fTrain.ModuleChain(clus)
        extract <- fExtract.ModuleChain(clus)
        map <- fMap.ModuleChain(clus)
        idcs_training <- benchmark$clustering.training_set
        bootstrap_idcs <- benchmark$bootstrap_indices
        
        n_iter <- benchmark$stability.n_iter
        n_cores <- benchmark$n_cores
        parallelise <- all(purrr::map_lgl(clus$modules, function(x) !x$wrapper_with_parameters$wrapper$prevent_parallel_execution))
        
        if (no_parallelisation)
          parallelise <- FALSE
        
        ## Get expression matrix and k-NNG (if needed)
        exprs <- GetExpressionMatrix(benchmark)
        knn <- if (clus$uses_knn_graph) GetkNNMatrix(benchmark) else NULL
        
        n_param_range <- seq_along(n_param_values)
        no_npar <- FALSE
        if (length(n_param_range) == 0) {
          no_npar <- TRUE
          n_param_range <- 1
        }
        
        ## Iterate over n-parameter values (if n-parameter is specified)
        purrr::walk(
          n_param_range,
          function(idx.n_param) {
            if (verbose) {
              .msg('\t-> evaluating clustering step: ')
              .msg_alt(GetNParameterIterationName_Clustering(benchmark, idx.subpipeline, idx.n_param))
            }
            
            n_param <- if (no_npar || npar_proj) NULL else n_param_values[idx.n_param]
            input <- GetClusteringInput(benchmark, idx.subpipeline, idx.n_param = if (no_npar) NULL else idx.n_param)
            this_exprs <- if (clus$uses_original_expression_matrix) exprs else NULL
            
            res <-
              if (benchmark$stability == 'single')
                DeployClustering_SingleRun(input, benchmark$subpipelines, train, extract, map, seed, idcs_training, knn, this_exprs, n_param, h5_path = benchmark$h5_path, idx.subpipeline = idx.subpipeline, idx.n_param = idx.n_param)
              else if (benchmark$stability == 'repeat')
                DeployClustering_Repeat(input, benchmark$subpipelines, train, extract, map, seed, idcs_training, knn, this_exprs, n_iter, n_param, n_cores, parallelise, h5_path = benchmark$h5_path, idx.subpipeline = idx.subpipeline, idx.n_param = idx.n_param)
              else if (benchmark$stability == 'bootstrap')
                DeployClustering_Bootstrap(input, benchmark$subpipelines, train, extract, map, seed, idcs_training, bootstrap_idcs, knn, this_exprs, n_param, n_iter, n_cores, parallelise, h5_path = benchmark$h5_path, idx.subpipeline = idx.subpipeline, idx.n_param = idx.n_param)
            
            if (verbose) {
              if (length(res$Timing) == 1) {
                .msg_alt_good(' done in ', round(res$Timing, 2), ' seconds\n')
              } else {
                .msg_alt_good(' ', length(res$Timing), ' runs done in an average of ', round(mean(res$Timing), 2), ' seconds each\n')
              }
            }
            
            res$ClusteringVector <- SeparateIntoSamples(res$ClusteringVector, benchmark)
            
            if (no_npar && !npar_proj)
              idx.n_param <- NULL
            .h5writeClusteringResult(res, benchmark, idx.subpipeline, idx.n_param)
            
            if (verbose)
              .msg('\t\t-> computing scores...')
            scores <- ScoreClustering(exprs, GetAnnotation(benchmark, concatenate = TRUE), res, benchmark$stability, bootstrap_idcs, benchmark$unassigned_labels, column_names = benchmark$column_names)
            if (verbose)
              .msg(' writing scores...')
            .h5writeClusteringScoring(scores, benchmark, idx.subpipeline, idx.n_param)
            if (verbose)
              .msg_alt_good(' done\n')
            
            ## Compute cluster medians
            
            exprs <- GetExpressionMatrix(benchmark, concatenate = TRUE)
            codes <- GetClustering(benchmark, idx.subpipeline, idx.n_param)
            codes <- as.factor(codes)
            
            meds <-
              matrix(
                apply(expand.grid(levels(codes), benchmark$column_names), 1, function(x) median(exprs[codes == x[1], x[2], drop = FALSE])),
                ncol = benchmark$column_count,
                dimnames = list(levels(codes), benchmark$column_names)
              )
            .h5writeClusterMedians(meds, benchmark, idx.subpipeline, idx.n_param)
              
          }
        )
      } # endif !is.null(clus)
    }
  )
  
  invisible(benchmark)
}
