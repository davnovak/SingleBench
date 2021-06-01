
HDF5_CreateHDF5AndWriteInputs <- function(
  benchmark,
  verbose
) {
  if (file.exists(benchmark$h5_path))
    file.remove(benchmark$h5_path)
  
  rhdf5::h5createFile(file = benchmark$h5_path)
  
  rhdf5::h5createGroup(file = benchmark$h5_path, group = 'Input')
  suppressMessages(.h5write(obj = benchmark$exprs, file = benchmark$h5_path, name = 'Input/ExpressionMatrix'))
  rhdf5::h5createGroup(file = benchmark$h5_path, group = 'Input/Layout')
  .h5write(obj = -1, file = benchmark$h5_path, name = 'Input/Layout/IsReferenceToSubpipeline')
  .h5write(obj = -1, file = benchmark$h5_path, name = 'Input/Layout/IsReferenceToNParamIteration')
  .h5write(obj = benchmark$column_names, file = benchmark$h5_path, name = 'Input/ColumnNames')
  .h5writeFactorVectorOrListOfThem(obj = benchmark$annotation, file = benchmark$h5_path, name = 'Input/Annotation')
  .h5write(obj = benchmark$row_count, file = benchmark$h5_path, name = 'Input/RowCount')
  .h5write(obj = benchmark$column_count, file = benchmark$h5_path, name = 'Input/ColumnCount')
  .h5write(obj = benchmark$bootstrap_indices, file = benchmark$h5_path, name = 'Input/BootstrapIndices')
}

HDF5_InitialiseEvaluationResults <- function(
  benchmark
) {
  slotname <- .h5_slotname() # 'EvaluationResults'
  
  ## If evaluation results from some previous unsuccesstul run are present, this will delete them first and create a new 'EvaluationResults' group
  g <- as.list(rhdf5::h5ls(benchmark$h5_path))$group
  g <- g[-grep('^/Input', g)]
  if (slotname %in% g)
    rhdf5::h5delete(benchmark$h5_path, slotname)
  rhdf5::h5createGroup(file = benchmark$h5_path, group = slotname)
  
  ## Save stability setting and random seeds
  .h5write(obj = benchmark$stability, file = benchmark$h5_path, name = .h5_slotname(suffix = 'Stability'))
  .h5write(obj = benchmark$seed.dimred, file = benchmark$h5_path, name = .h5_slotname(suffix = 'RandomSeed_Projection'))
  .h5write(obj = benchmark$seed.cluster, file = benchmark$h5_path, name = .h5_slotname(suffix = 'RandomSeed_Clustering'))
  
  ## Iterate over subpipelines
  purrr::walk(
    seq_len(benchmark$n_subpipelines),
    function(idx.subpipeline) {
      ## Create a slot for the subpipeline and a 'Projection' and 'Clustering' slot therein
      rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline))
      .h5write(obj = GetSubpipelineName(benchmark, idx.subpipeline), file = benchmark$h5_path, name = 'Name') # note down name of subpipeline
      rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection'))
      rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Clustering'))
      
      ## Extract the projection WrapperWithParameters object
      proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
      if (!is.null(proj)) {
        
        ## Iterate over the n-parameter values if specified
        no_npar <- FALSE
        n_param_range <- seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline))
        if (length(n_param_range) == 0) {
          n_param_range <- 1
          no_npar <- TRUE
        }
        
        purrr::walk(
          n_param_range,
          function(idx.n_param) {
            
            if (no_npar)
              idx.n_param <- NULL
            
            ## Create a 'Projection' group for every n-parameter iteration
            if (!is.null(idx.n_param))
              rhdf5::h5createGroup(
                file = benchmark$h5_path,
                group = .h5_slotname(
                  idx.subpipeline = idx.subpipeline,
                  tool_type = 'Projection',
                  idx.n_param = idx.n_param
                )
              )
            
            ## If multiple modules are chained within the projection step, create a group to store results for each module
            n_modules_proj <- GetProjectionModuleCount(benchmark, idx.subpipeline)
            if (!is.null(n_modules_proj)) {
              purrr::walk(
                seq_len(n_modules_proj - 1),
                function(idx.module)
                  rhdf5::h5createGroup(
                    file = benchmark$h5_path,
                    group = .h5_slotname(
                      idx.subpipeline = idx.subpipeline,
                      tool_type = 'Projection',
                      idx.n_param = idx.n_param,
                      idx.module = idx.module
                    )
                  )
              )
              rhdf5::h5createGroup(
                file = benchmark$h5_path,
                group = .h5_slotname(
                  idx.subpipeline = idx.subpipeline,
                  tool_type = 'Projection',
                  idx.n_param = idx.n_param,
                  suffix = 'Scores'
                )
              )
              if (!no_npar) {
                ## By default, reference index is set to -1 (~ no reference)
                .h5write(
                  obj = -1,
                  file = benchmark$h5_path,
                  name = .h5_slotname(
                    idx.subpipeline = idx.subpipeline,
                    tool_type = 'Projection',
                    idx.n_param = idx.n_param,
                    suffix = 'IsReferenceToSubpipeline'
                  )
                )
                .h5write(
                  obj = -1,
                  file = benchmark$h5_path,
                  name = .h5_slotname(
                    idx.subpipeline = idx.subpipeline,
                    tool_type = 'Projection',
                    idx.n_param = idx.n_param,
                    suffix = 'IsReferenceToNParamIteration'
                  )
                )
              }
            }
          })
        .h5write(
          obj = if (IsClone(proj)) proj$ref else -1,
          file = benchmark$h5_path,
          name = .h5_slotname(
            idx.subpipeline = idx.subpipeline,
            tool_type = 'Projection',
            suffix = 'IsReferenceToSubpipeline'
          )
        )
      }
      
      ## Extract the clustering WrapperWithParameters object
      clus <- benchmark$subpipelines[[idx.subpipeline]]$clustering
      if (!is.null(clus)) {
        
        ## Iterate over the n-parameter values if specified
        no_npar <- FALSE
        n_param_range <- seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline))
        if (length(n_param_range) == 0) {
          n_param_range <- 1
          no_npar <- TRUE
        }
        
        purrr::walk(
          n_param_range,
          function(idx.n_param) {
            
            if (no_npar)
              idx.n_param <- NULL
            
            ## Create a 'Clustering' group for every n-parameter iteration
            if (!is.null(idx.n_param))
              rhdf5::h5createGroup(
                file = benchmark$h5_path,
                group = .h5_slotname(
                  idx.subpipeline = idx.subpipeline,
                  tool_type = 'Clustering',
                  idx.n_param = idx.n_param
                )
              )
            
            ## Create slots to groups to store results, scores and (potential) reference pointers (for recycling results)
            # rhdf5::h5createGroup(
            #   file = benchmark$h5_path,
            #   group = .h5_slotname(
            #     idx.subpipeline = idx.subpipeline,
            #     tool_type = 'Clustering',
            #     idx.n_param = idx.n_param,
            #     suffix = 'ClusteringVector'
            #   )
            # )
            # rhdf5::h5createGroup(
            #   file = benchmark$h5_path,
            #   group = .h5_slotname(
            #     idx.subpipeline = idx.subpipeline,
            #     tool_type = 'Clustering',
            #     idx.n_param = idx.n_param,
            #     suffix = 'Timing'
            #   )
            # )
            rhdf5::h5createGroup(
              file = benchmark$h5_path,
              group = .h5_slotname(
                idx.subpipeline = idx.subpipeline,
                tool_type = 'Clustering',
                idx.n_param = idx.n_param,
                suffix = 'Scores'
              )
            )
          }
        )
      }
  }) # subpipeline iterations
  invisible(benchmark)
}
