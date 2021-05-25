
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
  .h5write(obj = -1, file = benchmark$h5_path, name = 'Input/Layout/IsReferenceTo')
  .h5write(obj = benchmark$column_names, file = benchmark$h5_path, name = 'Input/ColumnNames')
  .h5writeFactorVectorOrListOfThem(obj = benchmark$annotation, file = benchmark$h5_path, name = 'Input/Annotation')
  .h5write(obj = benchmark$row_count, file = benchmark$h5_path, name = 'Input/RowCount')
  .h5write(obj = benchmark$column_count, file = benchmark$h5_path, name = 'Input/ColumnCount')
  .h5write(obj = benchmark$bootstrap_indices, file = benchmark$h5_path, name = 'Input/BootstrapIndices')
}

HDF5_InitialiseEvaluationResults <- function(
  benchmark
) {
  slotname <- .h5_slotname() # 'EvaluationResults
  
  g <- as.list(rhdf5::h5ls(benchmark$h5_path))$group
  g <- g[-grep('^/Input', g)]
  if (slotname %in% g)
    rhdf5::h5delete(benchmark$h5_path, slotname)
  
  rhdf5::h5createGroup(file = benchmark$h5_path, group = slotname)
  .h5write(obj = benchmark$stability, file = benchmark$h5_path, name = .h5_slotname(suffix = 'Stability'))
  .h5write(obj = benchmark$seed.dimred, file = benchmark$h5_path, name = .h5_slotname(suffix = 'RandomSeed_Projection'))
  .h5write(obj = benchmark$seed.cluster, file = benchmark$h5_path, name = .h5_slotname(suffix = 'RandomSeed_Clustering'))
  purrr::walk(
    seq_len(benchmark$n_subpipelines),
    function(idx.subpipeline) {
      rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline))
      .h5write(obj = GetSubpipelineName(benchmark, idx.subpipeline), file = benchmark$h5_path, name = 'Name')
      rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection'))
      rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Clustering'))
      
      proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
      if (!is.null(proj) && !IsClone(proj)) {
        purrr::walk(
          seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline)),
          function(idx.n_param) {
            rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection', idx.n_param = idx.n_param))
            n_modules_proj <- GetProjectionModuleCount(benchmark, idx.subpipeline)
            if (!is.null(n_modules_proj)) {
              purrr::walk(
                seq_len(n_modules_proj - 1),
                function(idx.module)
                  rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection', idx.n_param = idx.n_param, idx.module = idx.module))
              )
              rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection', idx.n_param = idx.n_param, suffix = 'Result'))
              rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection', idx.n_param = idx.n_param, suffix = 'Scores'))
              .h5write(obj = -1, file = benchmark$h5_path, name = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection', suffix = 'IsReferenceTo'))
            }
            ref_val <- -1
            npar_val <- benchmark$n_params[[idx.subpipeline]]$projection[idx.n_param]
            if (idx.n_param %in% which(duplicated(benchmark$n_params[[idx.subpipeline]]$projection))) {
              ref_val <- which(benchmark$n_params[[idx.subpipeline]]$projection == npar_val)[1]
            }
            
            .h5write(obj = ref_val, file = benchmark$h5_path, name = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection', idx.n_param = idx.n_param, suffix = 'IsReferenceTo'))
          })
      } else {
        .h5write(obj = proj$ref, file = benchmark$h5_path, name = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Projection', suffix = 'IsReferenceTo'))
      }
      
      nparam_range <- seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline))
      if (IsClone(proj)) {
        nparam_range <- seq_len(GetNParameterIterationsCount(benchmark, proj$ref))
      }
      
      purrr::walk(
        nparam_range,
        function(idx.n_param) {
          rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Clustering', idx.n_param = idx.n_param))
          n_modules_clus <- GetClusteringModuleCount(benchmark, idx.subpipeline)
          if (!is.null(n_modules_clus)) {
            purrr::walk(
              seq_len(n_modules_clus - 1),
              function(idx.module)
                rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Clustering', idx.n_param = idx.n_param, idx.module = idx.module))
            )
            rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Clustering', idx.n_param = idx.n_param, suffix = 'Result'))
            rhdf5::h5createGroup(file = benchmark$h5_path, group = .h5_slotname(idx.subpipeline = idx.subpipeline, tool_type = 'Clustering', idx.n_param = idx.n_param, suffix = 'Scores'))
          }
          
      })
  })
}
