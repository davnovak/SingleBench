
.WhatIs.ValidityChecks <- function(e) {
  if (!is.null(e$cluster) && is.null(e$idx.subpipeline)) {
    stop('For retrieving a cluster, "idx.subpipeline" must be given')
  }
}

.Plot.ValidityChecks <- function(e) {
  if (class(e$benchmark) != 'Benchmark')
    stop('"benchmark" not a valid Benchmark object')
}

.PlotProjection.ValidityChecks <- function(e) {
  .Plot.ValidityChecks(e)
  if (!e$benchmark$score_projection)
    stop('In this benchmark, performance of projection steps was not scored in evaluation')
  if (is.null(e$benchmark$subpipelines[[e$idx.subpipeline]]$projection))
    stop(paste0('Subpipeline ', e$idx.subpipeline, ' does not contain a projection step'))
}

.PlotClustering.ValidityChecks <- function(e) {
  .Plot.ValidityChecks(e)
  if (is.null(e$benchmark$subpipelines[[e$idx.subpipeline]]$clustering))
    stop(paste0('Subpipeline ', e$idx.subpipeline, ' does not contain a clustering step'))
}

.WrapTool.ValidityChecks <- function(e) {
  if (!is.atomic(e$name) || !is.character(e$name))
    stop('Invalid "name"')
  if (!e$type[1] %in% c('projection', 'clustering'))
    stop('Invalid "type"')
}

.Fix.ValidityChecks <- function(e) {
  if (is.null(e$wrapper))
    stop(paste0('Neither "wrapper.projection.', e$tool_name, '" nor "wrapper.clustering.', e$tool_name, '" found in global namespace'))
  param_names <- names(e$params)
  for (param_name in param_names)
    if (!param_name %in% e$wrapper$args_to_train || param_name == 'input')
      stop(paste0('"', param_name, '" is not a valid parameter to fix for ', e$tool_name))
}

.Module.ValidityChecks <- function(e) {
  if (!is.null(e$n_param) && ((!e$n_param %in% e$wrapper_with_parameters$wrapper$args_to_train) || e$n_param == 'input'))
    stop(paste0('"', e$n_param, '" is not a valid n-parameter for ', e$tool_name))
}

.Chain.ValidityChecks <- function(e) {
  if (length(e$which_n_param) > 1)
    stop('When chaining multiple modules in a projection or clustering step, only one tool can have an adjustable n-parameter')
  if ('ProjectionModule' %in% e$type && 'ClusteringModule' %in% e$type)
    stop('You cannot chain projection and clustering modules in a single step')
}

.Subpipeline.ValidityChecks <- function(e) {
  
}

.AlignSubpipelines.ValidityChecks <- function(e) {
  
}

.Benchmark.ValidityChecks <- function(e) {
  
  if (length(e$input_class) > 1 || !e$input_class %in% c('fcs_paths', 'SummarizedExperiment', 'flowSet'))
    stop('Invalid input type')
  
  if (!TryWritePermission())
    stop('Write permission is needed in order to write the auxiliary HDF5 file')
  
  if (!is.character(e$h5_path) || length(e$h5_path) != 1 || nchar(e$h5_path) == 0 || !NameIsLegal(e$h5_path))
    stop('Invalid value of "h5_path": not a legal path name')
  if (file.exists(e$h5_path) && e$ask_overwrite) {
    cat(crayon::bgRed('     ')); .msg(' (?) '); cat(crayon::bgRed('     \n'))
    response <- readline(prompt = paste0('HDF5 file "', e$h5_path, '" already exists. Overwrite it? (Y/n) '))
    if (response != 'Y') { .msg_alt_bad('Aborting Benchmark object construction, returning NULL\n'); stop() }
  }
  
  if (e$input_class == 'fcs_paths') {
    fe <- file.exists(e$input)
    if (any(!fe)) stop(paste0('The following input files were not found:\n', paste0('\t', e$input[!fe], collapse = ',\n')))
    not_fcs <- sapply(e$input, function(f) !IsFormat(f, 'fcs'))
    if (any(!fe)) stop(paste0('The following input files are not FCS files:\n', paste0('\t', e$input[not_fcs], collapse = ', \n')))
  }
  
  if (!is.null(e$knn.algorithm) && !e$knn.algorithm %in% c('annoy', 'cover_tree', 'kd_tree', 'brute'))
    stop('Invalid k-NN algorithm')
}

.AddLayout.ValidityChecks <- function(e) {
  
  if (!e$benchmark$evaluated_previously)
    stop('No valid evaluation results available')
  if (!is.null(e$idx.subpipeline))
    if (round(e$idx.subpipeline) != e$idx.subpipeline || e$idx.subpipeline < 1 || e$idx.subpipeline > e$benchmark$n_subpipelines)
      stop('Invalid "idx.subpipeline" value')
  if (!file.exists(e$benchmark$h5_path))
    stop(paste0('Auxiliary HDF5 file ', e$benchmark$h5_path, ' not found'))
  
  if (!is.null(e$method) && e$method$type != 'projection')
    stop('"method" is not of type "projection"')
  
  if (e$benchmark$layout_available && e$ask_overwrite) {
    response <- readline(prompt = 'A 2-dimensional layout is available already. Overwrite existing layout? (Y/n) ')
    if (response != 'Y')
      stop()
  }
}

.PlotJaccardHeatmap.ValidityChecks <- function(e) {
  
  if (!e$benchmark$evaluated_previously)
    stop('No valid evaluation results available')
  if (round(e$idx.subpipeline) != e$idx.subpipeline || e$idx.subpipeline < 1 || e$idx.subpipeline > e$benchmark$n_subpipelines)
    stop('Invalid "idx.subpipeline" value')
}