
.h5write <- function(
  ## Handles NULLs
  obj, file, name
) {
  if (is.null(obj))
    rhdf5::h5write(NA, file, name)
  else
    rhdf5::h5write(obj, file, name)
}
  

.h5read <- function(
  file, name
) {
  obj <- rhdf5::h5read(file, name)
  if (is.atomic(obj) && is.na(obj))
    NULL
  else
    obj
}

.h5writeNamedMatrix <- function(
  obj, file, name
) {
  rhdf5::h5createGroup(file, name)
  
  rhdf5::h5write(obj, file, paste0(name, '/Matrix'))
  if (!is.null(colnames(obj)))
    rhdf5::h5write(colnames(obj), file, paste0(name, '/ColumnNames'))
  if (!is.null(rownames(obj)))
    rhdf5::h5write(rownames(obj), file, paste0(name, '/RowNames'))
}

.h5readNamedMatrix <- function(
  file, name
) {
  group <- rhdf5::h5read(file, name)
  obj <- group$Matrix
  if (!is.null(group$ColumnNames))
    colnames(obj) <- group$ColumnNames
  if (!is.null(group$RowNames))
    rownames(obj) <- group$RowNames
  obj
}

.h5writeFactorVectorOrListOfThem <- function( # recursive
  obj, file, name
) {
  rhdf5::h5createGroup(file, name)
  if (!is.list(obj)) {
    rhdf5::h5write(NA, file, paste0(name, '/NSamples'))
    rhdf5::h5write(as.integer(obj), file, paste0(name, '/Vector'))
    rhdf5::h5write(levels(obj), file, paste0(name, '/Levels'))
  } else {
    rhdf5::h5write(length(obj), file, paste0(name, '/NSamples'))
    for (idx_input in 1:length(obj)) {
      input_group_name <- paste0(name, '/Sample', formatC(idx_input, width = 4, flag = '0'))
      .h5writeFactorVectorOrListOfThem(obj[[idx_input]], file, input_group_name)
    }
  }
}

.h5readFactorVectorOrListOfThem <- function( # recursive
  file, name, concatenate = TRUE
) {
  nsamp <- rhdf5::h5read(file, paste0(name, '/NSamples'))
  if (is.na(nsamp)) {
    obj <- as.factor(rhdf5::h5read(file, paste0(name, '/Vector')))
    levels(obj) <- rhdf5::h5read(file, paste0(name, '/Levels'))
  } else {
      obj <- purrr::map(1:nsamp, function(idx_input) .h5readFactorVectorOrListOfThem(file, paste0(name, '/Sample', formatC(idx_input, width = 4, flag = '0'))))
      if (concatenate)
        obj <- as.factor(do.call(c, purrr::map(obj, as.character)))
  }
  obj
}

.h5_slotname <- function(
  idx.subpipeline = NULL,
  tool_type = NULL,
  idx.n_param = NULL,
  idx.module = NULL,
  prefix = 'EvaluationResults',
  suffix = NULL
) {
  
  prefix                     <- if (!is.null(prefix)) paste0('/', prefix) else NULL
  substring_subpipeline      <- if (!is.null(idx.subpipeline)) paste0('/SubPipeline', formatC(idx.subpipeline, width = 4, flag = '0')) else NULL
  substring_subpipeline_type <- if (!is.null(tool_type)) paste0('/', tool_type) else NULL
  substring_nparam           <- if (!is.null(idx.n_param)) paste0('/NParam', formatC(idx.n_param, width = 4, flag = '0')) else NULL
  substring_module           <- if (!is.null(idx.module)) paste0('/Module', formatC(idx.n_param, width = 4, flag = '0')) else NULL
  suffix                     <- if (!is.null(suffix)) paste0('/', suffix) else NULL
  
  s <- paste0(prefix, substring_subpipeline, substring_subpipeline_type, substring_nparam, substring_module, suffix)
  gsub('//', '/', s, fixed = TRUE)
}

.h5writeProjectionResult <- function(
  obj = NA, benchmark, idx.subpipeline, idx.n_param = NULL
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Projection', idx.n_param)
  suppressMessages(.h5writeNamedMatrix(obj$Projection, benchmark$h5_path, paste0(slotname, '/Result')))
  rhdf5::h5write(obj$Timing, benchmark$h5_path, paste0(slotname, '/Timing'))
}

.h5writeProjectionIntermediate <- function(
  obj, h5_path, idx.subpipeline, idx.n_param, idx.module
) {
  slotname <- .h5_slotname(
    idx.subpipeline = idx.subpipeline,
    tool_type       = 'Projection',
    idx.n_param     = idx.n_param,
    idx.module      = idx.module
  )
  suppressMessages(.h5writeNamedMatrix(obj, h5_path, paste0(slotname, '/Intermediate')))
}

.h5writeClusteringIntermediate <- function(
  obj, benchmark, idx.subpipeline, idx.n_param, idx.module
) {
  slotname <- .h5_slotname(
    idx.subpipeline = idx.subpipeline,
    tool_type       = 'Clustering',
    idx.n_param     = idx.n_param,
    idx.module      = idx_module,
    prefix          = 'EvaluationResults'
  )
  suppressMessages(.h5writeNamedMatrix(obj, benchmark$h5_path, paste0(slotname, '/ClusteringIntermediate')))
}

.h5writeProjectionReference <- function(
  benchmark, idx.subpipeline, idx.n_param = NULL, idx.subpipeline_ref, idx.n_param_ref = NULL
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Projection', idx.n_param)
  
  rhdf5::h5delete(benchmark$h5_path, paste0(slotname, '/IsReferenceToSubpipeline'))
  rhdf5::h5write(idx.subpipeline_ref, benchmark$h5_path, paste0(slotname, '/IsReferenceToSubpipeline'))
  
  if (!is.null(idx.n_param) && !is.null(idx.n_param_ref)) {
    rhdf5::h5delete(benchmark$h5_path, paste0(slotname, '/IsReferenceToNParamIteration'))
    rhdf5::h5write(idx.n_param_ref, benchmark$h5_path, paste0(slotname, '/IsReferenceToNParamIteration'))
  }
}

.h5readProjectionResult <- function(
  benchmark, idx.subpipeline, idx.n_param = NULL
) {
  
  idx_subpipeline_ref <- rhdf5::h5read(benchmark$h5_path, paste0(.h5_slotname(idx.subpipeline, 'Projection', idx.n_param), '/IsReferenceToSubpipeline'))
  idx_nparam_ref <- if (!is.null(idx.n_param)) rhdf5::h5read(benchmark$h5_path, paste0(.h5_slotname(idx.subpipeline, 'Projection', idx.n_param), '/IsReferenceToNParamIteration')) else NA
  
  if (is.na(idx_subpipeline_ref) || idx_subpipeline_ref < 0) {
    idx_subpipeline_ref <- idx.subpipeline
  
  } else if (idx_subpipeline_ref == 0) {
    return(GetExpressionMatrix(benchmark))
    
  }
  
  if (is.na(idx_nparam_ref) || idx_nparam_ref < 0) {
    idx_nparam_ref <- idx.n_param
    
  } else if (idx_nparam_ref == 0) {
    return(GetExpressionMatrix(benchmark))
    
  }
  
  slotname <- .h5_slotname(idx_subpipeline_ref, 'Projection', idx_nparam_ref)
  list(
    Projection = .h5readNamedMatrix(benchmark$h5_path, paste0(slotname, '/Result')),
    Timing = rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Timing'))
  )
}

.h5writekNNMatrix <- function(
  benchmark, obj
) {
  idx_knn <- which(rhdf5::h5ls(benchmark$h5_path)$group == '/Input' & rhdf5::h5ls(benchmark$h5_path)$name == 'kNN')
  if (length(idx_knn) == 0)
    rhdf5::h5createGroup(benchmark$h5_path, '/Input/kNN')
  suppressMessages(.h5write(obj$Indices, benchmark$h5_path, '/Input/kNN/Indices'))
  suppressMessages(.h5write(obj$Distances, benchmark$h5_path, '/Input/kNN/Distances'))
}

.h5writeDistanceMatrix <- function(
  benchmark, obj
) {
  suppressMessages(.h5write(obj, benchmark$h5_path, '/Input/DistanceMatrix'))
}

.h5readkNNMatrix <- function(
  benchmark
) {
  list(
    Indices = .h5read(benchmark$h5_path, '/Input/kNN/Indices'),
    Distances = .h5read(benchmark$h5_path, '/Input/kNN/Distances')
  )
}

.h5readDistanceMatrix <- function(
  benchmark
) {
    Indices = .h5read(benchmark$h5_path, '/Input/DistanceMatrix')
}

.h5writeProjectionScoring <- function(
  obj, benchmark, idx.subpipeline, idx.n_param
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Projection', idx.n_param)
  idx_ref <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/IsReferenceToSubpipeline'))
  if (is.na(idx_ref) || idx_ref < 0) {
    rhdf5::h5createGroup(benchmark$h5_path, paste0(slotname, '/Scoring'))
    .h5write(obj[['Layout k-NNG']]$Indices, benchmark$h5_path, paste0(slotname, '/Scoring/LayoutkNNGIndices'))
    .h5write(obj[['Layout k-NNG']]$Distances, benchmark$h5_path, paste0(slotname, '/Scoring/LayoutkNNGDistances'))
    .h5write(obj[['Collapsed']], benchmark$h5_path, paste0(slotname, '/Scoring/Collapsed'))
    .h5write(obj[['Co-Ranking Matrix']], benchmark$h5_path, paste0(slotname, '/Scoring/CoRankingMatrix'))
    .h5write(obj[['Local Continuity Meta-Criterion']], benchmark$h5_path, paste0(slotname, '/Scoring/LocalContinuityMetaCriterion'))
    .h5write(obj[['Relative Intrusiveness']], benchmark$h5_path, paste0(slotname, '/Scoring/RelativeIntrusiveness'))
    .h5write(obj[['Trustworthiness']], benchmark$h5_path, paste0(slotname, '/Scoring/Trustworthiness'))
    .h5write(obj[['Continuity']], benchmark$h5_path, paste0(slotname, '/Scoring/Continuity'))
    .h5write(obj[['Projection Neighbourhood']], benchmark$h5_path, paste0(slotname, '/Scoring/ProjectionNeighbourhood'))
  }
}

.h5readProjectionScoring <- function(
  benchmark, idx.subpipeline, idx.n_param = NULL
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Projection', idx.n_param)
  idx.subpipeline_ref <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/IsReferenceToSubpipeline'))
  idx.n_param_ref <-
    if (!is.null(idx.n_param))
      rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/IsReferenceToNParamIteration'))
    else
      NULL
  if (is.na(idx.n_param_ref) || idx.n_param_ref < 0)
    idx.n_param_ref <- NULL
  if (!is.na(idx.subpipeline_ref) && idx.subpipeline_ref > 0)
    slotname <- .h5_slotname(idx.subpipeline_ref, 'Projection', idx.n_param_ref)
  
  obj <- list()
  obj[['Layout k-NNG']] <- list(
    Indices = rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/LayoutkNNGIndices')),
    Distances = rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/LayoutkNNGDistances'))
  )
  obj[['Collapsed']] <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/Collapsed'))
  obj[['Co-Ranking Matrix']] <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/CoRankingMatrix'))
  obj[['Local Continuity Meta-Criterion']] <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/LocalContinuityMetaCriterion'))
  obj[['Relative Intrusiveness']] <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/RelativeIntrusiveness'))
  obj[['Trustworthiness']] <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/Trustworthiness'))
  obj[['Continuity']] <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/Continuity'))
  obj[['Projection Neighbourhood']] <- rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring/ProjectionNeighbourhood'))
  
  obj
}

.h5readClusteringInput <- function(
  benchmark, idx.subpipeline, idx.n_param = NULL, null_if_exprs = FALSE
) {
  proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
  if (!is.null(proj)) {
    res <- GetProjection(benchmark, idx.subpipeline, idx.n_param)
    if (is.list(res)) return(res$Projection) else return(res)
  } else {
    if (null_if_exprs)
      return(NULL)
    else
      return(GetExpressionMatrix(benchmark))
  }
}

.h5writeClusteringResult <- function(
  obj, benchmark, idx.subpipeline, idx.n_param
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Clustering', idx.n_param)
  rhdf5::h5write(obj$ClusteringVector, benchmark$h5_path, paste0(slotname, '/ClusteringVector'))
  rhdf5::h5write(obj$Timing, benchmark$h5_path, paste0(slotname, '/Timing'))
}

.h5readClusteringResult <- function(
  benchmark, idx.subpipeline, idx.n_param
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Clustering', idx.n_param)
  list(
    ClusteringVector = rhdf5::h5read(benchmark$h5_path, name = paste0(slotname, '/ClusteringVector')),
    Timing = rhdf5::h5read(benchmark$h5_path, name = paste0(slotname, '/Timing'))
  )
}

.h5writeAnnotationScoring <- function(
  obj, benchmark, suffix = ''
) {
  slotname <- paste0('/Input/AnnotationScore', suffix)
  rhdf5::h5write(obj, benchmark$h5_path, slotname)
}

.h5readAnnotationScoring <- function(
  benchmark, suffix = ''
) {
  slotname <- paste0('/Input/AnnotationScore', suffix)
  rhdf5::h5read(benchmark$h5_path, slotname)
}

.h5writeClusteringScoring <- function(
  obj, benchmark, idx.subpipeline, idx.n_param
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Clustering', idx.n_param, suffix = 'Scoring')
  rhdf5::h5write(obj, benchmark$h5_path, slotname)
}

.h5readClusteringScoring <- function(
  benchmark, idx.subpipeline, idx.n_param
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Clustering', idx.n_param)
  rhdf5::h5read(benchmark$h5_path, paste0(slotname, '/Scoring'))
}

.h5writeClusterMedians <- function(
  obj, benchmark, idx.subpipeline, idx.n_param
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Clustering', idx.n_param, suffix = 'MedianExpressionPerCluster')
  .h5writeNamedMatrix(obj, benchmark$h5_path, slotname)
}

.h5readClusterMedians <- function(
  benchmark, idx.subpipeline, idx.n_param
) {
  slotname <- .h5_slotname(idx.subpipeline, 'Clustering', idx.n_param)
  .h5readNamedMatrix(benchmark$h5_path, paste0(slotname, '/MedianExpressionPerCluster'))
}

.h5writeLabelMedians <- function(
  obj, benchmark
) {
  .h5writeNamedMatrix(obj, benchmark$h5_path, '/Input/MedianExpressionPerPopulation')
}

.h5readLabelMedians <- function(
  benchmark
) {
  .h5readNamedMatrix(benchmark$h5_path, '/Input/MedianExpressionPerPopulation')
}
