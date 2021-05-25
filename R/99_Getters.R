
#' Get \code{Benchmark} collapsed co-ranking matrix
#'
#' Extracts the collapsed co-ranking matrix for assessing quality of a projection.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer value: index of subpipeline that includes a projection step
#' @param idx.n_param optional integer value: index of subpipeline n-parameter iteration. Default value is \code{NULL}
#'
#' @seealso
#'
#' @export
GetCoRanking <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL
) {
  list(
    Matrix = .h5read(benchmark$h5_path, .h5_slotname(idx.subpipeline, 'Projection', idx.n_param, suffix = 'Scoring/CoRankingMatrix')),
    Collapsed = .h5read(benchmark$h5_path, .h5_slotname(idx.subpipeline, 'Projection', idx.n_param, suffix = 'Scoring/Collapsed')),
    K = .h5read(benchmark$h5_path, .h5_slotname(idx.subpipeline, 'Projection', idx.n_param, suffix = 'Scoring/ProjectionNeighbourhood'))
  )
}

#' Get \code{Benchmark} penalty scoring matrix
#'
#' Extracts the penalty scoring matrix, which defines a hierarchy penalty model, from a \code{Benchmark}-type object.
#' 
#' @param benchmark object of type \code{Benchmark} that uses hierarchical penalties
#'
#' @seealso
#' 
#' * **\code{CreatePenaltyScoringMatrix}**: creates a penalty scoring matrix for a given set of manually annotated populations
#'
#' @export
GetPenaltyScoringMatrix <- function(
  benchmark
) {
  .h5readNamedMatrix(benchmark, '/Input/Annotation/ScoringMatrix')
}

#' Get \code{Benchmark} manual annotation scoring
#'
#' Extracts values of unsupervised evaluation metrics applied to manual annotation of benchmark input data.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param suffix optional string (advanced): HDF5 slot name suffix
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetAnnotationScoring <- function(
  benchmark, suffix = ''
) {
  .h5readAnnotationScoring(benchmark, suffix)
}

#' Get \code{Benchmark} expression data
#'
#' Extracts expression matrix (or list of expression matrices) associated with inputs to a benchmark pipeline.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param concatenate logical value; if a list of layout matrices (per sample) is available, should it be concatenated into one vector? Default value is \code{FALSE}
#'
#' @seealso 
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetExpressionMatrix <- function(
  benchmark, concatenate = FALSE
) {
  obj <- .h5read(benchmark$h5_path, 'Input/ExpressionMatrix')
  if (is.list(obj)) {
    cn <- .h5read(benchmark$h5_path, 'Input/ColumnNames')
    purrr::walk(seq_len(obj), function(idx) colnames(obj[[idx]]) <- cn)
  } else {
    colnames(obj) <- .h5read(benchmark$h5_path, 'Input/ColumnNames')
  }
  if (concatenate && is.list(obj))
    obj <- do.call(rbind, obj)
  obj
}

#' Get \code{Benchmark} expression data \code{k}-NN matrix
#'
#' Extracts *k*-nearest-neighbours matrix associated with inputs to a benchmark pipeline.
#' 
#' @param benchmark object of type \code{Benchmark}
#'
#' @seealso 
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetkNNMatrix <- function(
  benchmark
) {
  .h5readkNNMatrix(benchmark)
}

#' Get \code{Benchmark} expression data distance matrix
#'
#' Extracts distance matrix associated with inputs to a benchmark pipeline.
#' 
#' @param benchmark object of type \code{Benchmark}
#'
#' @seealso 
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetDistanceMatrix <- function(
  benchmark
) {
  .h5readDistanceMatrix(benchmark)
}

#' Get \code{Benchmark} manual annotation
#'
#' Extracts manual annotation vector (or list of vectors) associated with inputs to a benchmark pipeline.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param concatenate logical value; if a list of annotation vectors (per sample) is available, should it be concatenated into one vector? Default value is \code{FALSE}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetAnnotation <- function(
  benchmark, concatenate = FALSE
) {
  .h5readFactorVectorOrListOfThem(benchmark$h5_path, 'Input/Annotation', concatenate = concatenate)
}

#' Get \code{Benchmark} sub-pipeline name
#'
#' Extracts name of a sub-pipeline in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a sub-pipeline of \code{benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetSubpipelineName <- function(
  benchmark, idx.subpipeline
) {
  GetSubpipelineTags(benchmark$subpipelines, idx.subpipeline)$subpipeline
}

#' Get \code{Benchmark} latent-space projection
#'
#' Extracts coordinate matrix of a latent-space projection generated by evaluating an object of type \code{Benchmark}.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a sub-pipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#' 
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#'
#' @export
GetProjection <- function(
  benchmark, idx.subpipeline, idx.n_param = NULL
) {
  .h5readProjectionResult(benchmark, idx.subpipeline, idx.n_param)$Projection
}

#' Get \code{Benchmark} 2-dimensional layout of input data
#'
#' Extracts coordinate matrix of a 2-dimensional layout associated with the input data for a \code{Benchmark} object (benchmark pipeline set-up).
#' This layout can be created by applying \code{AddLayout} to a previously evaluated benchmark pipeline and serves visualisation purposes.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param concatenate logical value; if a list of layout matrices (per sample) is available, should it be concatenated into one vector? Default value is \code{FALSE}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#' 
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#' 
#' * **\code{AddLayout}**: allows you to add a separate 2-dimensional layout of the input dataset or to use an existing projection (produced in the evaluation) as a visualisation layout.
#'
#' @export
GetLayout <- function(
  benchmark, concatenate = FALSE
) {
  if (!benchmark$layout_available) return(NULL)
  ref <- rhdf5::h5read(benchmark$h5_path, 'Input/Layout/IsReferenceTo')
  if (length(ref) == 2 && !any(is.na(ref)) && is.numeric(ref)) {
    obj <- GetProjection(benchmark, ref[1], ref[2])
  } else {
    obj <- rhdf5::h5read(benchmark$h5_path, 'Input/Layout/Coordinates')
  }
  if (is.list(obj)) {
    obj <- purrr::map(obj, function(x) { x[, 1:2]; colnames(x) <- c('Component1', 'Component2') })
  } else {
    obj <- obj[, 1:2]
    colnames(obj) <- c('Component1', 'Component2')
  }
    
  if (is.list(obj) && concatenate)
    obj <- do.call(rbind, obj)
  
  obj
}

#' Get \code{Benchmark} bootstrap indices
#'
#' Extracts vector of indices (or list of vectors of indices) used for taking bootstraps of input expression data for clustering (if this is done).
#' 
#' @param benchmark object of type \code{Benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#'
#' @export
GetBootstrapIndices <- function(
  benchmark
) {
  .h5read(benchmark$h5_path, 'Input/BootstrapIndices')
}

#' Get \code{Benchmark} clustering input
#'
#' Extracts the input to a clustering tool that is part of a benchmark sub-pipeline that was set up previously.
#' If a dimension-reduction tool precedes the clustering tool, the latent-space projection from that tool is retrieved.
#' Otherwise, the original expression data is retrieved.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a sub-pipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#' @param null_if_exprs logical: whether to return \code{NULL} instead of original expression data. Default value is \code{FALSE}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#' 
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#'
#' @export
GetClusteringInput <- function(
  benchmark, idx.subpipeline, idx.n_param = NULL, null_if_exprs = FALSE
) {
  n_npar <- GetNParameterIterationsCount(benchmark, idx.subpipeline)
  if (n_npar == 0)
    idx.n_param <- NULL
  .h5readClusteringInput(benchmark, idx.subpipeline, idx.n_param, null_if_exprs = null_if_exprs)
}

#' Get \code{Benchmark} clustering result
#'
#' Extracts the result of a clustering tool that is part of a benchmark sub-pipeline that was set up previously.
#' If stability analysis is turned on for the benchmark, either results of all runs of the clustering can be returned, or only the result of a single run on the full expression data can be returned.
#' 
#' If you choose to extract results of all runs of a clustering algorithm which clustered on bootstraps of the original data, the last slot in the list contains the single run on original data.
#' Furthermore, you can extract bootstrap indices using \code{GetBootstrapIndices}.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a sub-pipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#' @param all_runs logical: whether results of all runs should be returned (in a list). Default value is \code{FALSE}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#' 
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#'
#' @export
GetClustering <- function(
  benchmark, idx.subpipeline, idx.n_param = NULL, all_runs = FALSE
) {
  obj <- .h5readClusteringResult(benchmark, idx.subpipeline, idx.n_param)$ClusteringVector
  if (!all_runs) {
    if (benchmark$stability == 'repeat')
      obj <- obj[[1]]
    else if (benchmark$stability == 'bootstrap')
      obj <- obj[[length(obj)]]
  }
  obj
}

#' Get \code{Benchmark} *n*-parameter iteration name for a projection step
#'
#' Extracts name of an *n*-parameter iteration of the projection step of a sub-pipeline in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a subpipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#' @param with_tool_name logical: whether tool name, in addition to *n*-parameter value, should be included in the name. Default value is \code{TRUE}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetNParameterIterationName_Projection <- function(
  benchmark, idx.subpipeline, idx.n_param, with_tool_name = TRUE
) {
  
  proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
  
  if (is.null(proj))
    return(NULL)
  
  if (IsClone(proj)) {
    idx.subpipeline <- proj$ref
    proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
  }
  
  tool_name <- proj$name
  
  if (is.null(idx.n_param) || idx.n_param == 'NoNParameter')
    return(tool_name)
  
  n_param_name <- proj$modules[[proj$which_n_param]]$n_param
  n_param_val <- benchmark$n_params[[idx.subpipeline]]$projection[idx.n_param]
  
  string_n_param <- paste0(n_param_name, '=', n_param_val)
  
  if (with_tool_name)
    return(paste0(tool_name, ' [', string_n_param, ']'))
  else
    return(string_n_param)
}

GetProjectionModuleCount <- function(benchmark, idx.subpipeline) {
  if (!is.null(benchmark$subpipelines[[idx.subpipeline]]$projection))
    benchmark$subpipelines[[idx.subpipeline]]$projection$n_modules
  else
    NULL
}

GetClusteringModuleCount <- function(benchmark, idx.subpipeline) {
  if (!is.null(benchmark$subpipelines[[idx.subpipeline]]$clustering))
    benchmark$subpipelines[[idx.subpipeline]]$clustering$n_modules
  else
    NULL
}

#' Get \code{Benchmark} *n*-parameter iteration clustering name
#'
#' Extracts name of an *n*-parameter iteration of the clustering part of a subpipeline in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a subpipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#' @param with_tool_name logical: whether tool name, in addition to *n*-parameter value, should be included in the name. Default value is \code{TRUE}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetNParameterIterationName_Clustering <- function(
  benchmark, idx.subpipeline, idx.n_param, with_tool_name = TRUE
) {
  
  clus <- benchmark$subpipelines[[idx.subpipeline]]$clustering
  
  if (is.null(clus))
    return(NULL)
  
  tool_name <- clus$name
  
  if (is.null(idx.n_param) || is.null(clus$which_n_param))
    return(tool_name)
  
  n_param_name <- clus$modules[[clus$which_n_param]]$n_param
  n_param_val <- benchmark$n_params[[idx.subpipeline]]$clustering[idx.n_param]
  
  string_n_param <- paste0(n_param_name, '=', n_param_val)
  
  if (with_tool_name)
    return(paste0(tool_name, ' [', string_n_param, ']'))
  else
    return(string_n_param)
}

#' Get \code{Benchmark} *n*-parameter values for a subpipeline
#'
#' Extracts combinations of *n*-parameter values for each *n*-parameter iteration of a subpipeline in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a subpipeline of \code{benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetNParameterValues <- function(
  benchmark,
  idx.subpipeline
) {
  
  npar_proj <-
    if (IsClone(benchmark$subpipelines[[idx.subpipeline]]$projection))
      benchmark$n_params[[benchmark$subpipelines[[idx.subpipeline]]$projection$ref]]$projection
    else
      benchmark$n_params[[idx.subpipeline]]$projection
  npar_clus <- benchmark$n_params[[idx.subpipeline]]$clustering
  
  n <- max(length(npar_proj), length(npar_clus))
  if (is.null(npar_proj)) npar_proj <- rep(NA, n)
  if (is.null(npar_clus)) npar_clus <- rep(NA, n)
  
  data.frame(
    'idx.n_param' = seq_len(n),
    'npar_proj' = npar_proj,
    'npar_clus' = npar_clus
  )
}

GetNParameterNames <- function(
  benchmark,
  idx.subpipeline,
  with_tool_names = FALSE
) {
  
  proj <- benchmark$subpipelines[[idx.subpipeline]]$projection
  clus <- benchmark$subpipelines[[idx.subpipeline]]$clustering
  
  if (IsClone(proj))
    proj <- benchmark$subpipelines[[proj$ref]]$projection
  
  proj_n_param <- if (!is.null(proj) && !is.null(proj$which_n_param)) proj$modules[[proj$which_n_param]]$n_param else NULL
  proj_idx     <- if (!is.null(proj) && !is.null(proj$which_n_param)) proj$which_n_param else NULL
  proj_tool    <- if (!is.null(proj) && !is.null(proj$which_n_param)) proj$modules[[proj$which_n_param]]$name else NULL
  
  
  clus_n_param <- if (!is.null(clus) && !is.null(clus$which_n_param)) clus$modules[[clus$which_n_param]]$n_param else NULL
  clus_idx     <- if (!is.null(clus) && !is.null(clus$which_n_param)) clus$which_n_param else NULL
  clus_tool    <- if (!is.null(clus) && !is.null(clus$which_n_param)) clus$modules[[clus$which_n_param]]$name else NULL
  
  list(
    projection = if (is.null(proj_n_param)) { NULL } else { if (with_tool_names) paste0(proj_tool, ' (module ', proj_idx, ') ', proj_n_param) else proj_n_param },
    clustering = if (is.null(clus_n_param)) { NULL } else { if (with_tool_names) paste0(clus_tool, ' (module ', clus_idx, ') ', clus_n_param) else clus_n_param }
  )
}

#' Get \code{Benchmark} *n*-parameter iteration name
#'
#' Extracts name of an *n*-parameter iteration of a subpipeline in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a subpipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#' @param with_tool_names logical: whether tool name(s), in addition to *n*-parameter value(s), should be included in the name. Default value is \code{TRUE}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetNParameterIterationName <- function(
  benchmark, idx.subpipeline, idx.n_param, with_tool_names = TRUE
) {
  if (is.null(idx.n_param))
    return(GetSubpipelineName(benchmark, idx.subpipeline))
  name_proj <- GetNParameterIterationName_Projection(benchmark, idx.subpipeline, idx.n_param, with_tool_names)
  name_clus <- GetNParameterIterationName_Clustering(benchmark, idx.subpipeline, idx.n_param, with_tool_names)
  
  if (is.null(name_proj))
    name_clus
  else if (is.null(name_clus))
    name_proj
  else if (is.null(name_proj) && is.null(name_clus))
    NULL
  else
    paste0(name_proj, ' -> ', name_clus)
}

#' Get \code{Benchmark} *n*-parameter count for subpipeline
#'
#' Extract the number of *n*-parameter iterations for a chosen subpipeline in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a sub-pipeline of \code{benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' @export
GetNParameterIterationsCount <- function(
  benchmark, idx.subpipeline
) {
  if (!is.null(benchmark$n_params[[idx.subpipeline]]$clustering))
    length(benchmark$n_params[[idx.subpipeline]]$clustering)
  else if (!is.null(benchmark$n_params[[idx.subpipeline]]$projection))
    length(benchmark$n_params[[idx.subpipeline]]$projection)
  else if (IsClone(benchmark$subpipelines[[idx.subpipeline]]$projection)) {
    length(benchmark$n_params[[benchmark$subpipelines[[idx.subpipeline]]$projection$ref]]$projection)
  } else
    0
}

#' Get \code{Benchmark} dimension-reduction evaluation scores
#'
#' Extracts list of evaluation scores and other products of evaluating performance of a projection step of a sub-pipeline *n*-iteration that includes dimension reduction.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a subpipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#'
#' @export
GetProjectionScores <- function(
  benchmark, idx.subpipeline, idcs.n_param = NULL
) {
  
  if (is.null(idcs.n_param))
    idcs.n_param <- seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline))
  
  no_npar <- FALSE
  if (length(idcs.n_param) == 0) {
    no_npar <- TRUE
    idcs.n_param <- 1
  }
  
  proj_scores <- NULL
  if (!is.null(benchmark$subpipelines[[idx.subpipeline]]$projection)) {
    proj_scores <- purrr::map(
      idcs.n_param,
      function(idx.n_param) .h5readProjectionScoring(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
    names(proj_scores) <- purrr::map_chr(
      idcs.n_param,
      function(idx.n_param) GetNParameterIterationName_Projection(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
  }
  proj_scores
}

#' Get \code{Benchmark} clustering evaluation scores
#'
#' Extracts list of evaluation scores and other products of evaluating performance of a clustering tool used in a subpipeline *n*-iteration in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a subpipeline of \code{benchmark}
#' @param idx.n_param integer: index of an *n*-parameter iteration of a sub-pipeline of \code{benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#'
#' @export
GetClusteringScores <- function(
  benchmark, idx.subpipeline, idcs.n_param = NULL
) {
  if (is.null(idcs.n_param))
    idcs.n_param <- seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline))
  
  no_npar <- FALSE
  if (length(idcs.n_param) == 0) {
    no_npar <- TRUE
    idcs.n_param <- 1
  }
  
  clus_scores <- NULL
  if (!is.null(benchmark$subpipelines[[idx.subpipeline]]$clustering)) {
    clus_scores <- purrr::map(
      idcs.n_param,
      function(idx.n_param) .h5readClusteringScoring(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
    names(clus_scores) <- purrr::map_chr(
      idcs.n_param,
      function(idx.n_param) GetNParameterIterationName_Clustering(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
  }
  clus_scores
}

#' Get \code{Benchmark} evaluation scores of an *n*-parameter iteration
#'
#' Extracts list of evaluation scores and other products of evaluating performance of tools used in a subpipeline *n*-iteration in a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a sub-pipeline of \code{benchmark}
#' @param idcs.n_param integer or vector of integers: indices of *n*-parameter iterations of a subpipeline of \code{benchmark}
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#'
#' @export
GetScores <- function(
  benchmark, idx.subpipeline, idcs.n_param = NULL
) {
  if (is.null(idcs.n_param))
    idcs.n_param <- seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline))
    
  no_npar <- FALSE
  if (length(idcs.n_param) == 0) {
    no_npar <- TRUE
    idcs.n_param <- 1
  }
  
  proj_scores <- NULL
  if (!is.null(benchmark$subpipelines[[idx.subpipeline]]$projection) && benchmark$score_projections) {
    proj_scores <- purrr::map(
      idcs.n_param,
      function(idx.n_param) .h5readProjectionScoring(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
    names(proj_scores) <- purrr::map_chr(
      idcs.n_param,
      function(idx.n_param) GetNParameterIterationName_Projection(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
  }
  clus_scores <- NULL
  if (!is.null(benchmark$subpipelines[[idx.subpipeline]]$clustering)) {
    clus_scores <- purrr::map(
      idcs.n_param,
      function(idx.n_param) .h5readClusteringScoring(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
    names(clus_scores) <- purrr::map_chr(
      idcs.n_param,
      function(idx.n_param) GetNParameterIterationName_Clustering(benchmark, idx.subpipeline, if (no_npar) NULL else idx.n_param)
    )
  }
  list(
    Projection = proj_scores,
    Clustering = clus_scores
  )
}

#' Get \code{Benchmark} projection evaluation scores table
#'
#' Returns a \code{tibble} of projection evaluation metric values for specified subpipeline of a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a subpipeline of \code{benchmark} that includes a projection step
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' * **\code{Evaluate}**: runs all benchmark subpipelines and scores the performance of each tool
#'
#' @export
GetProjectionScoringTable <- function(
  benchmark,
  idx.subpipeline
) {
  scores <- GetProjectionScores(benchmark, idx.subpipeline)
  
  n_param_names  <- GetNParameterNames(benchmark, idx.subpipeline)
  n_param_values <- GetNParameterValues(benchmark, idx.subpipeline)
  
  lcmc <- purrr::map(scores, function(x) x$`Local Continuity Meta-Criterion`)
  Bnx <- purrr::map(scores, function(x) x$`Relative Intrusiveness`)
  eT <- purrr::map(scores, function(x) x$`Trustworthiness`)
  eC <- purrr::map(scores, function(x) x$`Continuity`)
  
  names(lcmc) <- names(Bnx) <- names(eT) <- names(eC) <- names(scores)
  
  full <- any(!is.na(unlist(eT)))
  
  if (full) {
    metrics_names <- c(
      'Local Continuity Meta-Criterion',
      'Relative Intrusiveness',
      'Trustworthiness',
      'Continuity'
    )
  } else {
    metrics_names <- c(
      'Local Continuity Meta-Criterion',
      'Relative Intrusiveness'
    )
  }
  
  if (nrow(n_param_values) == 0) {
    vals <- if (full) unlist(c(lcmc, Bnx, eT, eC)) else unlist(c(lcmc, Bnx))
    vals <- dplyr::tibble(
      'Evaluation Metric' = metrics_names,
      'Value' = vals
    )
    return(vals)
  } else {
    vals_by_npar <-
      if (full)
        dplyr::tibble(
          'npar_proj' = n_param_values$npar_proj,
          unlist(lcmc),
          unlist(Bnx),
          unlist(eT),
          unlist(eC)
        )
      else
        dplyr::tibble(
          'npar_proj' = n_param_values$npar_proj,
          unlist(lcmc),
          unlist(Bnx)
        )
    vals_by_npar <- vals_by_npar[, apply(vals_by_npar, 2, function(x) any(!is.na(x)))]
    colnames(vals_by_npar) <- c(n_param_names$projection, metrics_names)
    vals_by_npar <- tidyr::pivot_longer(vals_by_npar, cols = metrics_names, names_repair = 'minimal')
    colnames(vals_by_npar)[(ncol(vals_by_npar) - 1):ncol(vals_by_npar)] <- c('Evaluation Metric', 'Value')
    vals_by_npar$`Evaluation Metric` <- as.factor(vals_by_npar$`Evaluation Metric`)
    return(vals_by_npar)
  }
}

#' Get \code{Benchmark} clustering evaluation scores table
#'
#' Returns a \code{tibble} of clustering evaluation metric values for specified subpipeline of a \code{Benchmark} object (benchmark pipeline set-up).
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of a sub-pipeline of \code{benchmark} that includes a clustering step
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)
#'
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#'
#' @export
GetClusteringScoringTable <- function(
  benchmark,
  idx.subpipeline
) {
  scores <- GetClusteringScores(benchmark, idx.subpipeline)
  
  n_param_names  <- GetNParameterNames(benchmark, idx.subpipeline)
  n_param_values <- GetNParameterValues(benchmark, idx.subpipeline)
  
  db           <- purrr::map(scores, function(x) unlist(x$`Davies-Bouldin Index`))
  ari          <- purrr::map(scores, function(x) unlist(x$`Adjusted Rand Index`))
  #nmi          <- purrr::map(scores, function(x) unlist(x$`Normalised Mutual Information`))
  f1_bij       <- purrr::map(scores, function(x) unlist(x$`Mean F1 Across Matches (Bijective)`))
  f1_fixed_cl  <- purrr::map(scores, function(x) unlist(x$`Mean F1 Across Matches (Relaxed, Fixed Cluster)`))
  f1_fixed_lab <- purrr::map(scores, function(x) unlist(x$`Mean F1 Across Matches (Relaxed, Fixed Label)`))
  pr_bij       <- purrr::map(scores, function(x) unlist(x$`Mean Precision Across Matches (Bijective)`))
  pr_fixed_cl  <- purrr::map(scores, function(x) unlist(x$`Mean Precision Across Matches (Relaxed, Fixed Cluster)`))
  pr_fixed_lab <- purrr::map(scores, function(x) unlist(x$`Mean Precision Across Matches (Relaxed, Fixed Label)`))
  re_bij       <- purrr::map(scores, function(x) unlist(x$`Mean Recall Across Matches (Bijective)`))
  re_fixed_cl  <- purrr::map(scores, function(x) unlist(x$`Mean Recall Across Matches (Relaxed, Fixed Cluster)`))
  re_fixed_lab <- purrr::map(scores, function(x) unlist(x$`Mean Recall Across Matches (Relaxed, Fixed Label)`))
  names(db) <- names(ari) <- #names(nmi) <- 
    names(f1_bij) <- names(f1_fixed_cl) <- names(f1_fixed_lab) <-
    names(pr_bij) <- names(pr_fixed_cl) <- names(pr_fixed_lab) <-
    names(re_bij) <- names(re_fixed_cl) <- names(re_fixed_lab) <- names(scores)
  metrics_names <- c(
    'Davies-Bouldin Index', 'Adjusted Rand Index', #'Normalised Mutual Information',
    'Mean F1 Bijective', 'Mean F1 Fixed-Cluster', 'Mean F1 Fixed-Label',
    'Mean Precision Bijective', 'Mean Precision Fixed-Cluster', 'Mean Precision Fixed-Label',
    'Mean Recall Bijective', 'Mean Recall Fixed-Cluster', 'Mean Recall Fixed-Label'
  )
  multiple_runs <- length(f1_bij[[1]]) > 1
  if (nrow(n_param_values) == 0) {
    vals        <- c(db, ari, #nmi,
                     f1_bij, f1_fixed_cl, f1_fixed_lab, pr_bij, pr_fixed_cl, pr_fixed_lab, re_bij, re_fixed_cl, re_fixed_lab)
    if (!multiple_runs) {
      vals <- dplyr::tibble(
        'Evaluation Metric' = metrics_names,
        'Value' = unlist(vals)
      )
    } else {
      vals <- dplyr::tibble(
        'Evaluation Metric' = metrics_names,
        'Mean Value' = sapply(vals, mean),
        'Standard Deviation' = sapply(vals, sd)
      )
    }
    return(vals)
  } else {
    if (!multiple_runs) {
      vals_by_npar <- dplyr::tibble(
        'npar_proj' = n_param_values$npar_proj,
        'npar_clus' = n_param_values$npar_clus,
        unlist(db), unlist(ari), #unlist(nmi),
        unlist(f1_bij), unlist(f1_fixed_cl), unlist(f1_fixed_lab), unlist(pr_bij), unlist(pr_fixed_cl), unlist(pr_fixed_lab), unlist(re_bij), unlist(re_fixed_cl), unlist(re_fixed_lab)
      )
      vals_by_npar <- vals_by_npar[, apply(vals_by_npar, 2, function(x) any(!is.na(x)))]
      colnames(vals_by_npar) <- c(n_param_names$projection, n_param_names$clustering, metrics_names)
      vals_by_npar <- tidyr::pivot_longer(vals_by_npar, cols = metrics_names, names_repair = 'minimal')
      colnames(vals_by_npar)[(ncol(vals_by_npar) - 1):ncol(vals_by_npar)] <- c('Evaluation Metric', 'Value')
    } else {
      means_by_npar <- dplyr::tibble(
        'npar_proj' = n_param_values$npar_proj,
        'npar_clus' = n_param_values$npar_clus,
        sapply(db, mean), sapply(ari, mean),
        #sapply(nmi, mean),
        sapply(f1_bij, mean), sapply(f1_fixed_cl, mean), sapply(f1_fixed_lab, mean), sapply(pr_bij, mean), sapply(pr_fixed_cl, mean), sapply(pr_fixed_lab, mean), sapply(re_bij, mean), sapply(re_fixed_cl, mean), sapply(re_fixed_lab, mean)
      )
      sd_by_npar <- dplyr::tibble(
        'npar_proj' = n_param_values$npar_proj,
        'npar_clus' = n_param_values$npar_clus,
        sapply(db, sd), sapply(ari, sd), #sapply(nmi, mean),
        sapply(f1_bij, sd), sapply(f1_fixed_cl, sd), sapply(f1_fixed_lab, sd), sapply(pr_bij, sd), sapply(pr_fixed_cl, sd), sapply(pr_fixed_lab, sd), sapply(re_bij, sd), sapply(re_fixed_cl, sd), sapply(re_fixed_lab, sd)
      )
      
      npar_proj <- n_param_names$projection
      if (is.null(npar_proj))
        npar_proj <- 'npar_proj'
      npar_clus <- n_param_names$clustering
      if (is.null(npar_clus))
        npar_clus <- 'npar_clus'
      
      colnames(means_by_npar) <- colnames(sd_by_npar) <- c(npar_proj, npar_clus, metrics_names)
      means_by_npar <- tidyr::pivot_longer(means_by_npar, cols = metrics_names)
      sd_by_npar <- tidyr::pivot_longer(sd_by_npar, cols = metrics_names)
      
      vals_by_npar <- dplyr::tibble(
        'npar_proj' = means_by_npar[[npar_proj]],
        'npar_clus' = means_by_npar[[npar_clus]],
        'Evaluation Metric' = means_by_npar$name,
        'Mean Value' = means_by_npar$value,
        'Standard Deviation' = sd_by_npar$value
      )
      colnames(vals_by_npar)[1:2] <- c(npar_proj, npar_clus)
      idcs_keep <- apply(vals_by_npar, 2, function(x) !all(is.na(x)))
      vals_by_npar <- vals_by_npar[idcs_keep]
    }
    vals_by_npar$`Evaluation Metric` <- as.factor(vals_by_npar$`Evaluation Metric`)
    return(vals_by_npar)
  }
}

GetProjectionReference <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL
) {
  rhdf5::h5read(benchmark$h5_path, .h5_slotname(idx.subpipeline, 'Projection', idx.n_param, suffix = 'IsReferenceTo'))
}

GetRMSDPerCluster <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = NULL
) {
  cl <- GetClustering(benchmark, idx.subpipeline, idx.n_param, all_runs = TRUE)
  exprs <- GetExpressionMatrix(benchmark)
  if (benchmark$stability == 'repeat' && !is.null(idx.run)) {
    cl <- cl[[idx.run]]
    return(rmsd_per_cluster(exprs, cl))
  } else if (benchmark$stability %in% c('single', 'repeat')) {
    return(rmsd_per_cluster(exprs, cl))
  } else { # return results for multiple runs
    res <- purrr::map(cl, function(x) rmsd_per_cluster(exprs, x))
    names(res) <- paste0('Run', seq_along(res))
    return(res)
  }
}

GetLabelClusterMatching <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = NULL
) {
  s <- GetClusteringScores(b, idx.subpipeline, idx.n_param)
  bij <- s[[1]]$`Label-Cluster Matching (Bijective)`
  fc <- s[[1]]$`Label-Cluster Matching (Relaxed, Fixed Cluster)`
  fl <- s[[1]]$`Label-Cluster Matching (Relaxed, Fixed Label)`
  if (benchmark$stability == 'repeat' && !is.null(idx.run)) {
    bij <- bij[bij$Run == idx.run, -1]
    fc <- fc[fc$Run == idx.run, -1]
    fl <- fl[fl$Run == idx.run, -1]
  }
  if (ncol(bij) == 2) {
    bij <- data.frame(bij); colnames(bij) <- c('Population', 'Cluster')
    fc <- data.frame(fc); colnames(fc) <- c('Population', 'Cluster')
    fl <- data.frame(fl); colnames(fl) <- c('Population', 'Cluster')
  }
  list(
    'Bijective' = bij,
    'Fixed Cluster' = fc,
    'Fixed Label' = fl
  )
}

GetRMSDPerPopulation <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = NULL,
  match_type = 'bijective'
) {
  rmsd <- GetRMSDPerCluster(benchmark, idx.subpipeline, idx.n_param, idx.run)
  matching <- GetLabelClusterMatching(benchmark, idx.subpipeline, idx.n_param, idx.run)
  
  match_type <- match.arg(match_type, choices = c('bijective', 'fixed_cluster', 'fixed_label', 'real'))
  if (match_type == 'bijective') {
    matching <- matching$Bijective
  } else if (match_type == 'fixed_cluster') {
    matching <- matching$`Fixed Cluster`
  } else if (match_type == 'fixed_label') {
    matching <- matching$`Fixed Label`
  } else if (match_type == 'real') {
    annot <- GetAnnotation(b)
    res <- rmsd_per_cluster(GetExpressionMatrix(b), annot)
    names(res) <- levels(annot)
    if (!is.null(benchmark$unassigned_labels)) {
      res <- res[!names(res) %in% benchmark$unassigned_labels]
    }
    return(as.list(res))
  }
    
  lab_u <- levels(GetAnnotation(benchmark, concatenate = TRUE))
  lab_u <- lab_u[!lab_u %in% benchmark$unassigned_labels]
  if (benchmark$stability == 'single' || !is.null(idx.run)) {
    res <- purrr::map(lab_u, function(pop) {
      idx_clusters <- which(purrr::map_lgl(matching$Population, function(x) isTRUE(x == pop)))
      if (length(idx_clusters) == 0) {
        r <- NA
        names(r) <- 'NA'
        return(r)
      }
      r <- rmsd[as.numeric(matching$Cluster[idx_clusters])]
      names(r) <- idx_clusters
      r
    })
    names(res) <- lab_u
  } else {
    res <- vector(mode = 'list', length = max(matching$Run))
    names(res) <- paste0('Run', seq_along(res))
    for (idx.run in seq_along(res)) {
      m <- matching[matching$Run == idx.run, ]
      idx_clusters <- which(purrr::map_lgl(m$Population, function(x) isTRUE(x == pop)))
      if (length(idx_clusters) == 0) {
        r <- NA
        names(r) <- 'NA'
        return(r)
      }
      res[[idx.run]] <- purrr::map(lab_u, function(pop) {
        r <- rmsd[[idx.run]][as.numeric(m$Cluster[idx_clusters])]
        names(r) <- idx_clusters
        r
      })
      names(res[[idx.run]]) <- lab_u
    }
  }
  res
}

GetPopulationSizes <- function(
  benchmark, as_list = TRUE, include_unassigned = FALSE
) {
  annot <- GetAnnotation(benchmark, concatenate = TRUE)
  if (as_list) {
    res <- as.list(table(annot))
    if (!include_unassigned) {
      res <- res[!names(res) %in% benchmark$unassigned_labels]
    }
  } else {
    res <- as.data.frame(table(annot))
    colnames(res) <- c('Population', 'Size')
    if (!include_unassigned) {
      res <- res[!res$Population %in% benchmark$unassigned_labels, ]
    }
  }
  res
}

GetClusterSizes <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = NULL,
  as_list = TRUE
) {
  cl <- GetClustering(benchmark, idx.subpipeline, idx.n_param, all_runs = TRUE)
  
  if (benchmark$stability == 'repeat' && !is.null(idx.run)) {
    cl <- cl[[idx.run]]
  }
  
  if (is.list(cl)) {
    res <- vector(mode = 'list', length = length(cl))
    names(res) <- paste0('Run', seq_along(res))
    for (idx.run in seq_along(res)) {
      if (as_list) {
        res[[idx.run]] <- as.list(table(cl[[idx.run]]))
      } else {
        res[[idx.run]] <- as.data.frame(table(cl[[idx.run]]))
        colnames(res[[idx.run]]) <- c('Cluster', 'Size')
      }
    }
  } else {
    if (as_list) {
      res <- as.list(table(cl))
    } else {
      res <- as.data.frame(table(cl))
      colnames(res) <- c('Cluster', 'Size')
    }
  }
  
  res
}

GetMatchedClusterSizes <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = NULL,
  match_type = 'bijective'
) {
  sizes <- GetClusterSizes(benchmark, idx.subpipeline, idx.n_param, idx.run, as_list = FALSE)
  matching <- GetLabelClusterMatching(benchmark, idx.subpipeline, idx.n_param, idx.run)
  
  match_type <- match.arg(match_type, choices = c('bijective', 'fixed_cluster', 'fixed_label'))
  if (match_type == 'bijective') {
    matching <- matching$Bijective
  } else if (match_type == 'fixed_cluster') {
    matching <- matching$`Fixed Cluster`
  } else if (match_type == 'fixed_label') {
    matching <- matching$`Fixed Label`
  }
  
  lab_u <- levels(GetAnnotation(benchmark, concatenate = TRUE))
  lab_u <- lab_u[!lab_u %in% benchmark$unassigned_labels]
  if (benchmark$stability == 'single' || !is.null(idx.run)) {
    res <- purrr::map(lab_u, function(pop) {
      idx_clusters <- which(purrr::map_lgl(matching$Population, function(x) isTRUE(x == pop)))
      if (length(idx_clusters) == 0) return(0)
      idcs <- match(as.numeric(matching$Cluster[idx_clusters]), sizes$Cluster)
      r <- sizes$Size[idcs]
      names(r) <- idx_clusters
      r
    })
    names(res) <- lab_u
  } else {
    res <- vector(mode = 'list', length = max(matching$Run))
    names(res) <- paste0('Run', seq_along(res))
    for (idx.run in seq_along(res)) {
      m <- matching[matching$Run == idx.run, ]
      res[[idx.run]] <- purrr::map(lab_u, function(pop) {
        idx_clusters <- which(purrr::map_lgl(m$Population, function(x) isTRUE(x == pop)))
        if (length(idx_clusters) == 0) return(0)
        idcs <- match(as.numeric(m$Cluster[idx_clusters]), sizes[[idx.run]]$Cluster)
        r <- sizes[[idx.run]]$Size[idcs]
        names(r) <- idx_clusters
        r
      })
      names(res[[idx.run]]) <- lab_u
    }
  }
  res
}
