
SeparateIntoSamples <- function(
  obj, benchmark, single = FALSE
) {
  
  ## This function takes the concatenated output of a projection/clustering tool
  ## (transformed data or a vector of cluster indices) and converts it to a list
  ## of separate outputs per input sample (if multiple separate samples were given)
  ## as input
  
  if (benchmark$stability == 'single' || single) {
    if (length(benchmark$row_count) == 1)
      return(obj)
    ranges <- dplyr::tibble(
      From = c(1, benchmark$row_count + 1)[1:length(benchmark$row_count)],
      To = cumsum(benchmark$row_count)
    )
    if (is.matrix(obj))
      return(purrr::map2(ranges$From, ranges$To, function(from, to) obj[from:to, ]))  
    else if (is.vector(obj))
      return(purrr::map2(ranges$From, ranges$To, function(from, to) obj[from:to]))
    
  } else if (benchmark$stability == 'repeat') {
    if (is.list(obj))
      return(purrr::map(obj, function(x) SeparateIntoSamples(x, benchmark, single = TRUE)))
    else
      return(SeparateIntoSamples(obj, benchmark, single = TRUE))
  } else {
    return(obj)
  }
}
  

ResolveFeatureIndices <- function(
  idcs_features, channels, markers
) {
  if (is.character(idcs_features)) {
    return(
      if (all(idcs_features %in% channels)) which(channels %in% idcs_features)
      else if (all(idcs_features %in% markers)) which(markers %in% idcs_features)
      else unique(c(which(channels %in% idcs_features), which(markers %in% idcs_features)))
    )
  } else
    return(idcs_features)
}

GatherInput <- function(
  benchmark,
  input,
  input_class,
  input_features,
  compensation_matrix,
  transform,
  transform_features,
  transform_cofactor,
  input_labels,
  input_marker_types,
  batch_id,
  remove_labels,
  verbose
) {
  
  ## This function gathers and pre-processes input data to a benchmark pipeline.
  ## These data come either from paths to FCS files, a flowSet object or a
  ## SummarizedExperiment object
  
  if (!is.null(input_labels) && !is.list(input_labels))
    input_labels <- list(input_labels)
  if (input_class == 'fcs_paths')
    GatherInput_FCSPaths(benchmark, input, input_class, input_features, compensation_matrix, transform, transform_features, transform_cofactor, input_labels, input_marker_types, batch_id, remove_labels, verbose)
  else if (input_class == 'flowSet')
    GatherInput_flowSet(benchmark, input, input_class, input_features, compensation_matrix, transform, transform_features, transform_cofactor, input_labels, input_marker_types, batch_id, remove_labels, verbose)
  else if (input_class == 'SummarizedExperiment')
    GatherInput_SummarizedExperiment(benchmark, input, input_class, input_features, compensation_matrix, transform, transform_features, transform_cofactor, input_labels, input_marker_types, batch_id, remove_labels, verbose)
  benchmark$column_names    <- colnames(if (is.list(benchmark$exprs)) benchmark$exprs[[1]] else benchmark$exprs)
  cn_na <- is.na(benchmark$column_names)
  if (any(cn_na))
    benchmark$column_names[cn_na] <- paste0('Unititled', seq_len(sum(cn_na)))
  benchmark$row_count       <- if (is.list(benchmark$exprs)) sapply(benchmark$exprs, nrow) else nrow(benchmark$exprs)
  benchmark$column_count    <- if (is.list(benchmark$exprs)) ncol(benchmark$exprs[[1]]) else ncol(benchmark$exprs)
  benchmark$n_input_samples <- if (is.list(benchmark$exprs)) length(benchmark$exprs) else 1
  invisible(benchmark)
}

GatherInput_FCSPaths <- function(
  benchmark,
  input,
  input_class,
  input_features,
  compensation_matrix,
  transform,
  transform_features,
  transform_cofactor,
  input_labels,
  input_marker_types,
  batch_id,
  remove_labels,
  verbose
) {
  
  ## Resolve indices of feaures to use in analysis and to transform
  if (!is.null(input_features))
    input_features <- ResolveFeatureIndices(input_features, colnames(input), flowCore::markernames(input))
  if (!is.null(transform_features))
    transform_features <- ResolveFeatureIndices(transform_features, colnames(input), flowCore::markernames(input))
  else
    transform_features <- input_features
  
  ## Extract channel-marker pairs
  channels_and_markers <- GetChannelsAndMarkers(input, verbose)
  if (length(channels_and_markers$channels) == 0)
    stop('No compatible channels found across input FCS files')
  
  ## Extract and pre-process expression data
  benchmark$exprs <- purrr::map(
    input,
    function(fpath) {
      if (verbose) { .msg('Importing file '); .msg_name(input, '\n') }
      ff <- flowCore::read.FCS(fpath)
      if (!is.null(compensation_matrix)) {
        if (verbose) .msg('-> applying compensation\n')
        ff <- flowCore::compensate(ff, compensation_matrix)
      }
      if (!is.null(transform)) {
        if (verbose) { .msg('-> applying '); .msg_alt(transform); .msg(' transformation'); if (!is.null(transform_cofactor)) { .msg(' (cofactor = '); .msg_alt(transform_cofactor); .msg(')') }; .msg('\n') }
        if (is.null(transform_features)) transform_features <- seq_len(ncol(ff))
        if (transform == 'asinh') ff@exprs[, transform_features] <- asinh(ff@exprs[, transform_features] / transform_cofactor)
        else if (transform == 'estimateLogicle') {
          cols <- as.vector(colnames(ff@exprs)[transform_features])
          tf <- flowCore::estimateLogicle(ff, cols)
          ff <- flowCore::transform(ff, tf)
        }
      }
      if (is.null(input_features))
        input_features <- channels_and_markers$channels
      cn <- colnames(ff@exprs)
      cn[cn %in% channels_and_markers$channels] <- channels_and_markers$markers
      colnames(ff@exprs) <- cn
      ff@exprs[, input_features]
    }
  )
  
  ## Extract labels per event
  benchmark$annotation <- input_labels
  
  ## Filter events by labels
  if (!is.null(remove_labels)) {
    if (verbose) {
      .msg(paste0('-> removing all events with ', if (length(remove_labels) == 1) 'label ' else 'labels: '))
      .msg_alt(paste0(remove_labels, collapse = ', '), '\n')
    }
    if (is.list(benchmark$annotation)) {
      idcs_remove <- lapply(benchmark$annotation, function(x) x %in% remove_labels)
      for (idx.input in seq_along(benchmark$input)) {
        benchmark$exprs[[idx.input]] <- benchmark$exprs[[idx.input]][!idcs_remove[[idx.input]]]
        benchmark$annotation[[idx.input]] <- benchmark$annotation[[idx.input]][!idcs_remove[[idx.input]]]
        levels(benchmark$annotation[[idx.input]])[levels(benchmark$annotation[[idx.input]]) %in% remove_labels] <- NA
      }
    } else {
      idcs_remove <- benchmark$annotation %in% remove_labels
      benchmark$exprs[[1]] <- benchmark$exprs[[1]][!idcs_remove[[idx.input]]]
      levels(benchmark$annotation)[levels(benchmark$annotation) %in% remove_labels] <- NA
    }
  }
  
  ## Get relative transform and k-NN feature indices
  if (!is.null(benchmark$rel_idcs.transform_features))
    benchmark$rel_idcs.transform_features <- sapply(benchmark$rel_idcs.transform_features, function(x) which(input_features == x))
  if (!is.null(benchmark$rel_idcs.knn_features))
    benchmark$rel_idcs.knn_features <- sapply(benchmark$rel_idcs.knn_features, function(x) which(input_features == x))
  else
    benchmark$rel_idcs.knn_features <- benchmark$rel_idcs.transform_features
  
  ## Use single expression matrix and label vector if single input file is given
  if (is.list(benchmark$exprs) && length(benchmark$exprs) == 1) {
    benchmark$exprs <- benchmark$exprs[[1]]
    if (is.list(benchmark$annotation))
      benchmark$annotation <- benchmark$annotation[[1]]
  }
  
  benchmark$channels <- channels_and_markers$channels
  benchmark$markers <- channels_and_markers$markers
  invisible(benchmark)
}

GatherInput_flowSet <- function(
  benchmark,
  input,
  input_class,
  input_features,
  compensation_matrix,
  transform,
  transform_features,
  transform_cofactor,
  input_labels,
  input_marker_types,
  batch_id,
  remove_labels,
  verbose
) {
  
  ## Resolve indices of feaures to use in analysis and to transform
  if (!is.null(input_features))
    input_features <- ResolveFeatureIndices(input_features, colnames(input), flowCore::markernames(input))
  if (!is.null(transform_features))
    transform_features <- ResolveFeatureIndices(transform_features, colnames(input), flowCore::markernames(input))
  else
    transform_features <- input_features
  
  ## Extract channels
  channel_names <- colnames(as.list(input@frames)[[1]]@exprs)
  benchmark$channels <- channel_names
  helper_idx <- names(input@frames)[1]
  
  ## Extract labels per event
  if (is.null(input_labels)) {
    annotation_levels <- input@frames[[helper_idx]]@description$POPULATION_INFO[, 2]
    benchmark$annotation <- lapply(input@frames, function(x) {
      labs <- as.factor(x@exprs[, 'population_id'])
      levels(labs) <- annotation_levels[unique(labs)]
      labs
    })
  } else {
    benchmark$annotation <- input_labels
  }
  
  ## Extract and pre-process expression data
  data <- as.list(input@frames)
  if (!is.null(compensation_matrix)) {
    if (verbose) .msg('-> applying compensation\n')
    data <- lapply(data, function(x) flowCore::compensate(x, compensation_matrix))
  }
  if (is.null(input_features))
    input_features <- seq_len(ncol(data[[1]]))
  
  benchmark$exprs <- lapply(data, function(x) {
    out <- x@exprs
    
    cn <- colnames(out)
    mn <- flowCore::markernames(x)
    cn[cn %in% names(mn)] <- mn[cn %in% names(mn)]
    
    colnames(out) <- cn
    benchmark$markers <- colnames(out)
    out <- out[, input_features]
    
    if (!is.null(transform)) {
      if (verbose) {
        .msg('-> applying '); .msg_alt(transform); .msg(' transformation to flowFrame')
        if (!is.null(transform_cofactor)) { .msg(' (cofactor = '); .msg_alt(transform_cofactor); .msg(')') };  .msg('\n')
      }
      if (transform == 'asinh') {
        out[, transform_features] <- asinh(out[, transform_features] / transform_cofactor)
      }
      
      else if (transform == 'estimateLogicle')
        out <- flowCore::transform(x, flowCore::estimateLogicle(x, transform_features))
    }
    out
  })
  rm(data)
  
  ## Get relative transform and k-NN feature indices
  if (!is.null(benchmark$rel_idcs.transform_features))
    benchmark$rel_idcs.transform_features <- sapply(benchmark$rel_idcs.transform_features, function(x) which(input_features == x))
  if (!is.null(benchmark$rel_idcs.knn_features))
    benchmark$rel_idcs.knn_features <- sapply(benchmark$rel_idcs.knn_features, function(x) which(input_features == x))
  else
    benchmark$rel_idcs.knn_features <- benchmark$rel_idcs.transform_features
  
  if (!is.null(compensation_matrix) || (!is.null(transform) && transform == 'estimateLogicle'))
    benchmark$exprs <- ff@exprs
  
  
  ## Get relative transform and k-NN feature indices
  if (!is.null(benchmark$rel_idcs.transform_features))
    benchmark$rel_idcs.transform_features <- sapply(benchmark$rel_idcs.transform_features, function(x) which(input_features == x))
  if (!is.null(benchmark$rel_idcs.knn_features))
    benchmark$rel_idcs.knn_features <- sapply(benchmark$rel_idcs.knn_features, function(x) which(input_features == x))
  else
    benchmark$rel_idcs.knn_features <- benchmark$rel_idcs.transform_features
  
  ## Filter events by labels
  if (!is.null(remove_labels)) {
    if (verbose) {
      .msg('-> removing all events with ', if (length(remove_labels) == 1) 'label' else 'labels: ')
      .msg_alt(paste0(remove_labels, collapse = ', '), '\n')
    }
    if (is.list(benchmark$annotation)) {
      idcs_remove <- lapply(benchmark$annotation, function(x) x %in% remove_labels)
      for (idx.input in seq_along(benchmark$input)) {
        benchmark$exprs[[idx.input]] <- benchmark$exprs[[idx.input]][!idcs_remove[[idx.input]]]
        benchmark$annotation[[idx.input]] <- benchmark$annotation[[idx.input]][!idcs_remove[[idx.input]]]
        levels(benchmark$annotation[[idx.input]])[levels(benchmark$annotation[[idx.input]]) %in% remove_labels] <- NA
      }
    } else {
      idcs_remove <- benchmark$annotation %in% remove_labels
      benchmark$exprs[[1]] <- benchmark$exprs[[1]][!idcs_remove[[idx.input]]]
      levels(benchmark$annotation)[levels(benchmark$annotation) %in% remove_labels] <- NA
    }
  }
  
  ## Use single expression matrix and label vector if single input file is given
  if (is.list(benchmark$exprs) && length(benchmark$exprs) == 1) {
    benchmark$exprs <- benchmark$exprs[[1]]
    if (is.list(benchmark$annotation)) benchmark$annotation <- benchmark$annotation[[1]]
  }
  
  invisible(benchmark)
}

GatherInput_SummarizedExperiment <- function(
  benchmark,
  input,
  input_class,
  input_features,
  compensation_matrix,
  transform,
  transform_features,
  transform_cofactor,
  input_labels,
  input_marker_types,
  batch_id,
  remove_labels,
  verbose
) {
  
  ## Resolve indices of feaures to use in analysis and to transform
  if (!is.null(input_marker_types)) {
    if (!is.null(SummarizedExperiment::colData(input))) {
      types_to_keep <- as.character(unique(SummarizedExperiment::colData(input)$marker_class[SummarizedExperiment::colData(input)$marker_class %in% input_marker_types]))
      if (length(types_to_keep) == 0)
        stop('Requested SummarizedExperiment marker type(s) for filtering not found in data')
      if (verbose) {
        .msg(paste0('-> filtering markers by type', if (length(types_to_keep == 1)) ': ' else 's: '))
        .msg_alt(paste0(types_to_keep, collapse = ', '), '\n')
      }
      input_features <-
        SummarizedExperiment::colData(input)$marker_name[SummarizedExperiment::colData(input)$marker_class %in% types_to_keep]
    }
  }
  benchmark$channels <- as.character(SummarizedExperiment::colData(input)$channel_name)
  benchmark$markers <- as.character(SummarizedExperiment::colData(input)$marker_name)
  if (!is.null(input_features))
    input_features <- ResolveFeatureIndices(input_features, input$channel_name, input$marker_name)
  if (!is.null(transform_features))
    transform_features <- ResolveFeatureIndices(transform_features, input$channel_name, input$marker_class)
  else
    transform_features <- input_features
  
  ## Extract labels per event
  if (is.null(input_labels)) {
    rd <- SummarizedExperiment::rowData(input)
    if (!is.null(colnames(rd)) && ('population_id' %in% colnames(rd)))
      benchmark$annotation <- rd$population_id
    else
      benchmark$annotation <- SummarizedExperiment::rowData(input)[, 1]
  }
  
  ## Extract and pre-process expression data
  exprs <- SummarizedExperiment::assays(input)[[1]]
  colnames(exprs) <- input$channel_name
  if (!is.null(compensation_matrix) || (!is.null(transform) && transform == 'estimateLogicle'))
    ff <- flowCore::flowFrame(exprs)
  if (!is.null(compensation_matrix)) {
    if (verbose) .msg('-> applying compensation\n')
    ff <- flowCore::compensate(ff, compensation_matrix)
  }
  if (!is.null(transform)) {
    if (verbose) {
      .msg('-> applying '); .msg_alt(transform); .msg(' transformation')
      if (!is.null(transform_cofactor)) { .msg(' (cofactor = '); .msg_alt(transform_cofactor); .msg(')') }; .msg('\n')
    }
    if (transform == 'asinh')
      exprs[, transform_features] <- asinh(exprs[, transform_features] / transform_cofactor)
    else if (transform == 'estimateLogicle')
      ff <- flowCore::transform(ff, flowCore::estimateLogicle(ff, transform_features))
  }
  
  ## Get relative transform and k-NN feature indices
  if (!is.null(benchmark$rel_idcs.transform_features))
    benchmark$rel_idcs.transform_features <- sapply(benchmark$rel_idcs.transform_features, function(x) which(input_features == x))
  if (!is.null(benchmark$rel_idcs.knn_features))
    benchmark$rel_idcs.knn_features <- sapply(benchmark$rel_idcs.knn_features, function(x) which(input_features == x))
  else
    benchmark$rel_idcs.knn_features <- benchmark$rel_idcs.transform_features
  
  if (!is.null(compensation_matrix) || (!is.null(transform) && transform == 'estimateLogicle'))
    benchmark$exprs <- ff@exprs
  
  benchmark$exprs <- exprs
  rm(exprs)
  colnames(benchmark$exprs) <- input$marker_name
  
  ## Split into batches by row-data
  if (!is.null(batch_id)) {
    batch_idx_per_row <- SummarizedExperiment::rowData(input)[[batch_id]]
    if (!is.null(batch_idx_per_row) && length(unique(batch_idx_per_row)) > 1) {
      if (verbose) { .msg('-> dividing into distinct samples by '); .msg_alt(batch_id, '\n') }
      benchmark$exprs <- lapply(levels(batch_idx_per_row), function(x) benchmark$exprs[batch_idx_per_row == x, ])
      benchmark$annotation <- lapply(levels(batch_idx_per_row), function(x) benchmark$annotation[batch_idx_per_row == x])
    }
  }
  
  ## Filter features by type or indices
  if (!is.null(input_features)) {
    if (is.list(benchmark$exprs))
      benchmark$exprs <- lapply(benchmark$exprs, function(x) x[, input_features])
    else
      benchmark$exprs <- benchmark$exprs[, input_features]
  }
  
  ## Filter events by labels
  if (!is.null(remove_labels)) {
    if (verbose) { .msg('-> removing all events with ', if (length(remove_labels) == 1) 'label ' else 'labels: '); .msg_alt(paste0(remove_labels, collapse = ', '), '\n') }
    if (is.list(benchmark$annotation)) {
      idcs_remove <- lapply(benchmark$annotation, function(x) x %in% remove_labels)
      for (idx.input in seq_along(benchmark$input)) {
        benchmark$exprs[[idx.input]] <- benchmark$exprs[[idx.input]][!idcs_remove[[idx.input]]]
        benchmark$annotation[[idx.input]] <- benchmark$annotation[[idx.input]][!idcs_remove[[idx.input]]]
        levels(benchmark$annotation[[idx.input]])[levels(benchmark$annotation[[idx.input]]) %in% remove_labels] <- NA
      }
    } else {
      idcs_remove <- benchmark$annotation %in% remove_labels
      benchmark$exprs <- benchmark$exprs[!idcs_remove, ]
      benchmark$annotation <- benchmark$annotation[!idcs_remove]
      levels(benchmark$annotation)[levels(benchmark$annotation) %in% remove_labels] <- NA
    }
  }
  
  ## Use single expression matrix and label vector if single input file is given
  if (is.list(benchmark$exprs) && length(benchmark$exprs) == 1) {
    benchmark$exprs <- benchmark$exprs[[1]]
    if (is.list(benchmark$annotation)) benchmark$annotation <- benchmark$annotation[[1]]
  }
  
  invisible(benchmark)
}
