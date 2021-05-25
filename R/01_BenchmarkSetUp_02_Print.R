
PrintNParams <- function(n_params, offset = 1) {

  npars_proj <- n_params$projection
  npars_clus <- n_params$clustering
  
  if (!is.null(npars_proj) || !is.null(npars_clus)) {
    
    offset <- strrep('\t', offset)
    
    if (is.null(npars_clus)) {
      .msg_alt(offset, 'projection n-parameter values: ')
      for (idx_val in seq_along(npars_proj)) {
        .msg_val(abs(npars_proj[idx_val]))
        if (idx_val < length(npars_proj))
          .msg_alt(', ')
      }
      .msg('\n')
    } else if (is.null(npars_proj)) {
      .msg_alt(offset, 'clustering n-parameter values: ')
      for (idx_val in seq_along(npars_clus)) {
        .msg_val(abs(npars_clus[idx_val]))
        if (idx_val < length(npars_clus))
          .msg_alt(', ')
      }
      .msg('\n')
    } else {
      .msg_alt(offset, 'n-parameter value pairs: ')
      for (idx_val in seq_along(npars_proj)) {
        .msg_val(abs(npars_proj[idx_val]), '->', abs(npars_clus[idx_val]))
        if (idx_val < length(npars_clus))
          .msg_alt(', ')
      }
      .msg('\n')
    }
  }
}

print.Benchmark <- function(x, ...) {
  .msg('Benchmark object '); .msg_name(x$name, '\n')
  
  .msg_alt(length(x$n_input_samples)); .msg(' input ', if (length(x$n_input_samples) == 1) 'sample' else 'samples', '\n')
  .msg_alt(length(x$column_names)); .msg(if (length(x$column_names) == 1) ' feature:\n' else ' features:\n')
  .msg_alt('\t', paste0(x$column_names, collapse = ', '), '\n')
  .msg(if (length(x$row_count) == 1) 'Row count: ' else 'Row counts: ')
  .msg_alt(paste0(x$row_count, collapse = ', '), '\n')
  
  if (x$compute_knn) {
    if (!x$knn_available) {
      .msg('-> k-NN matrix will be computed using algorithm '); .msg_name(x$knn.algorithm); .msg(' with k='); .msg_alt(x$knn.k); .msg(' and dist='); .msg_alt(x$knn.distance, '\n')
    } else {
      .msg('-> pre-computed k-NN matrix with k='); .msg_alt(x$knn.k); .msg(' will be used\n')
    }
  }
  if (!x$executable) {
    .msg_alt_bad('No subpipelines set up'); .msg(' (you can extract processed input from this object, but not evaluate it)\n')
  } else {
    .msg_alt(x$n_subpipelines); .msg(if (x$n_subpipelines > 1) ' subpipelines' else ' sub-pipeline', ' set up:\n')
    
    for (idx_subpipeline in seq_along(x$subpipelines)) {
      print(x$subpipelines[[idx_subpipeline]], simple = TRUE, subpipelines_list = x$subpipelines, offset = 1, bullet_number = idx_subpipeline)
      PrintNParams(x$n_params[[idx_subpipeline]], offset = 2)
    }

    if (x$hierarchy)
      .msg_alt('This pipeline is set to use a hierarchical penalty system for clustering evaluation\n')
    if (isTRUE(x$uses_python))
      .msg_alt('This pipeline uses Python via reticulate\n')
    if (isTRUE(x$evaluated_previously))
      .msg_alt_good('Benchmark pipeline was evaluated previously\n')
    if (isTRUE(x$layout_available))
      .msg_alt_good('2-dimensional layout of input data is available\n')
  }
}
