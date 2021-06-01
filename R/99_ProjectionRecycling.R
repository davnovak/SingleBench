
FindProjectionResultIfAlreadyGenerated <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param,
  n_param
) {
  
  if (is.na(n_param))
    return(NULL)
  
  clone_group_idx <-
    if (length(benchmark$clone_groups) == 0)
      c()
    else
      which(sapply(benchmark$clone_groups, function(x) idx.subpipeline %in% x))
  if (length(clone_group_idx) == 0) {
    if (is.null(idx.n_param) || idx.n_param < 2)
      return(NULL)
    
    npar <- GetNParameterValues(benchmark, idx.subpipeline)$npar_proj
    npar <- npar[1:(idx.n_param - 1)]
    idx_previous <- which(npar == n_param)[1]
    if (length(idx_previous) > 0 && !is.na(idx_previous))
      return(list(
        idx.subpipeline = idx.subpipeline,
        idx.n_param = idx_previous
      ))
    
    return(NULL)
  }
  
  clone_groups <- benchmark$clone_groups[[clone_group_idx]]
  clone_groups <- clone_groups[clone_groups <= idx.subpipeline]
  clone_groups <- sort(clone_groups)
  
  for (subpipeline_idx in clone_groups) {
    n_param_idx <- which(benchmark$n_params[[subpipeline_idx]]$projection == n_param)[1]
    if (!is.na(n_param_idx) && length(n_param_idx) > 0 && ((subpipeline_idx < idx.subpipeline) || (subpipeline_idx == idx.subpipeline && n_param_idx < idx.n_param))) {
      return(list(
        idx.subpipeline = subpipeline_idx,
        idx.n_param = n_param_idx[1]
      ))
    }
  }
  return(NULL)
}