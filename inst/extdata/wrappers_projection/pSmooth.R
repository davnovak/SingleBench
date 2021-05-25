
wrapper.projection.pSmooth <- WrapTool(
  name = 'pSmooth',
  type = 'projection',
  r_packages = c(),
  use_knn_graph = TRUE,
  fun.build_model = function(input, knn, k = 50, n_iter = 1, mode = 'Uniform', phi2 = 1, lambda = 1) {
    if (n_iter == 0)
      return(input)
    k <- min(ncol(knn$Distances), k)
    w <-
      if (mode == 'Uniform')
        matrix(1, nrow = nrow(knn$Indices), ncol = ncol(knn$Indices))
      else if (mode == 'Strong')
        exp(-knn$Distances / median(knn$Distances))
      else if (mode == 'Local')
        exp(-knn$Distances^2 / quantile(knn$Distances, 0.95))
      else if (mode == 'Weak')
        t(apply(knn$Distances, 1, function(x) exp(-x^2 / max(x))))
    segs <- rep(1, nrow(input))
    res <- .Call('_SingleBench_pSmoothSingleIter', PACKAGE = 'SingleBench', input, knn$Indices, w, 1, phi2, lambda, segs, k)
    res
  }
)
