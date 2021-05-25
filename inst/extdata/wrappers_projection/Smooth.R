
wrapper.projection.smooth <- WrapTool(
  name = 'smooth',
  type = 'projection',
  r_packages = c(),
  use_knn_graph = TRUE,
  fun.build_model = function(input, knn, k = 50, n_iter = 1) {
    if (n_iter == 0)
      return(input)
    .Call('_SingleBench_denoise_matrix', PACKAGE = 'SingleBench', input, knn$Indices, k, n_iter)
  }
)
