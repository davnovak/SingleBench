
wrapper.clustering.densitycut <- WrapTool(
  name = 'densitycut',
  type = 'clustering',
  r_packages = c('densitycut'),
  fun.build_model =
    function(input, n_clusters, knn, alpha = 0.9, nu = seq(0, 1, by = 0.05), adjust = TRUE, maxit = 50, eps = 1e-05)
      densitycut::DensityCut(X = input, knn.index = knn$Indices, knn.dist = knn$Distances, alpha = alpha, nu = nu, adjust = adjust, maxit = maxit, eps = eps),
  fun.extract = function(model)
    model$cluster,
  use_knn_graph = TRUE
)