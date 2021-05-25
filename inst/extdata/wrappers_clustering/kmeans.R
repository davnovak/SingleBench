
wrapper.clustering.kmeans <- WrapTool(
  name = 'kmeans',
  type = 'clustering',
  r_packages = c(),
  fun.build_model =
    function(input, n_clusters)
      stats::kmeans(input, centers = n_clusters),
  fun.extract = function(model)
    model$cluster
)