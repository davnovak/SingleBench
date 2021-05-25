
wrapper.clustering.flowMeans <- WrapTool(
  name = 'flowMeans',
  type = 'clustering',
  r_packages = c('flowMeans', 'parallel'),
  fun.build_model =
    function(
      input, n_clusters, MaxN = NA, iter.max = 50, Mahalanobis = TRUE, Standardize = FALSE, Update = 'Mahalanobis', addNoise = TRUE
    )
      flowMeans::flowMeans(
        input, NumC = if (n_clusters == 0) NA else n_clusters, MaxN = MaxN, iter.max = iter.max, Mahalanobis = Mahalanobis,
        Standardize = Standardize, Update = Update, addNoise = addNoise
      ),
  fun.extract = function(model)
    model@Label
)