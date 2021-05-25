
wrapper.clustering.cytometree <- WrapTool(
  name = 'cytometree',
  type = 'clustering',
  r_packages = c('cytometree'),
  fun.build_model =
    function(input, n_clusters, minleaf = 1, t = 0.1, force_first_markers = NULL)
      cytometree::CytomeTree(input, minleaf = minleaf, t = t, force_first_markers = force_first_markers),
  fun.extract = function(model)
    model$labels
)