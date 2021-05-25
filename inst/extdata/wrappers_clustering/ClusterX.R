
wrapper.clustering.ClusterX <- WrapTool(
  name = 'ClusterX',
  type = 'clustering',
  r_packages = c('cytofkit'),
  fun.build_model =
    function(input, expression, n_clusters){
      suppressMessages(cytofkit::cytof_cluster(expression, input, method = 'ClusterX'))
    },
  fun.extract = function(model)
    model,
  use_original_expression_matrix = TRUE
)
