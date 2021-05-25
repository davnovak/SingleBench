
wrapper.clustering.Phenograph <- WrapTool(
  name = 'Phenograph',
  type = 'clustering',
  r_packages = c('Rphenograph', 'igraph'),
  use_knn_graph = TRUE,
  fun.build_model =
    function(input, n_clusters, knn, n_neighbours = 30) {
      nn <- knn$Indices[, 1:n_neighbours]
      links <- Rphenograph:::jaccard_coeff(nn)
      relations <- as.data.frame(links)
      colnames(relations) <- c('from', 'to', 'weight')
      
      g <- igraph::graph.data.frame(relations, directed = FALSE)
      community <- igraph::cluster_louvain(g)
      list(g, community, nrow(input))
    },
  fun.extract = function(model) {
    res <- as.numeric(igraph::membership(model[[2]]))
    res[1:model[[3]]]
  }
)