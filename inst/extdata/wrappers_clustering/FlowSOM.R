
wrapper.clustering.FlowSOM <- WrapTool(
  name = 'FlowSOM',
  type = 'clustering',
  r_packages = c('FlowSOM'),
  fun.build_model =
    function(input, n_clusters, grid_width = 10, grid_height = 10, max_metaclusters = 90) {
      som <- FlowSOM::SOM(data = input, xdim = grid_width, ydim = grid_height, silent = TRUE)
      clus <- som$mapping[, 1]
      meta <-
        if (n_clusters == 0)
          FlowSOM::MetaClustering(som$codes, method = 'metaClustering_consensus', max = max_metaclusters)
        else
          FlowSOM::metaClustering_consensus(som$codes, k = n_clusters)
      list(
        som = som,
        clustering = clus,
        metaclustering = meta
      )
    },
  fun.extract = function(model)
    model$metaclustering[model$clustering],
  fun.apply_model = function(model, input)
    model$meta[FlowSOM::MapDataToCodes(codes = model$som$codes, newdata = input)]
)
