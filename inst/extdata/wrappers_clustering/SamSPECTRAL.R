
wrapper.clustering.SamSPECTRAL <- WrapTool(
  name = 'SamSPECTRAL',
  type = 'clustering',
  r_packages = c('SamSPECTRAL'),
  fun.build_model =
    function(input, n_clusters, normal.sigma = 200, separation.factor = 0.4, maximum.number.of.clusters = 30, m = 3000)
      SamSPECTRAL::SamSPECTRAL(data.points = input, normal.sigma = normal.sigma, separation.factor = separation.factor, maximum.number.of.clusters = maximum.number.of.cllusters, m = m, talk = FALSE)
)
