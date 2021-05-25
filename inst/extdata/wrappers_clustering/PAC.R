
wrapper.clustering.PAC <- WrapTool(
  name = 'PAC',
  type = 'clustering',
  r_packages = c('PAC'),
  fun.build_model =
    function(input, n_clusters, maxlevel = 40, method = 'dsp', max.iter = 50)
      PAC::PAC(data = input, K = n_clusters, maxlevel = maxlevel, method = method, max.iter = max.iter)
)
