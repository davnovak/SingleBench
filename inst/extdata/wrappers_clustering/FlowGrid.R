
wrapper.clustering.FlowGrid <- WrapTool(
  name = 'FlowGrid',
  type = 'clustering',
  r_packages = c(),
  python_modules = c('FlowGrid'),
  fun.build_model =
    function(input, n_clusters, bin_n = 14, eps = 1.5, MinDenB = 3, MinDenC = 40) {
       fg <- reticulate::import('FlowGrid')
       model <- fg$FlowGrid(input, bin_n = as.integer(bin_n), eps = as.numeric(eps), MinDenB = as.integer(MinDenB), MinDenC = as.integer(MinDenC))
       cl <- model$clustering()
       cl[cl == -1] <- max(max(cl) + 1, 1)
       cl
    }
)