
wrapper.clustering.Depeche <- WrapTool(
  name = 'Depeche',
  type = 'clustering',
  r_packages = c('DepecheR'),
  fun.build_model =
    function(input, n_clusters, fixed_penalty = NULL, penalties = 2^seq(0, 5, by = 0.5), minARIImprovement = 0.01, optimARI = 0.95) {
      if (!is.null(fixed_penalty)) {
        penalties <- fixed_penalty
      }
      suppressMessages(res <- DepecheR::depeche(input,penalties = penalties, minARIImprovement = minARIImprovement, optimARI = optimARI, createOutput = FALSE))
      res
    },
  fun.extract = function(model)
    model$clusterVector
)