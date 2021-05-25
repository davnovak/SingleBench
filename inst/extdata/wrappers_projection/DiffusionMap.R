
wrapper.projection.DiffusionMap <- WrapTool(
  name = 'DiffusionMap',
  type = 'projection',
  r_packages = c('destiny'),
  fun.build_model = function(
    input, latent_dim, sigma = 'local', k = 150, n_eigs = 20, density_norm = TRUE, distance = 'euclidean'
  ) {
    dm <- destiny::DiffusionMap(data = input, sigma = sigma, k = k, n_eigs = n_eigs, density_norm = density_norm, distance = distance)
    ev <- destiny::eigenvectors(dm)
    list(
      model = dm,
      projection = ev[, 1:min(latent_dim, ncol(ev))]
    )
  },
  fun.extract = function(model)
    model$projection
)
