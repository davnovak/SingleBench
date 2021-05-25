
wrapper.projection.tSNE <- WrapTool(
  name = 'tSNE',
  type = 'projection',
  r_packages = c('Rtsne'),
  fun.build_model =
    function(input, latent_dim, initial_dim = 50, perplexity = 30, theta = 0.5)
      Rtsne::Rtsne(input, dims = latent_dim, check_duplicates = FALSE, initial_dims = initial_dim, perplexity = perplexity, theta = theta),
  fun.extract =
    function(model)
      model$Y
)