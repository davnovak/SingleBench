
wrapper.projection.PCA <- WrapTool(
  name = 'PCA',
  type = 'projection',
  r_packages = c('stats'),
  fun.build_model = function(
    input, latent_dim, retx = TRUE, center = TRUE, scale. = FALSE
  ) {
    model <- stats::prcomp(x = input, retx = retx, center = center, scale. = scale.)
    list(model, stats::predict(model, input)[, 1:latent_dim])
  },
  fun.extract = function(model)
    model[[2]]
)
