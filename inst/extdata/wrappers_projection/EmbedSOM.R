
## Not final, latent_dim not fixable

wrapper.projection.EmbedSOM <- WrapTool(
  name = 'EmbedSOM',
  type = 'projection',
  r_packages = c('EmbedSOM'),
  fun.build_model =
    function(input, latent_dim, xdim = 8, ydim = 8, zdim = 8) {
      map <- EmbedSOM::SOM(input, xdim = xdim, ydim = ydim, zdim = zdim)
      EmbedSOM::EmbedSOM(data = input, map = map)
    }
)
