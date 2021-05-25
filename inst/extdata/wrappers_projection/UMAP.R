
wrapper.projection.UMAP <- WrapTool(
  name = 'UMAP',
  type = 'projection',
  r_packages = c('uwot'),
  fun.build_model =
    function(input, latent_dim, n_neighbours = 15, scale = FALSE, metric = 'euclidean')
      uwot::umap(input, n_neighbors = n_neighbours, n_components = latent_dim, scale = scale, metric = metric, ret_model = TRUE),
  fun.extract =
    function(model)
      model$embedding,
  fun.apply_model =
    function(model, input)
    uwot::umap_transform(X = input, model = model)
)
