
wrapper.projection.ivis <- WrapTool(
  name = 'ivis',
  type = 'projection',
  r_packages = c(),
  python_modules = c('ivis'),
  use_knn_graph = TRUE,
  fun.build_model =
    function(
      input, latent_dim, knn, scale = TRUE, max_k = 150L, distance = 'pn', batch_size = 128L, epochs = 1000L, model = 'szubert',
      supervision_metric = 'mean_squared_error', supervision_weight = 0.5
    ) {
      if (scale == TRUE) input <- scale(input)
      ivis <- reticulate::import('ivis')
      knn <- kNNGTweak(knn, only_indices = TRUE, zero_index = TRUE, zeroth_neighbours = TRUE, new_k = min(ncol(knn$Indices), max_k))
      reticulate::py_capture_output(model <- ivis$Ivis(embedding_dims = as.integer(latent_dim), neighbour_matrix = knn, distance = distance, batch_size = batch_size, epochs = epochs, supervision_metric = supervision_metric, supervision_weight = supervision_weight))
      reticulate::py_capture_output(model <- model$fit(input))
      reticulate::py_capture_output(coordinates <- model$transform(input))
      list(
        model = model,
        coordinates = coordinates
      )
    },
  fun.extract = function(model)
    model$coordinates,
  fun.apply_model = function(model, input) {
    suppressMessages(res <- model$model$transform(input))
    res
  }
)