
wrapper.projection.VAE <- WrapTool(
  name = 'VAE',
  type = 'projection',
  r_packages = c('ruta', 'tensorflow', 'keras'),
  python_modules = c('tensorflow', 'keras'),
  fun.build_model =
    function(
      input, latent_dim, encoder_shape = c(256, 512, 1024), decoder_shape = c(1024, 512, 256), activation = 'relu',
      dropout = 0.1, validation_split = 0.1, epochs = 20
    ) {
      require(ruta)
      if (tensorflow::tf$executing_eagerly()) {
        tensorflow::tf$compat$v1$disable_eager_execution()
      }
      network <- input()
      for (d in as.integer(encoder_shape)) {
        network <- network + dense(units = d, activation = activation)
        if (!is.null(dropout)) {
          network <- network + dropout(rate = dropout)
        }
      }
      network <- network + variational_block(as.integer(latent_dim))
      for (d in as.integer(decoder_shape)) {
        network <- network + dense(units = d, activation = activation)
      }
      network <- network + output('sigmoid')
      learner <- autoencoder_variational(network, loss = 'mean_squared_error')
      
      n <- as.integer(nrow(input))
      
      idcs_train <- sample(1:n, round(n-validation_split*n))
      x_train <- input[idcs_train, ]
      x_test <- input[-idcs_train, ]
      model <- ruta::train(learner = learner, data = x_train, validation_data = x_test, epochs = as.integer(epochs))
      res <- encode(model, data.frame(input))
      list(model, res)
    },
  fun.extract = function(model)
    model[[2]],
  fun.apply_model = function(model, input) {
    require(ruta)
    encode(model, data.frame(input))
  }
)