
wrapper.clustering.immunoClust <- WrapTool(
  name = 'immunoClust',
  type = 'clustering',
  r_packages = c('immunoClust'),
  fun.build_model =
    function(input, n_clusters, N = NULL, min.count = -1, max.count = -1, min = NULL, max = NULL, I.buildup = 6, I.final = 4, I.trans = 6, modelName = 'mvt', tol = 1e-5, bias = 0.3, sub.tol = 1e-4, sub.bias = 0.3)
      immunoClust::cell.process(
        fcs = flowCore::flowFrame(input),
        N = N, min.count = min.count, max.count = max.count, min = min, max = max, I.buildup = I.buildup, I.final = I.final, I.trans = I.trans,
        modelName = modelName, tol = tol, bias = bias, sub.tol = sub.tol, sub.bias = sub.bias
      ),
  fun.extract = function(model)
    model@label
)