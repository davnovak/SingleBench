
wrapper.clustering.flowPeaks <- WrapTool(
  name = 'flowPeaks',
  type = 'clustering',
  r_packages = c('flowPeaks'),
  fun.build_model =
    function(input, tol = 0.1, h0 = 1, h = 1.5) {
      res <- flowPeaks::flowPeaks(x = input, tol = tol, h0 = h0, h = h)
      fp$peaks.cluster
    },
  fun.extract = function(model)
    model
)
