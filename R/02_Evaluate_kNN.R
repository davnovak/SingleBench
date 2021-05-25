
SavekNNMatrix <- function(
  benchmark,
  knn,
  verbose
) {
  if (verbose) .msg('Saving k-NN matrix...')
  .h5writekNNMatrix(benchmark, knn)
  if (verbose) .msg_alt_good(' done\n')
}

SaveDistanceMatrix <- function(
  benchmark,
  knn,
  verbose
) {
  if (verbose) .msg('Saving distance matrix...')
  .h5writeDistanceMatrix(benchmark, knn)
  if (verbose) .msg_alt_good(' done\n')
}


Evaluate_ComputekNNMatrix <- function(
  benchmark,
  verbose
) {
  exprs <- GetExpressionMatrix(benchmark, concatenate = TRUE)
  if (!is.null(benchmark$rel_idcs.knn_features))
    exprs <- exprs[, benchmark$rel_idcs.knn_features]
  if (verbose) .msg('Computing k-NN matrix...')
  
  systime <- NA
  res <- ComputekNNMatrix(exprs, benchmark$knn.k, benchmark$knn.algorithm, out.systime = systime)
  
  if (verbose) .msg_alt_good(' done in ', round(systime['elapsed'], 2), ' seconds\n')
  res
}

#' Compute a \code{k}-nearest-neighbour matrix for expression data
#'
#' Finds \code{k} closest neighbours to each point in high-dimensional space and the distances to these points.
#' You can choose from various exact or approximate algorithms to use: \code{annoy}, \code{kd_tree}, \code{cover_tree}, \code{CR} and \code{brute}.
#' 
#' @param exprs numeric matrix: a coordinate matrix of biological expression data (columns correspond to markers, rows correspond to cells)
#' @param k integer: number of nearest neighbours to find for each point
#' @param method string: \code{k}-NN algorithm to use. Default value is \code{annoy}
#' @param out.systime optional out-variable: if an object is passed as \code{out.systime}, a side-effect of executing this function is that this object will be assigned elapsed time (in seconds) needed to complete the \code{k}-NN search
#'
#' @return list with two slots: \code{Indices} contains a matrix of nearest neighbours to each point (per row) and \code{Distances} contains a matrix of corresponding Euclidean distances
#'
#' @seealso 
#'
#' * **\code{Denoise}**: denoises high-dimensional expression data to drive down unwanted technical variation
#'
#' @export
ComputekNNMatrix <- function(
  exprs,
  k,
  method = 'annoy',
  dist = 'euclidean', # only Euclidean supported now
  out.systime = NULL
) {
  if (method == 'annoy') {
    systime <- system.time(res <- knn_annoy(exprs, k, 300))
  } else {
    systime <- system.time(knn <- FNN::get.knn(exprs, k = k, algorithm = method))
    res <- list('Indices' = knn$nn.index, 'Distances' = knn$nn.dist)
  }
  if (!is.null(out.systime))
    eval.parent(substitute(out.systime <- systime['elapsed']))
  res
}

knn_annoy <- function(exprs, k, ntrees) {
  reticulate::source_python(file.path(system.file(package = 'SingleBench'), 'kNNAnnoy_ivis.py'))
  reticulate::source_python(file.path(system.file(package = 'SingleBench'), 'kNNAnnoy_GenerateMatrix.py'))
  tmp <- tempfile(tmpdir = '.', pattern = 'AnnoyIndex')
  reticulate::py_capture_output(res <- knn_annoy_python(exprs, nrow(exprs), ncol(exprs), tmpfile = tmp, ntrees = as.integer(ntrees), k = as.integer(k), distf = 'euclidean'))
  file.remove(tmp)
  Indices <- res[, 1, ][, -1] + 1
  Distances <- res[, 2, ][, -1]
  list(
    Indices = Indices,
    Distances = Distances
  )
}
  