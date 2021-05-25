
#' Extend gating to unassigned cells
#'
#' Uses a \code{k}-NN classifier to label cells that are marked as 'unassigned' in an existing manual annotation vector.
#' A majority vote is taken among nearest neighbours to an unassigned cells, excluding other unassigneds.
#' This procedure can be repeated in multiple consecutive iterations.
#' 
#' @param exprs numeric matrix: a coordinate matrix of biological expression data (columns correspond to markers, rows correspond to cells)
#' @param annotation string or factor vector: vector of labels per each row of \code{exprs}
#' @param unassigned_labels string or string vector: names of population labels to handle as 'unassigned' events (not belonging to a manually annotated population)
#' @param knn_indices numeric matrix: the \code{Indices} slot of a \code{k}-NN object produced by the function \code{ComputekNNMatrix}
#' @param k integer: number of nearest neighbours of each point to use for labelling. Default value is the maximum possible value (\code{k} of the \code{k}-NN matrix)
#' @param n_iter integer: number of labelling iterations. Default value is \code{1}
#' @param verbose logical: whether to display progress messages. Default value is \code{TRUE}
#'
#' @seealso 
#'
#' * **\code{ComputekNNMatrix}**: finds \code{k} nearest neighbours to each point in high-dimensional expression data
#'
#' @export
ExtendGates <- function(
  exprs,
  annotation,
  unassigned_labels,
  knn_indices,
  k = NULL,
  n_iter = 1,
  verbose = TRUE
) {
  if (is.null(k)) k <- ncol(knn_indices)
  res <- annotation
  for (idx_iter in seq_len(n_iter)) {
    if (verbose) { .msg('Gate-extension iteration '); .msg_alt(idx_iter); .msg(' of '); .msg_alt(n_iter, '\n') }
    idcs_unassigned <- which(annotation %in% unassigned_labels)
    count <- 0
    for (idx_label in idcs_unassigned) {
      count <- count + 1
      labels <- annotation[knn_indices[idx_label, seq_len(k)]]
      labels <- labels[!labels %in% unassigned_labels]
      if (length(labels) > 0) {
        plurality <- tail(names(sort(table(labels))), 1)
        res[idx_label] <- plurality
      }
      if (verbose && count %% 10000 == 0) { .msg('-> unassigned label '); .msg_alt(count); .msg(' of '); .msg_alt(length(idcs_unassigned), '\n') }
    }
    annotation <- res
  }
  res
}