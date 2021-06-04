
#' Apply unsupervised evaluation metrics to score manual annotation of data
#'
#' @param exprs numeric matrix: expression data, with measured markers in columns and events in rows
#' @param annotation factor vector: manual labels per row of \code{exprs}. Can be extracted via \code{GetAnnotation(benchmark, concatenate = TRUE)}
#' @param unassigned_labels character vector: names of label(s) given to cells that are left as unassigned by the annotation strategy. Default value is \code{c()} (empty vector)
#'
#' @seealso
#'
#' @export
ScoreAnnotation <- function(
  exprs,
  annotation,
  unassigned_labels = c()
) {
  
  if (is.list(exprs))
    exprs <- do.call(rbind, exprs)
  
  idcs_assigned <- which(!annotation %in% unassigned_labels)
  
  suppressWarnings(davies_bouldin <- clusterSim::index.DB(x = exprs[idcs_assigned, ], cl = as.integer(annotation[idcs_assigned]))$DB)
  list(
    'Davies-Bouldin Index' = davies_bouldin
  )
}