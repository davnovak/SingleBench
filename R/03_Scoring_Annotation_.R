
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