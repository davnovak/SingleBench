
rmsd_per_cluster <- function(
  exprs,
  clustering
) {
  if (is.list(exprs)) {
    exprs <- do.call(rbind, exprs)
  }
  if (is.list(clustering)) {
    clustering <- unlist(clustering)
  }
  res <- purrr::map_dbl(
    .x = sort(unique(clustering)),
    .f = function(cl) sqrt(mean(apply(exprs[clustering == cl, , drop = FALSE], MARGIN = 2, FUN = var)))
  )
  res[is.na(res)] <- 0
  res
}

