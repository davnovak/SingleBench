
rmsd_per_cluster <- function(
  exprs,
  clustering
) {
  
  ## This is a helper function which computes the RMSD per cluster given
  ## expression data and the clustering vector (integer vector). Alternatively,
  ## 'clustering' can also be a logical vector, in which case a single RMSD value
  ## (for the points for which clustering is TRUE) is returned
  
  if (is.list(exprs)) {
    exprs <- do.call(rbind, exprs)
  }
  if (is.list(clustering)) {
    clustering <- unlist(clustering)
  }
  if (is.logical(clustering)) { # 'one-hot' encoding
    res <- sqrt(mean(apply(exprs[clustering, , drop = FALSE], MARGIN = 2, FUN = var)))
  } else {
    res <- purrr::map_dbl(
      .x = sort(unique(clustering)),
      .f = function(cl) sqrt(mean(apply(exprs[clustering == cl, , drop = FALSE], MARGIN = 2, FUN = var)))
    )
    res[is.na(res)] <- 0
  }
  res
}

