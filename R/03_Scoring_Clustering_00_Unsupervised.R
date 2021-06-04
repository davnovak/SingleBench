
ari <- function(c0, c1) {
  
  ## This function computes the Adjusted Rand Index (Rand index)
  ## corrected for change, computed from the contingency table of
  ## two clusterings/label assignments
  
  counts <- table(c0, c1)
  if (nrow(counts) == 1 && ncol(counts) == 1)
    return(1)
  a <- sum(choose(counts, 2))
  b <- sum(choose(rowSums(counts), 2)) - a
  c <- sum(choose(colSums(counts), 2)) - a
  d <- choose(sum(counts), 2) - a - b - c
  (a-(a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
}
