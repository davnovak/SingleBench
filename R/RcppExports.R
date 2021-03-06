# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

pSmoothSingleIter <- function(exprs, idcs, weights, phi1, phi2, lambda, segs, k) {
    .Call('_SingleBench_pSmoothSingleIter', PACKAGE = 'SingleBench', exprs, idcs, weights, phi1, phi2, lambda, segs, k)
}

CollapsedCorankingMatrix <- function(idcs1, idcs2) {
    .Call('_SingleBench_CollapsedCorankingMatrix', PACKAGE = 'SingleBench', idcs1, idcs2)
}

denoise_matrix <- function(coords, kNN_idcs, K, n_iter) {
    .Call('_SingleBench_denoise_matrix', PACKAGE = 'SingleBench', coords, kNN_idcs, K, n_iter)
}

