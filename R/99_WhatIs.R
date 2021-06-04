
# WhatIs <- function(
#   benchmark, idx.subpipeline, idx.n_param = NULL, idx.run = 1, cluster = NULL, population = NULL
# ) {
#   .WhatIs.ValidityChecks(environment())
#   
#   if (!is.null(cluster)) {
#     cl <- GetClustering(benchmark, idx.subpipeline, idx.n_param, all_runs = TRUE, concatenate = TRUE)
#     if (is.list(cl)) {
#       cl <- cl[[if (is.null(idx.run)) 1 else idx.run]]
#     }
#     idcs <- which(cl == cluster)
#     if (is.null(idcs))
#       stop('Cluster not found')
#   } else {
#     ann <- GetAnnotation(benchmark, concatenate = TRUE)
#     idcs <- which(ann == population)
#     if (is.null(idcs))
#       stop('Population not found')
#   }
#   
#   
#   
# }