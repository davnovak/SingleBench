paths_projection <- dir(file.path(system.file(package = 'SingleBench'), 'extdata', 'wrappers_projection'), full.names = TRUE)
paths_clustering <- dir(file.path(system.file(package = 'SingleBench'), 'extdata', 'wrappers_clustering'), full.names = TRUE)

purrr::walk(paths_projection, source)
purrr::walk(paths_clustering, source)
