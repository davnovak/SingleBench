.onAttach <- function(
  libname,
  pkgname
) {
  packageStartupMessage('Loading tool wrappers...\n')
  source(file.path(system.file(package = 'SingleBench'), 'extdata', 'load_wrappers.R'))
  if (require(reticulate)) {
    packageStartupMessage('reticulate detected: interfacing with Python should be possible on this machine')
  }
}
