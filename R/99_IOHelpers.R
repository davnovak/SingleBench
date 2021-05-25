TryWritePermission <- function(path = '.') file.access(path, mode = 2) == 0

IsFormat <- function(
  fname,
  format,
  case_sensitive = FALSE
) {
  s <- strsplit(fname, '.')
  s <- s[[1]][length(s[[1]])]
  if (!case_sensitive) {
    format <- toupper(format)
    s <- toupper(format)
  }
  s == format
}

NameIsLegal <- function(
  fname = '.SingleBatch_LegalFnameTest'
) {
  ## Try write permissions first!
  if (file.exists(fname)) {
    return(TRUE)
  }
  tryCatch({
    file.create(fname)
    file.remove(fname)
  }, error = function(e) {
    stop(paste0('File name "', fname, '" is illegal'))
  })
  TRUE
}

SimplifyPath <- function(
  fpath
) {
  ## Cut cycles from file path
  ## Cut redundant './' from file path
  stroke <- .Platform$file.sep
  ss <- strsplit(fpath, stroke)[[1]]
  if (ss[1] == '') {
    ss <- ss[-1]
  }
  s <- as.list(ss)
  idx <- length(s)
  while(idx > 1) {
    if (s[[idx]] == '..' && s[[idx - 1]] != '..') {
      s[[idx]] <- NULL
      s[[idx - 1]] <- NULL
    } else if (s[[idx]] == '.') {
      s[[idx]] <- NULL
    }
    idx <- min(length(s), idx - 1)
  }
  s <- as.list(unlist(s))
  simple <- do.call(file.path, s)
  if (substr(fpath, 1, 1) == stroke) {
    simple <- paste0(stroke, simple)
  }
  ## Simplify any '//' to '/'
  gsub(paste0(stroke, stroke), stroke, simple, fixed = TRUE)
}

ShortenPath <- function(
  fpath,
  n
) {
  len <- nchar(fpath)
  if (len <= n) {
    return(fpath)
  }
  len_fname <- nchar(basename(fpath))
  if (len_fname >= (n - 4)) {
    d <- len_fname - (n - 4)
    fname <- basename(fpath)
    f <- substr(fname, d, len_fname)
    return(paste0('[..]', f))
  }
  d <- len - n
  s <- strsplit(fpath, '/')[[1]][-1]
  cs <- cumsum(sapply(rev(s)[-1], nchar)) + 1
  idx <- which(cs <= (n - len_fname - 7))
  if (length(idx) == 0) {
    return(paste0('[..]', basename(fpath)))
  }
  cutoff <- length(s) - max(idx)
  params <- c(as.list('[..]'), s[cutoff:(length(s) - 1)], basename(fpath))
  return(do.call(file.path, params))
}