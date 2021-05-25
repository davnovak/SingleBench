
GetChannelsAndMarkers <- function(
  fcs_paths,
  verbose
) {
  channels <- list()
  markers  <- list()
  for (idx in seq_along(fcs_paths)) {
    file <- fcs_paths[idx]
    ff   <- flowCore::read.FCS(file)
    desc <- na.omit(ff@parameters$desc)
    channels[[idx]] <- purrr::map_chr(names(desc), function(marker_slot) {
      channel_slot <- gsub('.{1}$', 'N', marker_slot)
      channel <- as.character(ff@description[[channel_slot]])
    })
    markers[[idx]] <- as.character(as.vector(desc))
  }
  names(channels) <- names(markers) <- n <- basename(fcs_paths)
  channels.block <- do.call(cbind, channels)
  idcs.c_i <- c() # indices of inconsistent channel values
  markers.block  <- do.call(cbind, channels)
  idcs.m_i <- c() # indices of inconsistent marker values
  for (idx.row in seq_len(nrow(channels.block))) {
    vals  <- channels.block[idx.row, ]
    if (length(unique(vals)) > 1)
      idcs.c_i <- c(idcs.c_i, idx.row)
    vals  <- markers.block[idx.row, ]
    if (length(unique(vals)) > 1)
      idcs.m_i <- c(idcs.m_i, idx.row)
  }
  idcs.exclude <- c()
  n.c_i <- length(idcs.c_i)
  n.m_i <- length(idcs.m_i)
  if ((n.c_i > 0 || n.m_i > 0)) {
    idcs.cm_i <- intersect(idcs.c_i, idcs.c_i)
    idcs.c_i <- idcs.c_i[!idcs.c_i %in% idcs.cm_i]
    idcs.m_i <- idcs.m_i[!idcs.m_i %in% idcs.cm_i]
    for (idx.row in idcs.cm_i) {
      pairs  <- purrr::map_chr(1:ncol(channels.block), function(idx.col) paste(markers.block[idx.row, idx.col], channels.block[idx.row, idx.col], sep = ' --- '))
      names(pairs) <- colnames(channels.block)
      tpairs <- table(pairs)
      if (verbose) {
        .msg_alt_bad('\t\t\t-> inconsistent channel-marker pair in position ', idx.row, '\n')
        for (idx.pair in seq_along(tpairs)) {
          pair <- names(tpairs)[idx.pair]
          .msg_alt_bad('\t\t\t  ---> pair "', pair, '" occurs in ', if (tpairs[idx.pair] > 1) { ' these files:\n' } else { ' this file:\n' })
          .msg_alt(paste0('\t\t\t\t', colnames(channels.block)[which(pairs == pair)], collapse = ',\n'), '\n')
        }
      }
    }
    for (idx.row in idcs.c_i) {
      ch <- channels.block[idx.row, ]
      tch <- table(ch)
      if (verbose) {
        .msg_alt_bad('\t\t\t-> inconsistent channel name in position ', idx.row, '\n')
        for (idx.ch in seq_along(tch)) {
          Ch <- names(tch)[idx.ch]
          .msg_alt_bad('\t\t\t  ---> channel "', Ch, '" occurs in ', if (tch[idx.ch] > 1) { ' these files:\n' } else { ' this file:\n' })
          .msg_alt(paste0('\t\t\t\t', colnames(channels.block)[which(ch) == Ch], collapse = ',\n'), '\n')
        }
      }
    }
    for (idx.row in idcs.m_i) {
      ma <- markers.block[idx.row, ]
      tma <- table(ma)
      if (verbose) {
        .msg_alt_bad('\t\t\t-> inconsistent marker name in position ', idx.row, '\n')
        for (idx.ma in seq_along(tma)) {
          Ma <- names(tma)[idx.ma]
          .msg_alt_bad('\t\t\t  ---> marker "', Ma, '" occurs in ', if (tma[idx.ma] > 1) { ' these files:\n' } else { ' this file:\n' })
          .msg_alt(paste0('\t\t\t\t', colnames(markers.block)[which(ma) == Ma], collapse = ',\n'), '\n')
        }
      }
    }
    idcs.exclude <- c(idcs.m_i, idcs.c_i, idcs.cm_i)
    markers <- markers[[1]]
    markers[idcs.exclude] <- NA
    channels <- channels[[1]]
    channels[idcs.exclude] <- NA
    markers <- na.omit(markers)
    channels <- na.omit(channels)
    if (length(channels) == 0)
      stop('No compatible channel marker pairing found in FCS files')
    if (verbose) {
      .msg_alt_bad('\t\t', length(idcs.exclude), ' out of the total ', nrow(channels.block), ' chanel-marker pairs are not compatible across all input FCS files\n')
      .msg_alt_bad('\t\t', if (length(idcs.exclude) == 1) { 'Position ' } else { 'Positions ' }, paste(idcs.exclude, collapse = ', '), ' will be excluded from analysis\n')
    }
  } else {
    markers <- markers[[1]]
    channels <- channels[[1]]
    if (verbose)
      .msg_alt_good('\t\tAll channel-marker name pairs are compatible across input files\n')
  }
  list(
    markers         = markers,
    channels        = channels,
    pairs           = paste(channels, markers, sep = ' --- '),
    idcs.mismatches = idcs.exclude
  )
}

AlignChannelMarkerPairsNames <- function(
  pairs
) {
  max_len <- max(
      purrr::map_int(pairs, function(pair) {
      s <- strsplit(pair, ' --- ')[[1]]
      nchar(s[1])
    })
  )
  purrr::map_chr(pairs, function(pair) {
    s <- strsplit(pair, ' --- ')[[1]]
    len <- nchar(s[1])
    paste0(paste0(rep(' ', max_len - len), collapse = ''), s[1], ' --- ', s[2])
  })
}