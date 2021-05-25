
BijectiveMatch <- function(
  c1, c2, obj, unassigned, scoring_matrix = NULL
) {
  ## Create integer mapping of each level of c1 and c2 with NAs for unassigned's
  
  levels(c1)[levels(c1) %in% unassigned] <- NA
  c2[is.na(c1)] <- NA
  
  l_c1 <- levels(c1)
  l_c2 <- levels(c2)
  
  c1 <- as.integer(c1)
  c2 <- as.integer(c2)
  
  ## Store group names and number of groups
  
  n_c1 <- length(l_c1)
  n_c2 <- length(l_c2)
  
  ## Create matrices for Precision, Recall and F1:
  ## columns correspond to c1 groups and rows correpond to c2 groups
  
  M_Precision <- M_Recall <- M_F1 <-
    matrix(NA, nrow = n_c2, ncol = n_c1, dimnames = list(l_c2, l_c1))
  
  for (idx_c1 in seq_along(l_c1)) {
    for (idx_c2 in seq_along(l_c2)) {
        TruePositive <- sum(c1 == idx_c1 & c2 == idx_c2, na.rm = TRUE)
        Positive     <- sum(c2 == idx_c2, na.rm = TRUE)
        True         <- sum(c1 == idx_c1, na.rm = TRUE)
        
        Precision <- M_Precision[idx_c2, idx_c1] <- if (Positive == 0) 0 else TruePositive / Positive
        Recall    <- M_Recall[idx_c2, idx_c1]    <- if (True == 0) 0 else TruePositive / True
        M_F1[idx_c2, idx_c1]                     <- if (Precision + Recall == 0) 0 else 2 * Precision * Recall / (Precision + Recall)
    }
  }
  
  ## Choose the metric to maximise in group matching and take the matrix with its values
  
  M <- switch(obj, f1 = M_F1, precision = M_Precision, recall = M_Recall)
  
  ## Match c1 groups to c2 groups by solving linear sum assignment problem
  ## (treat the clustering with more groups as reference when calling the LSAP solver)
  
  if (n_c1 <= n_c2) {
    matches <- clue::solve_LSAP(t(M), maximum = TRUE)
    matches <- l_c2[as.numeric(matches)]
    names(matches) <- l_c1
  } else {
    tmp_matches <- clue::solve_LSAP(M, maximum = TRUE)
    tmp_matches <- l_c1[as.numeric(tmp_matches)]
    names(tmp_matches) <- l_c2
    
    matches <- sapply(l_c1, function(l) if (l %in% tmp_matches) names(tmp_matches)[tmp_matches == l] else NA)
    names(matches) <- l_c1
  }
  
  ## Get Precision, Recall, F1 and number of matched data points (in c2) for each matched pair of groups
  
  Precision.PerMatch <-
    mapply(names(matches), matches, FUN = function(idx_c1, idx_c2) ifelse(is.na(matches[as.character(idx_c1)]), 0, M_Precision[as.character(idx_c2), as.character(idx_c1)]))
  Recall.PerMatch <-
    mapply(names(matches), matches, FUN = function(idx_c1, idx_c2) ifelse(is.na(matches[as.character(idx_c1)]), 0, M_Recall[as.character(idx_c2), as.character(idx_c1)]))
  F1.PerMatch <-
    mapply(names(matches), matches, FUN = function(idx_c1, idx_c2) ifelse(is.na(matches[as.character(idx_c1)]), 0, M_F1[as.character(idx_c2), as.character(idx_c1)]))
  NPoints.PerMatch <-
    sapply(matches, FUN = function(idx_c2) sum(c2 == idx_c2, na.rm = TRUE))
  
  names(Precision.PerMatch) <- names(Recall.PerMatch) <- names(F1.PerMatch) <- names(matches)
  
  ## Compute mean stats values across matches, unweighted and weighted by size of group in c1 (reference)
  
  Precision.UnweightedMean <- mean(Precision.PerMatch)
  Recall.UnweightedMean    <- mean(Recall.PerMatch)
  F1.UnweightedMean        <- mean(F1.PerMatch)
  
  ## If scoring matrix is available, return the hierarchical-fit score
  
  HierarchicalFitScore <- NULL
  
  if (!is.null(scoring_matrix)) {
    for (idx in seq_along(l_c2))
      l_c2[idx] <- if (l_c2[idx] %in% matches) names(matches)[which(matches == l_c2[idx])] else NA
    
    c2 <- factor(c2, labels = l_c2)
    c1 <- factor(c1, labels = l_c1)
    
    idcs_na <- which(is.na(c2) | is.na(c1))
    if (length(idcs_na)) {
      c2 <- c2[-idcs_na]
      c1 <- c1[-idcs_na]
    }
    
    fitscores <- sapply(seq_along(c1), function(idx) scoring_matrix[c1[idx], c2[idx]])
    HierarchicalFitScore <- mean(fitscores)
  }
  
  ## Return matching scheme, stats and which objective was used for matching
  
  list(
    Matches                  = matches,
    Objective                = obj,
    Precision.PerMatch       = Precision.PerMatch,
    Recall.PerMatch          = Recall.PerMatch,
    F1.PerMatch              = F1.PerMatch,
    Precision.UnweightedMean = Precision.UnweightedMean,
    Recall.UnweightedMean    = Recall.UnweightedMean,
    F1.UnweightedMean        = F1.UnweightedMean,
    HierarchicalFitScore     = HierarchicalFitScore
  )
}


RelaxedMatch <- function(
  c1, c2, obj, unassigned, reference_is_c2 = FALSE
) {
  ## Create integer mapping of each level of c1 and c2 with NAs for unassigned's
  
  if (reference_is_c2) {
    levels(c2)[levels(c2) %in% unassigned] <- NA
    c1[is.na(c2)] <- NA
  } else {
    levels(c1)[levels(c1) %in% unassigned] <- NA
    c2[is.na(c1)] <- NA
  }
  
  l_c1 <- levels(c1)
  l_c2 <- levels(c2)
  
  c1 <- as.integer(c1)
  c2 <- as.integer(c2)
  
  ## Store group names and number of groups
  
  n_c1 <- length(l_c1)
  n_c2 <- length(l_c2)
  
  ## Iteratively look for optimal match from c2 for each group in c1
  ## Compute Precision, Recall and F1 for each match that is made
  
  matches <- rep(NA, n_c1)
  names(matches) <- l_c1
  
  Scores <- matrix(NA, ncol = 3, nrow = n_c1, dimnames = list(l_c1, c('Precision', 'Recall', 'F1')))
  for (idx_c1 in seq_along(l_c1)) {
    stats <-
      t(sapply(seq_along(l_c2), function(idx_c2) {
        TruePositive <- sum(c1 == idx_c1 & c2 == idx_c2, na.rm = TRUE)
        Positive     <- sum(c2 == idx_c2, na.rm = TRUE)
        True         <- sum(c1 == idx_c1, na.rm = TRUE)
        
        Precision <- if (Positive == 0) 0 else TruePositive / Positive
        Recall    <- if (True == 0) 0 else TruePositive / True
        
        if (reference_is_c2) {
          tmp <- Precision
          Precision <- Recall
          Recall <- tmp
        }
        
        if (Precision == 0)
          Precision <- 1
        
        list(
          precision = Precision,
          recall    = Recall,
          #f1        = if (Precision + Recall == 0) 0 else 2 * Precision * Recall / (Precision + Recall)
          f1        = 2 * Precision * Recall / (Precision + Recall)
        )
      }))
    idx_match <- which.max(stats[, obj])[1]
    matches[idx_c1] <- if (stats[idx_match, obj] > 0) l_c2[idx_match] else NA
    Scores[idx_c1, ] <- unlist(stats[idx_match, ])
  }
  
  ## Get Precision, Recall, F1 and number of matched data points for each matched pair of groups
  
  Precision.PerMatch <- Scores[, 'Precision']
  Recall.PerMatch    <- Scores[, 'Recall']
  F1.PerMatch        <- Scores[, 'F1']
  NPoints.PerMatch   <- sapply(match(matches, l_c2), FUN = function(idx_c2) sum(c2 == idx_c2, na.rm = TRUE))
  
  ## Compute mean stats values across matches, unweighted and weighted by size of group in c1 (reference)
  
  Precision.UnweightedMean <- mean(Precision.PerMatch)
  Recall.UnweightedMean    <- mean(Recall.PerMatch)
  F1.UnweightedMean        <- mean(F1.PerMatch)
  
  ## Return matching scheme, stats and which objective was used for matching
  
  list(
    Matches                  = matches,
    Objective                = obj,
    Precision.PerMatch       = Precision.PerMatch,
    Recall.PerMatch          = Recall.PerMatch,
    F1.PerMatch              = F1.PerMatch,
    Precision.UnweightedMean = Precision.UnweightedMean,
    Recall.UnweightedMean    = Recall.UnweightedMean,
    F1.UnweightedMean        = F1.UnweightedMean
  )
}
