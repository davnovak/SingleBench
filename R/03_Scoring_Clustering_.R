
ScoreClustering <- function(
  exprs, annotation, result, stability, bootstrap_idcs, unassigned_labels
) {
  
  if (is.list(exprs))
    exprs <- do.call(rbind, exprs)
  
  idcs_assigned <- which(!annotation %in% unassigned_labels)
  
  if (stability == 'single') {
    ## Unsupervised evaluation metrics
    suppressWarnings(davies_bouldin <- clusterSim::index.DB(x = exprs[idcs_assigned, ], cl = result$ClusteringVector[idcs_assigned])$DB)
    
    ## Supervised evaluation metrics
    adjusted_rand_index <- ari(annotation[idcs_assigned], result$ClusteringVector[idcs_assigned])
    nmi <- aricode::NMI(c1 = annotation[idcs_assigned], c2 = result$ClusteringVector[idcs_assigned])
    cluster_matching <- Jacc(c1 = annotation, c2 = result$ClusteringVector, obj = 'f1', unassigned = unassigned_labels, generate_plot = FALSE, verbose = FALSE)
    
    matches_bijective <- data.frame(matrix(c(names(cluster_matching$Results.Bijective$Matches), cluster_matching$Results.Bijective$Matches), ncol = 2))
    colnames(matches_bijective) <- c('Population', 'Cluster')
    
    matches_argmaxf1_cluster <- data.frame(matrix(c(names(cluster_matching$Results.FixedC1$Matches), cluster_matching$Results.FixedC1$Matches), ncol = 2))
    colnames(matches_argmaxf1_cluster) <- c('Population', 'Cluster')
    
    matches_argmaxf1_gate <- data.frame(matrix(c(names(cluster_matching$Results.FixedC2$Matches), cluster_matching$Results.FixedC2$Matches), ncol = 2))
    colnames(matches_argmaxf1_gate) <- c('Cluster', 'Population')
    
    scores <- list(
      'Stability Analysis'                                    = 'single',
      
      'Label-Cluster Matching (Bijective)'                    = matches_bijective,
      'Mean F1 Across Matches (Bijective)'                    = cluster_matching$Results.Bijective$F1.UnweightedMean,
      'Mean Precision Across Matches (Bijective)'             = cluster_matching$Results.Bijective$Precision.UnweightedMean,
      'Mean Recall Across Matches (Bijective)'                = cluster_matching$Results.Bijective$Recall.UnweightedMean,
      'F1 Per Match (Bijective)'                              = cluster_matching$Results.Bijective$F1.PerMatch,
      'Precision Per Match (Bijective)'                       = cluster_matching$Results.Bijective$Precision.PerMatch,
      'Recall Per Match (Bijective)'                          = cluster_matching$Results.Bijective$Recall.PerMatch,
        
      'Label-Cluster Matching (Relaxed, Fixed Label)'          = matches_argmaxf1_cluster,
      'Mean F1 Across Matches (Relaxed, Fixed Label)'          = cluster_matching$Results.FixedC1$F1.UnweightedMean,
      'Mean Precision Across Matches (Relaxed, Fixed Label)'   = cluster_matching$Results.FixedC1$Precision.UnweightedMean,
      'Mean Recall Across Matches (Relaxed, Fixed Label)'      = cluster_matching$Results.FixedC1$Recall.UnweightedMean,
      'F1 Per Match (Relaxed, Fixed Label)'                    = cluster_matching$Results.FixedC1$F1.PerMatch,
      'Precision Per Match (Relaxed, Fixed Label)'             = cluster_matching$Results.FixedC1$Precision.PerMatch,
      'Recall Per Match (Relaxed, Fixed Label)'                = cluster_matching$Results.FixedC1$Recall.PerMatch,
        
      'Label-Cluster Matching (Relaxed, Fixed Cluster)'        = matches_argmaxf1_gate,
      'Mean F1 Across Matches (Relaxed, Fixed Cluster)'        = cluster_matching$Results.FixedC2$F1.UnweightedMean,
      'Mean Precision Across Matches (Relaxed, Fixed Cluster)' = cluster_matching$Results.FixedC2$Precision.UnweightedMean,
      'Mean Recall Across Matches (Relaxed, Fixed Cluster)'    = cluster_matching$Results.FixedC2$Recall.UnweightedMean,
      'F1 Per Match (Relaxed, Fixed Cluster)'                  = cluster_matching$Results.FixedC2$F1.PerMatch,
      'Precision Per Match (Relaxed, Fixed Cluster)'           = cluster_matching$Results.FixedC2$Precision.PerMatch,
      'Recall Per Match (Relaxed, Fixed Cluster)'              = cluster_matching$Results.FixedC2$Recall.PerMatch,
    
      'Davies-Bouldin Index'                                   = davies_bouldin,
      
      'Adjusted Rand Index'                                    = adjusted_rand_index#,
      #'Normalised Mutual Information'                          = nmi
    )
  } else if (stability == 'repeat') {
    ## Unsupervised evaluation metrics
    suppressWarnings(davies_bouldin <- sapply(result$ClusteringVector, function(cv) clusterSim::index.DB(x = exprs[idcs_assigned, ], cl = cv[idcs_assigned])$DB))
    
    ## Supervised evaluation metrics
    adjusted_rand_index <- sapply(result$ClusteringVector, function(cv) ari(annotation[idcs_assigned], cv[idcs_assigned]))
    nmi <- sapply(result$ClusteringVector, function(cv) aricode::NMI(c1 = annotation[idcs_assigned], c2 = cv[idcs_assigned]))
      
    cluster_matching <- purrr::map(result$ClusteringVector, function(cv) Jacc(c1 = annotation, c2 = cv, obj = 'f1', unassigned = unassigned_labels, generate_plot = FALSE, verbose = FALSE))
    
    aggres <- function(data, family) { # aggregate results
      d <- purrr::map(data, function(x) x[[family]])
      subfamilies <- names(d[[1]])
      res <- list()
      for (x in subfamilies) {
        res[[x]] <- purrr::map(d, function(y) y[[x]])
        names(res[[x]]) <- paste0('Run', 1:length(res[[x]]))
      }
      res
    }
    
    cluster_matching.bijective <- aggres(cluster_matching, 'Results.Bijective')
    cluster_matching.argmax_f1_cluster <- aggres(cluster_matching, 'Results.FixedC1')
    cluster_matching.argmax_f1_gate <- aggres(cluster_matching, 'Results.FixedC2')
    
    n_runs <- length(cluster_matching.bijective$Matches)
    matches_bijective <- purrr::map(1:n_runs, function(idx_run) cbind(names(cluster_matching.bijective$Matches[[idx_run]]), cluster_matching.bijective$Matches[[idx_run]]))
    n_pops <- purrr::map_int(matches_bijective, nrow)
    matches_bijective <- data.frame(cbind(rep(1:n_runs, times = n_pops), do.call(rbind, matches_bijective)))
    colnames(matches_bijective) <- c('Run', 'Population', 'Cluster')
    rownames(matches_bijective) <- NULL
    
    matches_argmax_f1_cluster <- purrr::map(1:n_runs, function(idx_run) cbind(names(cluster_matching.argmax_f1_cluster$Matches[[idx_run]]), cluster_matching.argmax_f1_cluster$Matches[[idx_run]]))
    n_pops <- purrr::map_int(matches_argmax_f1_cluster, nrow)
    matches_argmax_f1_cluster <- data.frame(cbind(rep(1:n_runs, times = n_pops), do.call(rbind, matches_argmax_f1_cluster)))
    colnames(matches_argmax_f1_cluster) <- c('Run', 'Population', 'Cluster')
    rownames(matches_argmax_f1_cluster) <- NULL
    
    matches_argmax_f1_gate <- purrr::map(1:n_runs, function(idx_run) cbind(cluster_matching.argmax_f1_gate$Matches[[idx_run]], names(cluster_matching.argmax_f1_gate$Matches[[idx_run]])))
    n_pops <- purrr::map_int(matches_argmax_f1_gate, nrow)
    matches_argmax_f1_gate <- data.frame(cbind(rep(1:n_runs, times = n_pops), do.call(rbind, matches_argmax_f1_gate)))
    colnames(matches_argmax_f1_gate) <- c('Run', 'Population', 'Cluster')
    rownames(matches_argmax_f1_gate) <- NULL
    
    scores <- list(
      'Stability Analysis'                                               = 'repeat',
        
      'Label-Cluster Matching (Bijective)'                               = matches_bijective,
      'Mean F1 Across Matches (Bijective)'                               = as.vector(cluster_matching.bijective$F1.UnweightedMean),
      'Mean Precision Across Matches (Bijective)'                        = as.vector(cluster_matching.bijective$Precision.UnweightedMean),
      'Mean Recall Across Matches (Bijective)'                           = as.vector(cluster_matching.bijective$Recall.UnweightedMean),
      'F1 Per Match (Bijective)'                                         = cluster_matching.bijective$F1.PerMatch,
      'Precision Per Match (Bijective)'                                  = cluster_matching.bijective$Precision.PerMatch,
      'Recall Per Match (Bijective)'                                     = cluster_matching.bijective$Recall.PerMatch,
        
      'Label-Cluster Matching (Relaxed, Fixed Label)'                     = matches_argmax_f1_cluster,
      'Mean F1 Across Matches (Relaxed, Fixed Label)'                     = as.vector(cluster_matching.argmax_f1_cluster$F1.UnweightedMean),
      'Mean Precision Across Matches (Relaxed, Fixed Label)'              = as.vector(cluster_matching.argmax_f1_cluster$Precision.UnweightedMean),
      'Mean Recall Across Matches (Relaxed, Fixed Label)'                 = as.vector(cluster_matching.argmax_f1_cluster$Recall.UnweightedMean),
      'F1 Per Match (Relaxed, Fixed Label)'                               = cluster_matching.argmax_f1_cluster$F1.PerMatch,
      'Precision Per Match (Relaxed, Fixed Label)'                        = cluster_matching.argmax_f1_cluster$Precision.PerMatch,
      'Recall Per Match (Relaxed, Fixed Label)'                           = cluster_matching.argmax_f1_cluster$Recall.PerMatch,
                         
      'Label-Cluster Matching (Relaxed, Fixed Cluster)'                   = matches_argmax_f1_gate,
      'Mean F1 Across Matches (Relaxed, Fixed Cluster)'                   = as.vector(cluster_matching.argmax_f1_gate$F1.UnweightedMean),
      'Mean Precision Across Matches (Relaxed, Fixed Cluster)'            = as.vector(cluster_matching.argmax_f1_gate$Precision.UnweightedMean),
      'Mean Recall Across Matches (Relaxed, Fixed Cluster)'               = as.vector(cluster_matching.argmax_f1_gate$Recall.UnweightedMean),
      'F1 Per Match (Relaxed, Fixed Cluster)'                             = cluster_matching.argmax_f1_gate$F1.PerMatch,
      'Precision Per Match (Relaxed, Fixed Cluster)'                      = cluster_matching.argmax_f1_gate$Precision.PerMatch,
      'Recall Per Match (Relaxed, Fixed Cluster)'                         = cluster_matching.argmax_f1_gate$Recall.PerMatch,
                         
      'Davies-Bouldin Index'                                              = davies_bouldin,
      
      'Adjusted Rand Index'                                               = adjusted_rand_index#,
      #'Normalised Mutual Information'                                     = nmi
    )
  } else if (stability == 'bootstrap') {
    idx_full <- length(result$ClusteringVector) # index of run on original data
    
    ## Unsupervised evaluation metrics
    suppressWarnings(davies_bouldin <- sapply(
      1:(idx_full - 1),
      function(idx_iter) {
        cv <- result$ClusteringVector[[idx_iter]]
        ann <- annotation
        tmp_exprs <- exprs
        if (!is.list(tmp_exprs)) {
          tmp_exprs <- tmp_exprs[bootstrap_idcs[[idx_iter]][[1]], ]
          ann <- ann[bootstrap_idcs[[idx_iter]][[1]]]
        } else {
          for (idx_input in seq_along(tmp_exprs)) {
            tmp_exprs[[idx_input]] <- tmp_exprs[[idx_input]][bootstrap_idcs[[idx_iter]][[idx_input]], ]
            ann[[idx_input]] <- ann[[idx_input]][bootstrap_idcs[[idx_iter]][[idx_input]], ]
          }
          cv <- do.call(rbind, cv)
          tmp_exprs <- do.call(rbind, tmp_exprs)
          ann <- do.call(rbind, ann)
        }
        idcs_assigned <- which(!ann %in% unassigned_labels)
        clusterSim::index.DB(x = exprs[idcs_assigned, ], cl = cv[idcs_assigned])$DB
      }
    ))
    
    ## Supervised evaluation metrics
    adjusted_rand_index <- sapply(
      1:(idx_full - 1),
      function(idx_iter) {
        cv <- result$ClusteringVector[[idx_iter]]
        ann <- annotation
        tmp_exprs <- exprs
        if (!is.list(tmp_exprs)) {
          tmp_exprs <- tmp_exprs[bootstrap_idcs[[idx_iter]][[1]], ]
          ann <- ann[bootstrap_idcs[[idx_iter]][[1]]]
        } else {
          for (idx_input in seq_along(tmp_exprs)) {
            tmp_exprs[[idx_input]] <- tmp_exprs[[idx_input]][bootstrap_idcs[[idx_iter]][[idx_input]], ]
            ann[[idx_input]] <- ann[[idx_input]][bootstrap_idcs[[idx_iter]][[idx_input]], ]
          }
          cv <- do.call(rbind, cv)
          tmp_exprs <- do.call(rbind, tmp_exprs)
          ann <- do.call(rbind, ann)
        }
        idcs_assigned <- which(!ann %in% unassigned_labels)
        ari(annotation[idcs_assigned], cv[idcs_assigned])
      }
    )
    nmi <- sapply(
      1:(idx_full - 1),
      function(idx_iter) {
        cv <- result$ClusteringVector[[idx_iter]]
        ann <- annotation
        tmp_exprs <- exprs
        if (!is.list(tmp_exprs)) {
          tmp_exprs <- tmp_exprs[bootstrap_idcs[[idx_iter]][[1]], ]
          ann <- ann[bootstrap_idcs[[idx_iter]][[1]]]
        } else {
          for (idx_input in seq_along(tmp_exprs)) {
            tmp_exprs[[idx_input]] <- tmp_exprs[[idx_input]][bootstrap_idcs[[idx_iter]][[idx_input]], ]
            ann[[idx_input]] <- ann[[idx_input]][bootstrap_idcs[[idx_iter]][[idx_input]], ]
          }
          cv <- do.call(rbind, cv)
          tmp_exprs <- do.call(rbind, tmp_exprs)
          ann <- do.call(rbind, ann)
        }
        idcs_assigned <- which(!ann %in% unassigned_labels)
        aricode::nmi(annotation[idcs_assigned], cv[idcs_assigned])
      }
    )
    
    
    cluster_matching <-
      purrr::map(
        1:(idx_full - 1),
        function(idx_iter) {
          cv <- result$ClusteringVector[[idx_iter]]
          tmp_annotation <- annotation
          if (!is.list(tmp_annotation))
            tmp_annotation <- tmp_annotation[bootstrap_idcs[[idx_iter]][[1]]]
          else
            for (idx_input in seq_along(tmp_annotation))
              tmp_annotation[[idx_input]] <- tmp_annotation[[idx_input]][bootstrap_idcs[[idx_iter]][[idx_input]]]
          Jacc(c1 = tmp_annotation, c2 = cv, obj = 'f1', unassigned = unassigned_labels, generate_plot = FALSE, verbose = FALSE)
        }
      )
    
    cluster_matching.full <- Jacc(
      c1 = annotation, c2 = result$ClusteringVector[[idx_full]], obj = 'f1', unassigned = unassigned_labels, generate_plot = FALSE, verbose = FALSE
    )
    
    sample.cluster_matching.bijective <- cluster_matching.full$Results.Bijective
    sample.cluster_matching.argmax_f1_cluster <- cluster_matching.full$Results.FixedC1
    sample.cluster_matching.argmax_f1_gate <- cluster_matching.full$Results.FixedC2
    
    scores <- list(
      'Stability Analysis'                                                  = 'bootstrap',

      'Label-Cluster Matching (Bijective)'                                  = matrix(c(names(sample.cluster_matching.bijective$Matches), sample.cluster_matching.bijective$Matches), ncol = 2),
      'Mean F1 Across Matches Per Bootstrap (Bijective)'                    = sapply(cluster_matching, function(x) x$Results.Bijective$F1.UnweightedMean),
      'Mean Precision Across Matches Per Bootstrap (Bijective)'             = sapply(cluster_matching, function(x) x$Results.Bijective$Precision.UnweightedMean),
      'Mean Recall Across Matches Per Bootstrap (Bijective)'                = sapply(cluster_matching, function(x) x$Results.Bijective$Recall.UnweightedMean),
      'F1 Per Match (Bijective)'                                            = sample.cluster_matching.bijective$F1.PerMatch,
      'Precision Per Match (Bijective)'                                     = sample.cluster_matching.bijective$Precision.PerMatch,
      'Recall Per Match (Bijective)'                                        = sample.cluster_matching.bijective$Recall.PerMatch,

      'Label-Cluster Matching (Relaxed, Fixed Label)'                        = matrix(c(names(sample.cluster_matching.argmax_f1_cluster$Matches), sample.cluster_matching.argmax_f1_cluster$Matches), ncol = 2),
      'Mean F1 Across Matches Per Bootstrap (Relaxed, Fixed Label)'          = sapply(cluster_matching, function(x) x$Results.FixedC1$F1.UnweightedMean),
      'Mean Precision Across Matches Per Bootstrap (Relaxed, Fixed Label)'   = sapply(cluster_matching, function(x) x$Results.FixedC1$Precision.UnweightedMean),
      'Mean Recall Across Matches Per Bootstrap (Relaxed, Fixed Label)'      = sapply(cluster_matching, function(x) x$Results.FixedC1$Recall.UnweightedMean),
      'F1 Per Match (Relaxed, Fixed Label)'                                  = sample.cluster_matching.argmax_f1_cluster$F1.PerMatch,
      'Precision Per Match (Relaxed, Fixed Label)'                           = sample.cluster_matching.argmax_f1_cluster$Precision.PerMatch,
      'Recall Per Match (Relaxed, Fixed Label)'                              = sample.cluster_matching.argmax_f1_cluster$Recall.PerMatch,

      'Label-Cluster Matching (Relaxed, Fixed Cluster)'                      = matrix(c(names(sample.cluster_matching.argmax_f1_gate$Matches), sample.cluster_matching.argmax_f1_gate$Matches), ncol = 2),
      'Mean F1 Across Matches Per Bootstrap (Relaxed, Fixed Cluster)'        = sapply(cluster_matching, function(x) x$Results.FixedC2$F1.UnweightedMean),
      'Mean Precision Across Matches Per Bootstrap (Relaxed, Fixed Cluster)' = sapply(cluster_matching, function(x) x$Results.FixedC2$Precision.UnweightedMean),
      'Mean Recall Across Matches Per Bootstrap (Relaxed, Fixed Cluster)'    = sapply(cluster_matching, function(x) x$Results.FixedC2$Recall.UnweightedMean),
      'F1 Per Match (Relaxed, Fixed Cluster)'                                = sample.cluster_matching.argmax_f1_gate$F1.PerMatch,
      'Precision Per Match (Relaxed, Fixed Cluster)'                         = sample.cluster_matching.argmax_f1_gate$Precision.PerMatch,
      'Recall Per Match (Relaxed, Fixed Cluster)'                            = sample.cluster_matching.argmax_f1_gate$Recall.PerMatch,

      'Davies-Bouldin Index'                                                 = davies_bouldin,
      
      'Adjusted Rand Index'                                                  = adjusted_rand_index#,
      #'Normalised Mutual Information'                                        = nmi
    )
  }
  
  scores
}