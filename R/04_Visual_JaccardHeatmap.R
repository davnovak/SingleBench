
.JaccardSimilarityHeatmap <- function(
  jt, counts, eval_bijective, eval_fixedc1, eval_fixedc2, unassigned, title, c1_name = 'c1', c2_name = 'c2'
) {
  
  ## Change Jaccard-table column order to get a nice diagonal from the 1:1 matching
  ## Move rows corresponding to unassigned's to the bottom
  
  matches_bijective   <- as.vector(stats::na.exclude(eval_bijective$Matches))
  unassigneds_present <- if (is.null(unassigned)) c() else intersect(unassigned, rownames(jt))
  
  order_c1 <- c(names(eval_bijective$Matches), unassigneds_present)
  order_c2 <- c(matches_bijective, names(eval_fixedc2$Matches)[!names(eval_fixedc2$Matches) %in% matches_bijective])
  
  row_ordering <- match(order_c1, rownames(jt))
  col_ordering <- match(order_c2, colnames(jt))
  
  if (length(col_ordering) > length(matches_bijective)) {
    ## Change Jaccard-table column order of right submatrix (for clusters that were not matched bijectively)
    
    matches_fixedc2       <- eval_fixedc2$Matches
    
    ## Matches that were not made bijectively:
    additional_c2_matches <- names(matches_fixedc2)[!names(matches_fixedc2) %in% matches_bijective]
    additional_c2_matches <- additional_c2_matches[!additional_c2_matches %in% names(matches_fixedc2)[is.na(matches_fixedc2)]]
    
    ## Matches that were not made either way:
    unmatched_c2          <- colnames(jt)[!colnames(jt) %in% matches_bijective & !colnames(jt) %in% additional_c2_matches]
    
    matched_c1 <- matches_fixedc2[match(additional_c2_matches, names(matches_fixedc2))]
    
    ranking <- match(matched_c1, rownames(jt))
    
    r <- names(matched_c1)[order(ranking)]
    o <- match(r, colnames(jt))
    
    others <- which(!colnames(jt) %in% r & !colnames(jt) %in% matches_bijective)
    
    col_ordering[(length(matches_bijective) + 1):length(col_ordering)] <- c(o, others)
  }
  
  idcs_bijectively_matched <- which(!is.na(eval_bijective$Matches))
  others <- which(is.na(eval_bijective$Matches))
  if (length(others) > 0) {
    row_ordering <- c(idcs_bijectively_matched, others)
    if (length(unassigneds_present) > 0)
      row_ordering <- c(row_ordering, which(rownames(jt) %in% unassigneds_present))
  }
    
  row_ordering <- row_ordering[!is.na(row_ordering)]
  
  jt     <- jt[row_ordering, col_ordering]
  counts <- counts[row_ordering, col_ordering]
  
  ## Retrieve Precision, Recall and F1 scores for each match for each matching approach
  ## Format them as sidebar annotations (data frames accepted by pheatmap)
  
  annot_row <- data.frame(
    Precision = eval_fixedc1$Precision.PerMatch,
    Recall = eval_fixedc1$Recall.PerMatch,
    F1 = eval_fixedc1$F1.PerMatch,
    'Precision_Bijective' = eval_bijective$Precision.PerMatch,
    'Recall_Bijective' = eval_bijective$Recall.PerMatch,
    'F1_Bijective' = eval_bijective$F1.PerMatch
  )
  annot_row <- annot_row[rowSums(annot_row) > 0, ]
  if (length(unassigneds_present) > 0)
    for (count in seq_along(unassigneds_present)) annot_row <- rbind(annot_row, rep(0, times = 6))
  rownames(annot_row) <- rownames(jt)
  
  annot_col <- data.frame(
    Precision = eval_fixedc2$Precision.PerMatch,
    Recall = eval_fixedc2$Recall.PerMatch,
    F1 = eval_fixedc2$F1.PerMatch
  )
  rownames(annot_col) <- names(eval_fixedc2$Matches)
  
  ## Produce the basic heatmap
  
  ph <- pheatmap::pheatmap(
    jt, main = title, legend = FALSE, display_numbers = counts, fontsize_number = 6, cluster_cols = FALSE, cluster_rows = FALSE,
    annotation_row = annot_row, annotation_col = annot_col, annotation_legend = FALSE, border_color = 'grey',
    gaps_row = min(nrow(jt) - length(unassigneds_present), sum(!is.na(matches_bijective))),
    gaps_col = min(length(matches_bijective), nrow(jt) - length(unassigneds_present)),
    number_color = 'black', silent = TRUE,
    annotation_colors =
      list(
        Precision           = RColorBrewer::brewer.pal(6, 'Greens'),
        Recall              = RColorBrewer::brewer.pal(6, 'Blues'),
        F1                  = RColorBrewer::brewer.pal(6, 'Purples'),
        Precision_Bijective = RColorBrewer::brewer.pal(6, 'Greens'),
        Recall_Bijective    = RColorBrewer::brewer.pal(6, 'Blues'),
        F1_Bijective        = RColorBrewer::brewer.pal(6, 'Purples')
      )
  )
  
  ## Access and tweak graphical parameters of the grob corresponding to cells of the heatmap
  
  grob_classes <- lapply(ph$gtable$grobs, class)
  idx_maingrob <- which(sapply(grob_classes, function(cl) 'gTree' %in% cl))[1]
  grob_names   <- names(ph$gtable$grobs[[idx_maingrob]]$children)
  rects_name   <- grob_names[grep('rect', grob_names)]
  text_name   <- grob_names[grep('text', grob_names)]
  gp           <- ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp # graphical parameters of the rectangles
  
  ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp$col <- gp$fill
  ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp$lwd <-
    matrix(0, nrow = nrow(gp$fill), ncol = ncol(gp$fill), dimnames = list(rownames(gp$fill), colnames(gp$fill)))
  
  ## Highlight cells corresponding to columns chosen by fixed-c1 matching by drawing pink frames around them
  ## To that end, create a new 'overlay' grob with the pink frames to the grob tree
  
  overlay_grob         <- ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]
  overlay_grob$gp$fill <- NA

  m                    <- data.frame(Group.c1 = names(eval_fixedc1$Matches), Group.c2 = eval_fixedc1$Matches)
  m                    <- m[!is.na(m$Group.c2),]
  for (idx_match in seq_len(nrow(m))) {
    if (m$Group.c1[idx_match] %in% rownames(overlay_grob$gp$col) && m$Group.c2[idx_match] %in% colnames(overlay_grob$gp$col)) {
      overlay_grob$gp$col[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- '#ff96e1'
      overlay_grob$gp$lwd[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- 2
    }
  }
  ph$gtable$grobs[[idx_maingrob]]$children$overlay_fixedc1 <- overlay_grob
  
  ## Highlight cells corresponding to rows chosen by fixed-c2 matching by drawing light-green frames around them
  ## To that end, create a new 'overlay' grob with the light-green frames to the grob tree (smaller than the pink ones, to avoid overlaps)
  
  overlay_grob         <- ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]
  overlay_grob$gp$fill <- NA
  overlay_grob$height  <- overlay_grob$height - grid::unit(0.6, 'strheight', as.character(overlay_grob$height))
  overlay_grob$width   <- overlay_grob$width - grid::unit(0.02, 'strwidth', as.character(overlay_grob$width))
  m                    <- data.frame(Group.c2 = names(eval_fixedc2$Matches), Group.c1 = eval_fixedc2$Matches)
  m                    <- m[!is.na(m$Group.c1),]
  for (idx_match in seq_len(nrow(m))) {
    if (m$Group.c1[idx_match] %in% rownames(overlay_grob$gp$col) && m$Group.c2[idx_match] %in% colnames(overlay_grob$gp$col)) {
      overlay_grob$gp$col[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- '#85ffa5'
      overlay_grob$gp$lwd[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- 2
    }
  }
  ph$gtable$grobs[[idx_maingrob]]$children$overlay_fixedc2 <- overlay_grob
  
  text_grob_name    <- ph$gtable$grobs[[idx_maingrob]]$childrenOrder['text']
  other_grobs_names <- ph$gtable$grobs[[idx_maingrob]]$childrenOrder[ph$gtable$grobs[[idx_maingrob]]$childrenOrder != text_grob_name]
  ph$gtable$grobs[[idx_maingrob]]$childrenOrder <- c(other_grobs_names, 'overlay_fixedc1', 'overlay_fixedc2', text_grob_name)
  ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp$lwd <- 1
  
  ## Make the labels for row and column annotation (with stats for each match) prettier
  
  idx_grob <- which(sapply(ph$gtable$grobs, function(x) !is.null(x$label) && 'Precision_Bijective' %in% x$label))
  ph$gtable$grobs[[idx_grob]]$label   <- c(paste0('Precision (Fixed ', c1_name, ')'), paste0('Recall (Fixed ', c1_name, ')'), paste0('F1 (Fixed ', c1_name, ')'), 'Precision (bijective)', 'Recall (bijective)', 'F1 (bijective)')
  ph$gtable$grobs[[idx_grob]]$gp$font <- c(2L, 2L, 2L, 3L, 3L, 3L)
  ph$gtable$grobs[[idx_grob]]$rot <- 290
  
  idx_grob <- which(sapply(ph$gtable$grobs, function(x) !is.null(x$label) && 'Precision' %in% x$label))
  ph$gtable$grobs[[idx_grob]]$label   <- c(paste0('Precision (Fixed ', c2_name, ')'), paste0('Recall (Fixed ', c2_name, ')'), paste0('F1 (Fixed ', c2_name, ')'))
  
  ## The labels of those groups in c1 which belong to unassigned will be in italics
  
  if (length(unassigneds_present) > 0) {
    idx_grob <- which(sapply(ph$gtable$grobs, function(x) all(rownames(jt) %in% x$label)))
    font_vec <- rep(1L, length(ph$gtable$grobs[[idx_grob]]$label))
    font_vec[length(font_vec):(length(font_vec - length(unassigneds_present)) + 1)] <- 3L
    ph$gtable$grobs[[idx_grob]]$gp$font <- font_vec
  }
  
  ## Put c2 groups used in bijective matching in bold and rotate the labels of c2 groups
  
  idx_grob <- which(sapply(ph$gtable$grobs, function(x) all(colnames(jt) %in% x$label)))
  if (length(idx_grob) > 1)
    idx_grob <- which(sapply(ph$gtable$grobs, function(x) all(rownames(jt) %in% x$label)))
  font_vec <- rep(1L, length(ph$gtable$grobs[[idx_grob]]$label))
  font_vec[1:length(eval_fixedc1$Matches)] <- 2L
  
  ph$gtable$grobs[[idx_grob]]$gp$font <- font_vec
  ph$gtable$grobs[[idx_grob]]$rot     <- 0
  ph$gtable$grobs[[idx_grob]]$hjust   <- 0.5
  ph$gtable$grobs[[idx_grob]]$vjust   <- 1.5
  
  hmap <- cowplot::plot_grid(NULL, ph[[4]], NULL, nrow = 3, rel_heights = c(0.2, 15, 1))
  explanation <-
    cowplot::plot_grid(
      NULL,
      cowplot::plot_grid(
        NULL, grid::textGrob(paste0('One-to-one (bijective) mapping corresponds to matrix diagonal.\nMapping of groups from ', c2_name, ' to groups from ', c1_name, ' corresponds to green frames.\nMapping of groups from ', c1_name, ' to groups from ', c2_name, ' corresponds to pink frames.'), gp = grid::gpar(fontsize = 10)), ncol = 2, rel_widths = c(0.3, 0.7)
      ), nrow = 2, rel_heights = c(5, 1)
    )
  cowplot::ggdraw(hmap) + cowplot::draw_plot(explanation)
  
}

#' Create a comprehensive similarity heatmap plot
#'
#' Creates a comprehensive plot showing how the results of clustering map onto manually labelled populations using different matching strategies.
#'
#' Using this, you can look at which cell populations were identified correctly or incorrectly by automated clustering and what kinds of mistakes the clustering set-up made.
#' Moreover, the heatmap shows eventual differences between matching clusters to populations bijectively (one-to-one) and taking the best cluster for each population (and vice versa).
#' A comparison of precision, recall and F1 scores for each match is also provided.
#'
#' To create the plot, you need to specify which sub-pipeline and *n*-parameter iteration you want to look at.
#' If you choose multiple *n*-parameter iterations, a list of plots is returned.
#'
#' @param benchmark object of type \code{Benchmark}, as generated by the constructor \code{Benchmark} and evaluated using \code{Evaluate}
#' @param idx.subpipeline integer value: index of sub-pipeline that includes a clustering step
#' @param idx.n_param integer: index of *n*-parameter iteration of interest. Default value is \code{NULL}
#' @param idx.run integer: if clustering was run repeatedly for stability analysis, which run should be used to plot the heatmap. Default value is \code{1}
#'
#' @export
PlotSimilarityHeatmap <- PlotJaccardHeatmap <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = 1
) {
  
  .PlotJaccardHeatmap.ValidityChecks(environment())
  .PlotClustering.ValidityChecks(environment())
  
  annotation <- GetAnnotation(benchmark, concatenate = TRUE)
  
  no_npar <- FALSE
  if (length(idx.n_param) == 0) {
    no_npar <- TRUE
    idx.n_param <- 1
  }
  
  if (no_npar)
    idx.n_param <- NULL
  
  clustering <- GetClustering(benchmark, idx.subpipeline, idx.n_param, all_runs = TRUE)
  if (benchmark$stability %in% c('repeat', 'bootstrap')) {
    clustering <- clustering[[idx.run]]
  }
  name       <- GetNParameterIterationName(benchmark, idx.subpipeline, idx.n_param)
  scores     <- GetClusteringScores(benchmark, idx.subpipeline, idx.n_param)
  if (length(scores) == 1) scores <- scores[[1]]
  
  matches <- scores$`Label-Cluster Matching (Bijective)`
  f1 <- scores$`F1 Per Match (Bijective)`
  pr <- scores$`Precision Per Match (Bijective)`
  re <- scores$`Recall Per Match (Bijective)`
  if (benchmark$stability %in% c('repeat', 'bootstrap')) {
    matches <- matches[matches$Run == idx.run, -1]
    f1 <- f1[[idx.run]]
    pr <- pr[[idx.run]]
    re <- re[[idx.run]]
  }
  eval_bijective <- list(
    Matches            = matches,
    F1.PerMatch        = f1,
    Precision.PerMatch = pr,
    Recall.PerMatch    = re
  )
  m <- eval_bijective$Matches$Cluster
  names(m) <- eval_bijective$Matches$Population
  eval_bijective$Matches <- m
  
  matches <- scores$`Label-Cluster Matching (Relaxed, Fixed Label)`
  f1 <- scores$`F1 Per Match (Relaxed, Fixed Label)`
  pr <- scores$`Precision Per Match (Relaxed, Fixed Label)`
  re <- scores$`Recall Per Match (Relaxed, Fixed Label)`
  if (benchmark$stability %in% c('repeat', 'bootstrap')) {
    matches <- matches[matches$Run == idx.run, -1]
    f1 <- f1[[idx.run]]
    pr <- pr[[idx.run]]
    re <- re[[idx.run]]
  }
  eval_fixedc1 <- list(
    Matches            = matches,
    F1.PerMatch        = f1,
    Precision.PerMatch = pr,
    Recall.PerMatch    = re
  )
  m <- eval_fixedc1$Matches$Cluster
  names(m) <- eval_fixedc1$Matches$Population
  eval_fixedc1$Matches <- m
  
  matches <- scores$`Label-Cluster Matching (Relaxed, Fixed Cluster)`
  f1 <- scores$`F1 Per Match (Relaxed, Fixed Cluster)`
  pr <- scores$`Precision Per Match (Relaxed, Fixed Cluster)`
  re <- scores$`Recall Per Match (Relaxed, Fixed Cluster)`
  if (benchmark$stability %in% c('repeat', 'bootstrap')) {
    matches <- matches[matches$Run == idx.run, -1]
    f1 <- f1[[idx.run]]
    pr <- pr[[idx.run]]
    re <- re[[idx.run]]
  }
  eval_fixedc2 <- list(
    Matches            = matches,
    F1.PerMatch        = f1,
    Precision.PerMatch = pr,
    Recall.PerMatch    = re
  )
  m <- eval_fixedc2$Matches$Population
  names(m) <- eval_fixedc2$Matches$Cluster
  eval_fixedc2$Matches <- m
  
  jt <- JaccardTable(c1 = annotation, c2 = as.factor(clustering))
  counts <- table(annotation, clustering)
  
  idcs_keep <- rowSums(counts) > 0
  jt <- jt[idcs_keep, ]
  counts <- counts[idcs_keep, ]
  
  if (nrow(jt) == 1)
    stop('Only one population present')
  
  res <- .JaccardSimilarityHeatmap(
    jt, counts, eval_bijective, eval_fixedc1, eval_fixedc2, unassigned = benchmark$unassigned_labels,
    title = paste0('Jaccard similarity heatmap: ', name),
    c1_name = 'Labels', c2_name = 'Clusters'
  )
  
  res
}
