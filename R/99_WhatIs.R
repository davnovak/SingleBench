
#' Inspect a specific cluster or annotated population
#'
#' For a select group of cells (cluster or annotated population), signal intensiry for each parameter (marker) is plotted, as well as the heatmap of inter-parameter Spearman correlations.
#' If there are multiple clusters or populations, these plots are generated for all of them.
#' 
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of subpipeline (if a cluster is to be retrieved). Default value is \code{NULL}
#' @param idx.n_param integer: *n*-parameter iteration index (if a cluster is to be retrieved and *n*-parameter was specified in the benchmark set-up). Default value is \code{NULL}
#' @param idx.run integer: which clustering run to use if clustering was run repeatedly. Default value is \code{1}
#' @param cluster integer or integer vector: index (indices) of cluster of interest (if a cluster is to be retrieved). Default value is \code{NULL}
#' @param population string or string vector: name (names) of population of interest (if a population is to be retrieved). Default value is \code{NULL}
#' @param pointsize numeric: size of points in expression scatter plot. Default value is \code{0.1}
#'
#' @export
WhatIs <- function(
  benchmark, idx.subpipeline = NULL, idx.n_param = NULL, idx.run = 1, cluster = NULL, population = NULL, pointsize = 0.1
) {
  .WhatIs.ValidityChecks(environment())

  # title <- GetNParameterIterationName(benchmark, idx.subpipeline, idx.n_param)
  
  palette <- c(
    RColorBrewer::brewer.pal(8, 'Dark2'),
    RColorBrewer::brewer.pal(12,'Paired'),
    RColorBrewer::brewer.pal(12,'Set3')
  )
  
  exprs <- GetExpressionMatrix(benchmark, concatenate = TRUE)
  
  f_plot <- function(cell_idcs, title, subtitle) {
    d <- tidyr::pivot_longer(data.frame(exprs[cell_idcs, ]), cols = tidyselect::everything())
    popnum <- as.integer(as.factor(d$name))
    popcol <- palette[popnum]
    set.seed(1); popnum <- popnum + rnorm(nrow(d), 0, 0.1)
    
    p_par <- ggplot() +
      scattermore::geom_scattermost(xy = cbind(popnum, d$value), color = popcol, pointsize = pointsize) +
      scale_x_discrete('Parameter', limits = unique(d$name)) +
      theme_grey() +
      theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
    
    cor <- Hmisc::rcorr(exprs[cell_idcs, ], type = 'spearman')
    r <- cor$r
    diag(r) <- NA
    r[which(upper.tri(r))] <- NA
    p_cor <- pheatmap::pheatmap(r, silent = TRUE, cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 9)
    
    grob_classes <- purrr::map(p_cor$gtable$grobs, class)
    idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
    grob_names <- names(p_cor$gtable$grobs[[idx_grob]]$children)
    idx_rect <- grob_names[grep('rect', grob_names)][1]
    p_cor$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- 'black'
    p_cor$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0.1
    
    cowplot::plot_grid(
    #  cowplot::plot_grid(
    #    grid::textGrob(title, gp = grid::gpar(fontface = 'bold', fontsize = 10, fill = 'yellow')),
        grid::textGrob(subtitle, gp = grid::gpar(lty = 2, fontsize = 10)),
    #    nrow = 2
    #  ),
      cowplot::plot_grid(p_par, p_cor[[4]]),
      
      nrow = 2,
      rel_heights = c(1, 9)
    )
  }
  
  res <- vector(mode = 'list', length = length(cluster) + length(population))
  idx_res <- 1
  
  for (cl in cluster) {
    cv <- GetClustering(benchmark, idx.subpipeline, idx.n_param, all_runs = TRUE, concatenate = TRUE)
    if (is.list(cl)) {
      cv <- cv[[if (is.null(idx.run)) 1 else idx.run]]
    }
    idcs <- which(cv == cl)
    subtitle <- paste0('Cluster: ', cl, ' ... ', length(idcs), if (length(idcs) > 1) ' cells' else ' cell')
    if (is.null(idcs))
      stop(paste0('Cluster ', cl, ' not found'))
    
    res[[idx_res]] <- f_plot(idcs, title, subtitle)
    idx_res <- idx_res + 1
  }
  
  if (!is.null(population))
    ann <- GetAnnotation(benchmark, concatenate = TRUE)
  for (pop in population) {
    idcs <- which(ann == pop)
    subtitle <- paste0('Population: ', pop, ' ... ', length(idcs), if (length(idcs) > 1) ' cells' else ' cell')
    if (is.null(idcs))
      stop('Population "', pop, '" not found')
    
    res[[idx_res]] <- f_plot(idcs, title, subtitle)
    idx_res <- idx_res + 1
  }

  cowplot::plot_grid(plotlist = res, ncol = 1)
}
