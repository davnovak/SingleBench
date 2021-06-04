
#' Create a labels-overlay plot
#'
#' Creates a plot showing separation of manually annotated populations using previously generated 2-dimensional layout of input expression data.
#' For large numbers of data points, this function uses raster graphics to visualise the layout faster.
#'
#' @param benchmark object of type \code{Benchmark}
#' @param exclude_unassigned logical: if \code{TRUE}, data points that are considered unassigned per manual annotation are omitted. Default value is \code{FALSE}
#' @param raster_threshold integer: maximum number of data points for which vector graphics should be used. Default value is \code{5000}
#' @param plot_title string: title of the plot. Default value is '*Manual population labels*'
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#' 
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#' 
#' * **\code{AddLayout}**: allows you to add a separate 2-dimensional layout of the input dataset or to use an existing projection (produced in the evaluation) as a visualisation layout.
#'
#' @export
PlotLabelsOverlay <- function(
  benchmark,
  exclude_unassigned = FALSE,
  raster_threshold = 5000,
  plot_title = 'Manual population labels'
) {
  if (class(benchmark) != 'Benchmark') stop('"benchmark" is not of class "Benchmark"')
  if (!benchmark$layout_available) stop('No 2-dimensional layout available')
  
  layout <- as.data.frame(GetLayout(benchmark), concatenate = TRUE)
  annotation <- GetAnnotation(benchmark, concatenate = TRUE)
  
  if (exclude_unassigned) {
    levels(annotation)[levels(annotation) %in% benchmark$unassigned_labels] <- NA
    idcs <- is.na(annotation)
    if (sum(idcs) > 0) {
      annotation <- annotation[!idcs]
      layout <- layout[!idcs, ]
    }
  }
  
  palette <- c(
    RColorBrewer::brewer.pal(8, 'Dark2'),
    RColorBrewer::brewer.pal(12,'Paired'),
    RColorBrewer::brewer.pal(12,'Set3')
  )
  
  Component1 <- Component2 <- NULL # avoid build check warning
  
  ggplot(layout, aes(x = Component1, y = Component2, col = annotation)) +
    labs(col = 'label') + ggtitle(plot_title) + theme_dark() + theme(plot.margin = unit(c(1, 1, 1, 1), 'cm')) +
    if (nrow(layout) > raster_threshold)
      scattermore::geom_scattermore(pointsize = 0.005)
    else
      geom_point(size = 0.5, alpha = 0.85)
}

#' Create a clusters-overlay plot
#'
#' Creates a plot showing separation of generated clusters using previously generated 2-dimensional layout of input expression data.
#' For large numbers of data points, this function uses raster graphics to visualise the layout faster.
#' 
#' You need to specify a single sub-pipeline and one or more *n*-parameter iteration by index.
#' If multiple *n*-parameter iterations are chosen, a list of plots will be generated.
#' 
#' You can choose a custom title for your plot.
#' The plot subtitle is always the *n*-parameter iteration name.
#'
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of sub-pipeline that includes a clustering step
#' @param idcs.n_param integer or vector of integers: indices of *n*-parameter (cluster count) of interest. Default value is \code{NULL}, which translates to all *n*-parameter iterations of the given sub-pipeline
#' @param exclude_unassigned logical: if \code{TRUE}, data points that are considered unassigned per manual annotation are omitted. Default value is \code{FALSE}
#' @param raster_threshold integer: maximum number of data points for which vector graphics should be used. Default value is \code{5000}
#' @param plot_title string: title of the plot. Default value is '*Clustering*'
#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#' 
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#' 
#' * **\code{AddLayout}**: allows you to add a separate 2-dimensional layout of the input dataset or to use an existing projection (produced in the evaluation) as a visualisation layout.
#'
#' @export
PlotClustersOverlay <- function(
  benchmark,
  idx.subpipeline,
  idcs.n_param = NULL,
  exclude_unassigned = FALSE,
  raster_threshold = 5000L,
  plot_title = 'Clustering'
) {
  if (class(benchmark) != 'Benchmark') stop('"benchmark" is not of class "Benchmark"')
  if (!benchmark$layout_available) stop('No 2-dimensional layout available')
  if (!benchmark$evaluated_previously) stop('Benchmark pipeline has not been evaluated')
  
  if (is.null(idcs.n_param))
    idcs.n_param <- seq_len(GetNParameterIterationsCount(benchmark, idx.subpipeline))
  
  layout <- as.data.frame(GetLayout(benchmark))
  
  if (exclude_unassigned) {
    annotation <- GetAnnotation(benchmark)
    levels(annotation)[levels(annotation) %in% benchmark$unassigned_labels] <- NA
    idcs <- is.na(annotation)
    if (sum(idcs) > 0)
      layout <- layout[!idcs, ]
  }
  
  palette <- c(
    RColorBrewer::brewer.pal(8, 'Dark2'),
    RColorBrewer::brewer.pal(12,'Paired'),
    RColorBrewer::brewer.pal(12,'Set3')
  )
  
  Component1 <- Component2 <- NULL # avoid build check warning
  
  res <- purrr::map(
    idcs.n_param, function(idx.n_param) {
      clustering <- GetClustering(benchmark, idx.subpipeline, idx.n_param)
      if (exists('idcs') && sum(idcs) > 0)
        clustering <- clustering[!idcs]
      subpipeline_name <- GetNParameterIterationName(benchmark, idx.subpipeline, idx.n_param)
      ggplot(layout, aes(x = Component1, y = Component2, col = as.factor(clustering))) +
        labs(title = plot_title, subtitle = subpipeline_name, col = 'cluster') +
        theme_dark() + theme(plot.margin = unit(c(1, 1, 1, 1), 'cm')) +
        if (nrow(layout) > raster_threshold)
          scattermore::geom_scattermore(pointsize = 0.005)
        else
          geom_point(size = 0.5, alpha = 0.85)
    }
  )
  if (length(idcs.n_param) == 1) res[[1]] else res
}
