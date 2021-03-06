% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_Visual_Overlays.R
\name{PlotLabelsOverlay}
\alias{PlotLabelsOverlay}
\title{Create a labels-overlay plot}
\usage{
PlotLabelsOverlay(
  benchmark,
  exclude_unassigned = FALSE,
  raster_threshold = 5000,
  plot_title = "Manual population labels"
)
}
\arguments{
\item{benchmark}{object of type \code{Benchmark}}

\item{exclude_unassigned}{logical: if \code{TRUE}, data points that are considered unassigned per manual annotation are omitted. Default value is \code{FALSE}}

\item{raster_threshold}{integer: maximum number of data points for which vector graphics should be used. Default value is \code{5000}}

\item{plot_title}{string: title of the plot. Default value is '\emph{Manual population labels}'}
}
\description{
Creates a plot showing separation of manually annotated populations using previously generated 2-dimensional layout of input expression data.
For large numbers of data points, this function uses raster graphics to visualise the layout faster.
}
\seealso{
\itemize{
\item \strong{\code{Benchmark}}: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
\item \strong{\code{Evaluate}}: runs all benchmark sub-pipelines and scores the performance of each tool
\item \strong{\code{AddLayout}}: allows you to add a separate 2-dimensional layout of the input dataset or to use an existing projection (produced in the evaluation) as a visualisation layout.
}
}
