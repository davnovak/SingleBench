<img src="logo.png" alt="SingleBench" width="350"/>

This is a pre-release of `SingleBench`.
It is an `R` framework and API for

* applying projection (de-noising/dimension-reduction) and clustering to cytometry (and other single-cell) datasets,

* benchmarking performance,

* optimising parameters.

Since there is little consensus regarding the best quality measures for clustering, `SingleBench` gives runs a battery of tests and shows which cell populations are better or worse detected.

`SingleBench` is easy to extend with new projection and clustering tools.
It uses `HDF5` files to store intermediate results while saving memory and preventing loss of data due to interruption.
Repeated runs of clustering (for stability analysis) can be parallelised.

## Installation

Install `SingleBench` using the package `devtools` from your `R` console:

```
devtools::install_github('davnovak/SingleBench')
```

This will install required dependencies.
However, to run individual projection/clustering tools with `SingleBench`, you will need to install them also.

## Usage

After installation, the package needs to be loaded.
This exports projection/clustering method wrappers to the global namespace.

```
library(SingleBench)
```

A start-up message lists available wrappers and potential missing packages.

We start by setting up a pipeline: each pipeline is made up of subpipelines, which are combinations of tools and their parameters.
Each subpipeline can have one or both of two modules: *projection* and *clustering*.

A single subpipeline for data smoothing (a simple denoising algorithm), followed by clustering is created using this syntax:

```
subpipelines <- list()
subpipelines[[1]] <- 
    Subpipeline(
        projection = Module(Fix('pSmooth', k = 50, mode = 'Local'), n_param = 'n_iter'),
        clustering = Module(Fix('FlowSOM', grid_width = 10, grid_height = 10), n_param = 'n_clusters')
    )
```
 
`Fix` associates a method with input parameter values (`pSmooth` with parameter values of `k` and `mode`).
We plug this into a `Module`, which *may* specify an *n*-parameter: this is a numeric parameter that is kept variable for executing parameter sweeps over a range of values (`pSmooth` has `n_iter` here, and `FlowSOM` has `n_clusters`).

Since we might want to transform our data in multiple steps to make it more amenable to clustering, the projection can chain multiple modules back-to-back (instead of a single `Module`, use a `list` of them).

This is how we set up the range of *n*-parameter values for both of our modules:
 
```
n_params <- list()
n_params[[1]] <- list(
    projection = rep(c(NA, 1, 2, 3), each = 2),
    clustering = c(30, 40)
)
```

If *n*-parameter values for both projection and for clustering are given, they will be aligned (like any two vectors would be with `cbind`).
 
Setting an *n*-parameter value of projection to `NA` will cause omitting the projection step in that iteration.
If an *n*-parameter value is re-used in the projection step, results get recycled.
 
To create a second subpipeline that uses the same projection method, use `CloneFrom`:

```
subpipelines[[2]] <-
    Subpipeline(
        projection = CloneFrom(1),
        clustering = Module(Fix('Depeche'), n_param = 'fixed_penalty')
    )      
n_params[[2]] <- list(
	 projection = rep(c(NA, 1, 2, 3), times = 3),
    clustering = rep(c(0.1, 0.3, 0.5), each = 4),
    expand = TRUE
)
```

If any projection result is needed multiple times (the same projection module and, if *n*-parameter is specified, the same *n*-parameter value is set), it is only produced once and recycled later.

Finally, we can construct a `Benchmark` object with these settings.
For demonstration, we use the package `HDCytoData` to retrieve a cytometry dataset as a `SummarizedExperiment` object for use in our benchmark.
We will only want to keep *type* markers and use an *asinh* (cofactor = 5) transformation.
We also specify that there exists a label for unannotated cells ('*unassigned*').
To evaluate stability, we will run each clustering step 5 times.

```
library(HDCytoData)
d <- Levine_32dim_SE() # for speed-up, use square-bracket subsetting to downsample the input
b <- Benchmark(
  input              = d,
  transform          = 'asinh',
  transform_cofactor = 5,
  input_marker_types = 'type',
  unassigned_labels  = 'unassigned',
  subpipelines       = subpipelines,
  n_params           = n_params,
  stability_repeat   = 5
)
```

This will set up our pipeline, load input data and create an auxiliary HDF5 to store all intermediate products of evaluation.

The next step is to check if the pipeline set-up is correct and evaluate (run) our benchmark pipeline. 

```
print(b)
Evaluate(b)
```

As the benchmark pipeline is being evaluated, you will get progress messages.
By default, unsupervised and supervised evaluation scores are computed for the clustering step.
By setting `score_projections = TRUE` in `Benchmark`, the projection step is scored as a dimension reduction step (measuring structure preservation in the lower-dimensional embedding).

Optionally, you can generate a 2-dimensional layout of your data for visualisation purposes.

```
AddLayout(b, method = Fix('UMAP', latent_dim = 2))
```

To extract results of your benchmark, you can use the getters `GetAnnotation`, `GetClustering`, `GetProjection` and `GetLayout`.

```
print(b)                                 # get a read-out of setup and evalution metrics
gates <- GetAnnotation(b) # get vector of manual gates
clus <- GetClustering(b, idx.subpipeline = 1, idx.n_param = 2)
                                         # get vector of cluster indices (for n_clusters=40)
layout <- GetLayout(b)    # extract the 2-d layout

cs <- GetClusterSizes(b, 1, 2)
ps <- GetPopulationSizes(b)
cr <- GetRMSDPerCluster(b, 1, 2)
pr <- GetRMSDPerPopulation(b, 1, 2)
```

You can also view the values of evaluation metrics or create informative plots.

```
PlotJaccardHeatmap(b, 1)
PlotLabelsOverlay(b)
PlotClustersOverlay(b, 1, 3)
PlotPopulationHeatmap(b)
PlotClusterHeatmap(b, 1, 2)
PlotCompositionMap(b, 1)
PlotNParameterMap_Clustering(b, 1)
ShowReadout_Clustering(b, 1)
PlotRMSDToSize(b, 1, 1)
```

## Documentation

The easiest way to learn about other features of `SingleBench` (stability analysis, parallelisation, HPC use and others) is to read the documentation for each of the user-level functions shown here.

To see all functions exported from `SingleBench`, navigate to tab *Packages* in RStudio, click on `SingleBench` and look through them in your *Help* tab.

## Method wrappers

Projection and clustering tools are specified in the benchmark set-up by **wrappers**. These are functions with a prescribed signature that allow `SingleBench` to interact with your tools of interest. They are loaded into your global namespace when you load `SingleBench`. The recipes for wrappers themselves are `R` scripts located in `inst/extdata/wrappers_projection` and `inst/extdata/wrappers_clustering`. To generate a wrapper, use the function `WrapTool`. You will find everything you need to know about writing a new wrapper in the documentation for `WrapTool`. You can also look through existing wrappers for inspiration.

## References & acknowledgements

The approximate implementation of *k*-nearest-neighbour graph building using the annoy algorithm is adopted from the *ivis* dimension-reduction tool.

*Szubert, B., Cole, J.E., Monaco, C. et al. Structure-preserving visualisation of high dimensional single-cell datasets. Sci Rep 9, 8914 (2019). https://doi.org/10.1038/s41598-019-45301-0*

The bijective matching of clusters and populations is based on this 2016 clustering benchmark:

*Weber, L.M. and Robinson, M.D. (2016), Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data. Cytometry, 89: 1084-1096. https://doi.org/10.1002/cyto.a.23030*

## Future plans

I am working on extending and describing various de-noising techniques for cytometry data and on using variational autoencoders (and useful priors) for creating more clustering-amenable projections of high-dimensional data and producing semi-synthetic data.

Sampling parameters and sorting results using a multi-perspective assessment via Pareto-efficiency could be implemented into `SingleBench` also; from:

*Putri GH, Koprinska I, Ashhurst TM, King NJC, Read MN. Using single-cell cytometry to illustrate integrated multi-perspective evaluation of clustering algorithms using Pareto fronts. Bioinformatics. 2021 Jan 28:btab038. doi: 10.1093/bioinformatics/btab038.*

