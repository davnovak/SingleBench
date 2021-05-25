<img src="logo.png" alt="SingleBench" width="350"/>

`SingleBench` is an `R` framework and API for

* applying projection (denoising/dimension-reduction) and clustering to single-cell data (flow cytometry, mass cytometry, CITE-seq, single-cell RNA-seq),

* benchmarking performance of these pipelines on various datasets,

* optimising parameters to boost performance of automated detection of cell states and cell types.

`SingleBench` is easy to extend. It includes a toolbox to create wrappers for new projection and clustering algorithms rapidly.
It makes use of auxialiary `HDF5` files to store intermediate results of evaluation while minimising memory requirements.
Furthermore, it allows users to take advantage of multiple CPU cores or submit jobs to a `gridengine` high-performance cluster.

## Motivation

`SingleBench` is being created for a new extensive benchmark study.
It is somewhat inspired by the Weber et al. (2016) benchmarking study (*DOI: 10.1002/cyto.a.23030*) and other papers that have come out since.
It aims to provide a more powerful evaluation framework (similar in some regards to `snakemake` or `pipeComp`) and explore the notion of denoising and feature extraction prior to clustering.

## Tutorial

After installation, the `SingleBench` package needs to be loaded.
(Upon loading, method wrappers for projection and clustering tools are loaded into the global namespace.)

```
library(SingleBench)
```

The initial step is to set up a benchmark pipeline.
A single pipeline is made up of subpipelines, which are combinations of tools and their parameters.
Each subpipeline can have one or both of two modules: *projection* and *clustering*.
For instance, a single subpipeline for data smoothing (simple denoising algorithm), followed by `ivis` dimension reduction and `FlowSOM` clustering is created using this syntax:

```
subpipelines <- list()
subpipelines[[1]] <- 
    Subpipeline(
        projection = list(
            Module(Fix('smooth', n_iter = 1, k = 50)),
            Module(Fix('ivis'), n_param = 'latent_dim')
        ),
        clustering = Module(Fix('FlowSOM', grid_width = 25, grid_height = 25), n_param = 'n_clusters')
    )
 ```
 
Here, `Fix` takes a method (for which a wrapper exists in the global namespace) and fixes its parameter values.
The wrapper-with-parameters is plugged into a `Module`, which *may* specify an *n*-parameter: this is a single parameter of interest, such as latent-space dimensionality or target number of clusters.
A single subpipeline (or parts thereof) can be run multiple times, iterating over different *n*-parameter values (to see how they affect performance).

Both the projection step and the clustering step of a subpipeline are made up of a module or a list of modules (applied sequentially).
The projection step, as well as the clustering step, may only have a single *n*-parameter each.
To iterate over latent-space dimensionality of 20, 15 and 10 and target cluster count of 35, 40 and 45 (for each of the projections), use the following set-up:
 
```
n_params <- list()
n_params[[1]] <- list(
    projection = rep(c(20, 15, 10), each = 3),
    clustering = rep(c(35, 40, 45), times = 3)
)
```
 
(Alternatively, setting an *n*-parameter value of projection to `NA` will result in omitting the projection step in that iteration.)
 
To create a second subpipeline that re-uses the projection step (and its results), simply clone it from the first one:

```
subpipelines[[2]] <-
    Subpipeline(
        projection = CloneFrom(1),
        clustering = Module(Fix('Depeche'), n_param = 'fixed_penalty')
    )      
n_params[[2]] <- list(
    clustering = rep(c(0.1, 0.3, 0.5), times = 3)
)
```

Finally, we can construct a `Benchmark` object with these settings.
For demonstration, we use the package `HDCytoData` to retrieve a cytometry dataset as a `SummarizedExperiment` object for use in our benchmark.
We will only want to keep *type* markers (per annotation of the dataset) and use an *asinh* (cofactor = 5) transformation.
Furthermore, we specify that there exists a label for unannotated cells ('*unassigned*'), which need to be handled differently when calculating supervised evaluation metrics.

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
  n_params           = n_params
)
```

This will set up our pipeline, load input data and create an auxiliary HDF5 file to store any large chunks of data.

The next step is to evaluate (run) our benchmark pipeline. 

```
Evaluate(b)
```

As the benchmark pipeline is being evaluated, you will get progress messages.
By default, unsupervised and supervised evaluation scores for the clustering step 

Optionally, you can generate a 2-dimensional layout of your data for visualisation purposes.

```
AddLayout(b, method = Fix('UMAP', latent_dim = 2))
```

To extract results of your benchmark, you can `print` your `Benchmark` object or use the getters `GetAnnotation`, `GetClustering`, `GetProjection` and `GetLayout`.

```
print(b)                                 # get a read-out of setup and evalution metrics
gates <- GetAnnotation(b) # get vector of manual gates
clus <- GetClustering(b, idx.subpipeline = 1, idx.n_param = 2)
                                         # get vector of cluster indices (for n_clusters=40)
layout <- GetLayout(b)    # extract the 2-d layout
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
```

## Documentation

The easiest way to learn about other features of `SingleBench` (stability analysis, parallelisation, HPC use and others) is to read the documentation for each of the user-level functions shown here.

To see all functions exported from `SingleBench`, navigate to tab *Packages* in RStudio, click on `SingleBench` and look through them in your *Help* tab.

## Method wrappers

Projection and clustering tools are specified in the benchmark set-up by **wrappers**. These are functions with a prescribed signature that allow `SingleBench` to interact with your tools of interest. They are loaded into your global namespace when you load `SingleBench`. The recipes for wrappers themselves are `R` scripts located in `inst/extdata/wrappers_projection` and `inst/extdata/wrappers_clustering`. To generate a wrapper, use the function `WrapTool`. You will find everything you need to know about writing a new wrapper in the documentation for `WrapTool`. You can also look through existing wrappers for inspiration.

