
#' \code{Benchmark} object constructor
#'
#' Creates an object of type \code{Benchmark}, used to configure a benchmarking pipeline and input data.
#' This object allows to run projection (dimension-reduction/denoising) and clustering pipelines on high-dimensional single-cell data (such as **flow cytometry**, **mass cytometry**, **CITE-seq**, **scRNA-seq**).
#' 
#' Once set up, use \code{Evaluate} to run the benchmark.
#'
#' # Input data and tools
#'
#' Input to the pipeline (expression data matrix/matrices) is passed as a path to an FCS file (or multiple FCS files) or a \code{flowSet} or \code{SummarizedExperiment} object, using the \code{input} argument.
#' Projection and clustering tools (that should be applied to the data upon evaluation) are passed to the constructor using modules (see the \code{Module} constructor).
#'
#' # Manual labelling of data points
#' 
#' Each input data point needs to have a label, assigning it to a population.
#' Labels can either be included in the input data object or supplied separately in a vector (\code{input_labels}).
#' 
#' Including labels directly is only possible with a \code{flowSet} or \code{SummarizedExperiment} (see convention in package \code{HDCytoData}).
#' With a \code{flowSet}, include a '\code{population_id}' channel with numeric flags and mapping of sorted numeric flags to population names, respectively, as columns of a data frame accessible for each \code{flowFrame} via \code{@description$POPULATION_INFO}.
#' With a \code{SummarizedExperiment}, include a vector of population labels directly as the first column of \code{rowData}.
#' In addition, \code{colData} should contain a \code{channel_name}, \code{marker_name} and (optionally) a \code{marker_class} slot.
#'
#' # Subpipelines and *n*-parameters
#'
#' A pipeline in made up of subpipelines, which are combinations of tools and their parameters.
#' Each subpipeline can have one or both of two modules: projection and clustering.
#' For instance, a single subpipeline for data smoothing, \code{ivis} latent-space projection and subsequent \code{FlowSOM} clustering is created using:
#' 
#'      subpipelines <- list()
#'      subpipelines[[1]] <- 
#'          Subpipeline(
#'              projection = list(
#'                  Module(Fix('smooth', n_iter = 1, k = 50)),
#'                  Module(Fix('ivis'), n_param = 'latent_dim')
#'              ),
#'              clustering = Module(Fix('FlowSOM', grid_width = 25, grid_height = 25), n_param = 'n_clusters')
#'          )
#' 
#' Up to two *n*-parameters (one per projection step and one per clustering step) can be specified.
#' This means that, within the subpipeline, we can iterate over multiple combinations of values of these parameters, gathering results for each of them.
#' These values need to be specified in the \code{Benchmark} constructor.
#' For instance, to try out multiple latent-space dimensions and target cluster counts, we can specify:
#' 
#'       n_params <- list()
#'       n_params[[1]] <- list(
#'           projection = rep(c(20, 15, 10), each = 3),
#'           clustering = rep(c(35, 40, 45), times = 3)
#'      )
#'
#' If you want to set up a second subpipeline that re-uses the projection step, you can simply clone results from the first subpipeline:
#' 
#'      subpipeline[[2]] <-
#'          Subpipeline(
#'              projection = CloneFrom(1),
#'              clustering = Module(Fix('Depeche'), n_param = 'k')
#'          )
#'      
#'      n_params[[2]] <- list(
#'          clustering = rep(c(20, 30, 40), times = 3)
#'      )
#'
#' # Stability analysis of clustering
#'   
#' If you are interested in analysing stability of a clustering algorithm, run the tool multiple times (with a different random seed) or apply it multiple times to a bootstrap of the input data.
#' That way, for your clustering evaluation metrics (scores) you will get a range of values (instead of a single values) for each subpipeline and *n*-parameter iteration.
#' For repeated runs, use the \code{stability_repeat} parameter.
#' For bootstraps, use the \code{stability_bootstrap} parameter.
#' 
#' # Train-and-map
#' 
#' Some projection and clustering methods allow you to train a model on a subset of your input data and map the results onto the rest.
#' If your input dataset consists of multiple FCS files, \code{flowFrame}s or \code{SummarizedExperiment} batches, you can specify the training set indices in \code{projection.training_set} and \code{clustering.training_set}.
#' 
#' # Hierarchical labelling
#' 
#' If the manual labels of your input data are derived from a hierarchy of populations (ie. **gating hierarchy** for cytometry data), you can make use of the entire hierarchy for evaluation purposes.
#' For instance, instead of using a '**CD4+ T cell**' label, you can use '**Lymphocyte/T cell/CD4+ T cell**' (using a path-like label with '\code{/}' as separator).
#' Then, if you apply a clustering tool and match each cluster to a population present in the data, \code{SingleBench} can evaluate the quality of clustering more carefully.
#' Specifically, instead of distinguishing between *match* versus *mismatch*, a scoring matrix is produced which penalises mismatches with different severity.
#' For instance, to misclassify '**Lymphocyte/T cell/CD4+ T cell**' as '**Lymphocyte/T cell**' can be better than misclassifying it as '**Lymphocyte/T cell/CD8+ T cell**', which is still better than misclassifying it as '**Lymphocyte/B cell/Alpha-Beta Mature B Cell**'.
#' 
#' The scoring of each potential mismatch is based on the route from the true label to the predicted label through the label hierarchy tree.
#' To parametrise the hierarchical penalty model, you can set 4 custom values.
#' Firstly, the '*constant generalisation penalty*' \code{g_c} penalises the first step taken in the direction of the tree root and the '*additive generalisation penalty*' \code{g_a} penalises each step in that direction.
#' Secondly, the '*constant misidentification penalty*' \code{m_a} panelises the first step taken in the direction of the tree leaves and the '*additive misidentification penalty*' \code{m_a} penalises each step in that direction.
#' The values of these penalties are positive values, and the sum of penalties for a misclassification get subtracted from \code{1}, which is the score for correct match.
#' By default, \code{g_c = 0}, \code{g_a = 0.2}, \code{m_c = 0.4}, \code{m_a = 0}.
#'
#' If you want to use hierarchical penalties, make sure you use the full labels (paths using '\code{/}' as separator), set the argument \code{hierarchy} to \code{TRUE}.
#' Optionally, set \code{hierarchy.params} to a named list with slots \code{g_c}, \code{g_a}, \code{m_c} and \code{m_a}.
#' Alternatively, instead of computing a scoring matrix for the penalties automatically, you can include a slot \code{scoring_matrix}, which contains the names of all labels as column and row names and a match/mismatch score for each combination of true and predicted label.
#' This scoring matrix will then be used directly and will override any other paramters of the model.
#' 
#' # \code{k}-NN matrix construction
#' 
#' Various dimension-reduction or denoising algorithms may require the construction of a *k*-nearest-neighbour graph.
#' If any of them are included in the pipeline, a single \code{k}-NN matrix will be constructed at the start of evaluation.
#' You can specify the \code{k}-NN building algorithm, the value of \code{k} and the distance metric to use by specifying \code{knn.algorithm}, \code{knn.k}, \code{knn.distance}
#' 
#' # Distance matrix construction
#' 
#' You can decide to score the results of projection steps using full distance matrices (distances of each point to each other point) of the original input matrix and each produced projection.
#' This can be set up in the \code{Evaluate} function.
#' (Beware that construction of a full distance matrix has quadratic complexity with regard to number of points.)
#' To speed things up, you can load a pre-computed distance matrix of your original input data by specifying the value of parameter \code{precomputed_dist}.
#' To compute such a distance matrix you can use \code{as.matrix(stats::dist(.))} or \code{coRanking:::euclidean_C(.)}.
#' (However, know that even if you do provide a pre-computed distance matrix of the original data, distance matrices for each projection produced in the benchmark pipeline will still need to be computed, and this can take a long time for larger datasets.)
#' 
#' # HDF5 integration
#' 
#' All large vectors and matrices generated during benchmark set-up and evaluation are stored in an auxiliary HDF5 file.
#' Erasing it amounts to discarding results of the benchmark.
#' 
#' # Creating new wrappers
#' 
#' To create new wrappers for projection and clustering tools, use \code{WrapTool}.
#' 
#' @param input string vector, \code{flowSet} or \code{SummarizedExperiments}: input dataset given as vector of FCS file paths, \code{flowSet} object or \code{SummarizedExperiment} object
#' @param input_labels optional string or factor vector (or a list of these vectors): if manual labels of data points are not included in the input FCS file(s), \code{flowSet} or \code{SummarizedExperiment}, specify them here. For a single FCS file or a \code{flowSet} or \code{SummarizedExperiment}, use a single vector (of length equal to number of expression matrix rows).
#' If multiple FCS files are used as input, use a list of vectors, one for each of the files (in the same order as the FCS file paths).
#' @param hierarchy logical: whether to use hierarchical penalties for mismatches in scoring clustering results. Default value is \code{FALSE}
#' @param hierarchy.params optional list: parameter values for hierarchical penalties model
#' @param input_features optional numeric or string vector: which columns of the input data should be used in analysis. Default value is \code{NULL}, which causes all columns to be used
#' @param remove_labels optional string vector: names of population labels to exclude from analysis. Default value is \code{NULL}
#' @param unassigned_labels optional string or string vector: names of population labels to handle as 'unassigned' events (not belonging to a manually annotated population) for purposes of computing evaluation metrics. Default value is '\code{unassigned}'
#' @param compensation_matrix optional numeric matrix: compensation matrix to use if input is uncompensated flow cytometry data. See \code{flowCore::compensate}. Default value is \code{NULL}, whereby compensation is not applied
#' @param transform optional string: name of transformation function to apply to input expression data. If specified, it must be one of \code{asinh}, \code{estimateLogicle}. Default value is \code{NULL}, whereby transformation is not applied
#' @param transform_cofactor optional numeric value: cofactor for a transformation function specified in \code{transform}. For example, \code{asinh} with cofactor \code{5} is usually appropriate for mass cytometry data. Default value is \code{NULL}
#' @param transform_features optional numeric or string vector: channels or indices of columns to transform (if \code{transform} is specified). Intersection between the indices given here and the indices retrieved from \code{idcs.input_features} or using \code{input_marker_types} are taken (if either of those two parameters is specified). Default value is \code{NULL}, which translates to 'all features that are used in the analysis'
#' @param input_marker_types optional string: if using \code{SummarizedExperiment} as input, which type(s) of markers (specified in vector \code{marker_type} of \code{colData}) should be used in analysis. This overrides \code{idcs.input_features}. Default value is \code{NULL}
#' @param batch_id optional string: if using \code{SummarizedExperiment} as input, which \code{rowData} column should be used to separate input data into batches if available. Default value is \code{sample_id}
#' @param projection.training_set optional vector of integers: if input dataset comprises multiple samples, which ones should be used as the training set for building a dimension-reduction model. Default value is \code{NULL} (DR is applied to the entire dataset)
#' @param clustering.training_set optional vector of integers: if input dataset comprises multiple samples, which ones should be used as the training set for building a clustering model. Default value is \code{NULL} (clustering is applied to the entire dataset)
#'
#' @param precomputed_knn optional list: a previously computed \code{k}-NN matrix for input data as a list containing slots \code{Indices} and \code{Distances}. Overrides \code{knn.algorithm} Default value is \code{NULL}
#' @param knn.k integer: if data denoising or the evaluation of DR or clustering tools requires generation of a \code{k}-NN matrix, what value of \code{k} should be used. Default value is \code{100}
#' @param knn.algorithm string: if \code{k}-NN matrix needs to be constructed, which algorithm should be used. The options are: '\code{annoy}' '\code{cover_tree}', '\code{kd_tree}' and '\code{brute}'. Default value is '\code{kd_tree}', but it is recommended to use '\code{annoy}' if \code{Python} and \code{reticulate} is available
#' @param knn.distance string: if \code{k}-NN matrix needs to be constructed, which distance metric should be used. The only option currently is '\code{euclidean}'. Default value is '\code{euclidean}'
#' @param knn.features optional numeric or string vector: channels or indices of columns to use in construction of \code{k}-NN matrix. Default value is \code{NULL}, which causes the \code{k}-NN algorithm to use \code{idcs.transform_features}
#' 
#' @param precomputed_dist optional numeric matrix: a previously computed distance matrix. Default value is \code{NULL}
#' 
#' @param subpipelines object of class \code{BenchmarkSubpipeline} generated using the \code{Subpipeline} constructor
#' @param n_params list of named *n*-parameter lists for each subpipeline (eg. \code{list(list(projection = c(20, 10), clustering = c(30, 30)))})
#'
#' @param stability_repeat optional integer value between \code{2} and \code{1000}: number of repeated runs of each clustering step on full input. Default value is \code{NULL}
#' @param stability_bootstrap optional integer value between \code{2} and \code{1000}: number of repeated runs of each clustering step on different bootstraps of input, in addition to one final run on full data. Default value is \code{NULL}. Overrides \code{cluster.stability_randomstart}
#'
#' @param benchmark_name string: benchmark name. Default value is '\code{SingleBench_Benchmark}'
#' @param h5_path string: name of auxiliary HDF5 file to generate. Default value is \code{benchmark_name} (with an '\code{.h5}' extension)
#' @param ask_overwrite logical: if \code{h5_path} is specified and an HDF5 file with the specified name exists already, should the user be asked before the file is overwritten? Default value is \code{TRUE}
#' @param verbose logical value: specifies if informative messages should be produced throughout construction of the \code{Benchmark} object. Default value is \code{TRUE}
#'
#' @seealso
#' 
#' * **\code{Evaluate}**: runs all benchmark subpipelines and scores the performance of each tool
#' 
#' * **\code{AddLayout}**: allows you to add a separate 2-dimensional layout of the input dataset or to use an existing projection (produced in the evaluation) as a visualisation layout.
#'
#' @export
Benchmark <- function(
  input,
  input_labels                  = NULL,
  hierarchy                     = FALSE,
  hierarchy.params              = list('g_c' = 0, 'g_a' = 0.2, 'm_c' = 0.4, 'm_a' = 0),
  input_features                = NULL,
  remove_labels                 = NULL,
  unassigned_labels             = 'unassigned',
  compensation_matrix           = NULL,
  transform                     = NULL,
  transform_cofactor            = NULL,
  transform_features            = NULL,
  input_marker_types            = NULL,
  batch_id                      = 'sample_id',
  projection.training_set       = NULL,
  clustering.training_set       = NULL,
  
  precomputed_knn               = NULL,
  knn.k                         = 100,
  knn.algorithm                 = 'annoy',
  knn.distance                  = 'euclidean',
  knn.features                  = NULL,
  
  precomputed_dist              = NULL,
  
  subpipelines                  = NULL,
  n_params                      = NULL,
  
  stability_repeat              = NULL,
  stability_bootstrap           = NULL,
  
  benchmark_name                = 'SingleBench_Benchmark',
  h5_path                       = paste0(benchmark_name, '.h5'),
  ask_overwrite                 = TRUE,
  verbose                       = TRUE
) {
    ## Resolve input paths
    input_class <- class(input)
    if (length(input_class) == 1 && input_class == 'character')
      input_class <- 'fcs_paths'
    
    ## Check validity of arguments
    .Benchmark.ValidityChecks(environment())
    
    if (input_class == 'fcs_path') {
      input <- sapply(input, function(fname) file.path(getwd(), fname))
      input <- sapply(input, SimplifyPath)
    }

    ## Create benchmark-pipeline S3 object
    b            <- new.env(hash = TRUE)
    b$executable <- FALSE
    b$name       <- benchmark_name
    b$h5_path    <- h5_path
    if (is.null(transform_features))
      transform_features <- input_features
    b$rel_idcs.knn_features       <- knn.features
    b$rel_idcs.transform_features <- transform_features
    
    ## Gather and pre-process input data
    GatherInput(b, input, input_class, input_features, compensation_matrix, transform, transform_features, transform_cofactor, input_labels, input_marker_types, batch_id, remove_labels, verbose)
    b$unassigned_labels <- unassigned_labels[unassigned_labels %in% unlist(b$annotation)]
    
    ## Align subpipelines and n-parameter iterations to resolve re-cycling of intermediates
    AlignSubpipelines(b, subpipelines, n_params)
    
    ## Set up train-and-map
    b$projection.training_set  <- projection.training_set
    b$clustering.training_setp <- clustering.training_set
    
    ## Set up clustering stability analysis
    b$stability <- 'single'
    if (!is.null(stability_bootstrap)) {
      stability_repeat <- NULL
      b$stability <- 'bootstrap'; b$stability.n_iter <- stability_bootstrap
    }
    if (!is.null(stability_repeat)) {
      stability_bootstrap <- NULL
      b$stability <- 'repeat'; b$stability.n_iter <- stability_repeat
    }
    
    ## Set up k-NN graph construction
    b$knn.k         <- knn.k
    b$knn.algorithm <- knn.algorithm
    b$knn.distance  <- knn.distance
    
    b$compute_knn   <- FALSE
    for (idx_subpipeline in seq_along(b$subpipelines)) {
      for (idx_module_proj in seq_along(b$subpipelines[[idx_subpipeline]]$projection$modules)) {
        if (b$subpipelines[[idx_subpipeline]]$projection$modules[[idx_module_proj]]$wrapper_with_parameters$wrapper$uses_knn_graph)
          b$compute_knn <- TRUE
      }
      for (idx_module_clus in seq_along(b$subpipelines[[idx_subpipeline]]$clustering$modules)) {
        b$subpipelines[[idx_subpipeline]]$clustering$modules[[idx_module_clus]]$wrapper_with_parameters$wrapper$uses_knn_graph
      }
    }
    
    b$knn_available <- FALSE
    if (!is.null(precomputed_knn)) {
      if (length(precomputed_knn) != 2 || any(!c('Indices', 'Distances') %in% names(precomputed_knn)))
        stop('Invalid "precomputed_knn" format')
      k <- ncol(precomputed_knn$Indices)
      if (verbose) { .msg('Pipeline will use pre-computed k-NN matrix with k='); .msg_alt(k, '\n') }
      b$knn_available <- TRUE
      b$knn.k <- k
    }
    
    ## Full distance matrix construction (eventually triggered later by the evaluator)
    b$compute_dist <- FALSE
    b$dist_available <- FALSE
    if (!is.null(precomputed_dist)) {
      if (nrow(precomputed_dist) != benchmark$row_count || ncol(precomputed_dist) != benchmark$row_count)
        stop('Invalid "precomputed_dist" dimensions')
      if (verbose) { .msg('Pipeline will use pre-computed distance matrix') }
      b$dist_available <- TRUE
    }
    
    ## Generate bootstrap indices for clustering stability analysis
    bs <- NA
    if (!is.null(stability_bootstrap)) {
      set.seed(1)
      bs <- purrr::map(seq_len(stability_bootstrap), function(i) purrr::map(b$row_count, function(row_count) sample(seq_len(row_count), size = row_count, replace = TRUE)))
      names(bs) <- paste0('BootstrapIndices', formatC(seq_along(bs), width = 4, flag = '0'))
      b$bootstrap_indices <- bs
    }
    
    ## Create auxiliary HDF5 file
    if (verbose) {
      .msg('Creating HDF5 file '); .msg_name(b$h5_path, '\n')
    }
    HDF5_CreateHDF5AndWriteInputs(b, verbose)
    
    if (!is.null(precomputed_knn)) {
      SavekNNMatrix(b, precomputed_knn, verbose)
    }
    if (!is.null(precomputed_dist)) {
      SaveDistanceMatrix(b, precomputed_dist, verbose)
    }
    
    b$exprs <- b$annotation <- b$bootstrap_indices <- NULL
    gc(verbose = FALSE)
    
    b$evaluated_previously <- FALSE
    b$layout_available     <- FALSE
    
    ## Create hierarchical penalties scoring matrix
    b$hierarchy <- hierarchy
    if (b$hierarchy) {
      if (!is.null(hierarchy.params[['scoring_matrix']])) {
        if (verbose) .msg('Creating a hierarchical penalty model for clustering evaluation\n')
        scoring_matrix <- hierarchy.params[['scoring_matrix']]
      } else {
        if (verbose) .msg('Adopting a pre-set scoring matrix with hierarchical penalties for clustering evaluation\n')
        labs <- levels(GetAnnotation(b, concatenate = TRUE))
        if (!is.null(b$unassigned_labels)) labs <- labs[!labs %in% b$unassigned_labels]
        scoring_matrix <- CreatePenaltyScoringMatrix(labs, hierarchy.params[['g_c']], hierarchy.params[['g_a']], hierarchy.params[['m_c']], hierarchy.params[['m_a']])
      }
      SavePenaltyScoringMatrix(b, scoring_matrix)
    }
    
    ## Determine if Python is required
    b$uses_python <-
      any(unlist(c(
          purrr::map(b$wrappers.project, function(x) if (IsNonwrapper(x)) FALSE else x$uses_python),
          purrr::map(b$wrappers.cluster, function(x) if (IsNonwrapper(x)) FALSE else x$uses_python)
      )))
    
    ## Compute median expression per population
    if (verbose)
      .msg('Computing median expression per each labelled population\n')
    exprs <- GetExpressionMatrix(b, concatenate = TRUE)
    annotation <- GetAnnotation(b, concatenate = TRUE)
    
    meds <-
      matrix(
        apply(expand.grid(levels(annotation), b$column_names), 1, function(x) median(exprs[annotation == x[1], x[2], drop = FALSE])),
        ncol = b$column_count,
        dimnames = list(levels(annotation), b$column_names)
      )
    .h5writeLabelMedians(meds, b)
    
    if (!b$executable && verbose) {
      .msg_alt_bad('No subpipelines set up'); .msg(' (you can extract processed input from this object, but not evaluate it)\n')
    }
    
    gc(verbose = FALSE)
    structure(b, class = 'Benchmark')
}
