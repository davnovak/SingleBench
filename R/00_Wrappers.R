
IsWrapper <- function(x) 'BenchmarkToolWrapper' %in% class(x)
IsNonwrapper <- function(x) is.null(x) || (is.atomic(x) && is.na(x))
IsProjectionWrapper <- function(x) 'BenchmarkToolWrapper_Reduction' %in% class(x)
IsClusteringWrapper <- function(x) 'BenchmarkToolWrapper_Clustering' %in% class(x)
IsListOfProjectionWrappers <- function(x) is.list(x) && !IsProjectionWrapper(x) && all(sapply(x, function(xx) is.null(xx) || (is.atomic(xx) && is.na(xx)) || IsProjectionWrapper(xx))) && 'BenchmarkToolWrapper_Projection' %in% unique(unlist(sapply(x, class)))
IsListOfClusteringWrappers <- function(x) is.list(x) && !IsClusteringWrapper(x) && all(sapply(x, function(xx) is.null(xx) || (is.atomic(xx) && is.na(xx)) || IsClusteringWrapper(xx))) && 'BenchmarkToolWrapper_Clustering' %in% unique(unlist(sapply(x, class)))
IsClone <- function(x) length(class(x)) == 1 && class(x) == 'Clone' || isTRUE(attr(x, 'IsClone'))

#' Tweak format of a \code{k}-nearest-neighbour graph object
#'
#' Changes format of a \code{k}-nearest-neighbour graph (\code{k}-NNG) object produced by \code{SingleBench}.
#' You will typically want to use this inside \code{WrapTool} to change the format of a pre-computed \code{k}-NNG object for use with a projection or clustering tool.
#' To only keep the \code{Indices} slot of \code{knn} (only keeping a single matrix instead of a list of two matrices), set parameter \code{only_indices} to \code{TRUE}.
#' To modify the \code{Indices} slot of \code{knn} to use zero-indexing, set parameter \code{zero_index} to \code{TRUE}.
#' To include the zero-th neighbour (self) for each point in \code{knn}, set parameter \code{zeroth_neighbour} to \code{TRUE}.
#' If you want to limit the number of neighbours (to a value smaller than the original \code{k}), specify the new value using parameter \code{new_k}.
#'
#' @param knn list: \code{knn} object created either directly via \code{ComputekNNMatrix} or by including a \code{k}-NNG generation step in a benchmark pipeline and evaluating it
#' @param only_indices logical value: whether to only keep indices of nearest neighbours, discarding the distance matrix. Default value is \code{TRUE}
#' @param zero_index logical value: whether to use zero indexing for nearest-neighbour indices. Default value is \code{TRUE}
#' @param zeroth_neighbours logical value: whether to include the zero-th neighbour to each point (self). Default value is \code{TRUE}
#' @param new_k optional integer value: new \code{k} for \code{k}-NN. Default value is \code{NULL}
#'
#' @export
kNNGTweak <- function(
  knn,
  only_indices = TRUE,
  zero_index = TRUE,
  zeroth_neighbours = TRUE,
  new_k = NULL
) {
  if (!is.null(new_k)) {
    knn$Indices <- knn$Indices[, seq_len(new_k)]
    knn$Distances <- knn$Distances[, seq_len(new_k)]
  }
  if (zeroth_neighbours) {
    knn$Indices <- cbind(seq_len(nrow(knn$Indices)), knn$Indices)
    knn$Distances <- cbind(0, knn$Distances)
  }
  if (zero_index)
    knn$Indices <- knn$Indices - 1
  if (only_indices)
    knn <- knn$Indices
  knn
}

#' Load projection and clustering tool wrappers
#'
#' (Re-)loads wrappers defined in \code{inst/extdata/wrappers_projection} and \code{inst/extdata/wrappers_clustering} in the \code{SingleBench} package directory.
#' This can be used if you changed any of them or you removed the wrapper functions from your global namespace.
#'
#' @export
LoadWrappers <- function() {
  source(file.path(system.file(package = 'SingleBench'), 'extdata', 'load_wrappers.R'))
}

#' Generate a projection or clustering tool wrapper
#'
#' This function lets you create wrappers of projection or clustering tools.
#' Then, you can include them in benchmark pipelines.
#' 
#' # Basic components of a tool wrapper
#' 
#' To create a wrapper, you need to specify a handful of components (as arguments to \code{WrapTool}).
#' \code{name} is a unique string identifier. This is also included in the name of the wrapper (for example, \code{FlowSOM} will have \code{wrapper.clustering.FlowSOM}).
#' \code{type} specifies whether it is a projection tool (for dimension reduction or denoising) or clustering tool.
#' The string vector \code{r_packages} specifies names all required \code{R} packages and \code{python_modules} specifies names of required \code{Python} modules (that will be accessed via \code{reticulate}: the \code{R}/\code{Python} interface).
#' 
#' # Modelling functions
#' 
#' Modelling functions are the ones that do the work: transform input data.
#' At least one of them (\code{fun.build_model}) needs to be specified.
#' 
#' \code{fun.build_model.single_input} takes a single coordinate matrix of data and returns a model.
#' The model is an object from which the desired result (projection coordinate matrix or vector of cluster indices per data point) can be extracted.
#' \code{fun.build_model.batch_input}, instead, takes a list of multiple coordinate matrices (one per sample) as input and returns a model.
#'
#' If the tool does not distinguish between a single input matrix and multiple input matrices (it would just concatenate the inputs and apply \code{fun.build_model.single_input}), \code{fun.build_model.batch_input} can be left unspecified and it will be auto-generated.
#' In that case, you can  specify the function summarily as \code{fun.build_model}.
#'
#' \code{fun.extract} is a function that takes a model object (generated by \code{fun.build_model...}) as input and extracts results of the model applied to the original input data.
#' \code{fun.apply_model.single_input} takes a model object and a *new* coordinate matrix as input.
#' It returns the result of applying the *previously* trained model on *new* data.
#' \code{fun.apply_model.batch_input} takes a list of coordinate matrices as input and applies the model to new data.
#' 
#' Results of the \code{...batch_input} functions should not be split up into lists according to the sizes of the original inputs: they always return a single coordinate matrix or cluster vector (the splitting per sample is implemented automatically).
#'
#' ## Minimal function signatures
#' 
#' The minimal signature of a \code{fun.build_model...} function is \code{function(input)}.
#' Other arguments, with their default values, can (and should) be included: that way, changes in other parameters can be tested.
#'
#' For example, a simple signature of a \code{fun.build_model...} function for the dimension-reduction tool \code{t-SNE} might be \code{function(input, latent_dim = 2, perplexity = 2)}, allowing the user to alter target dimensionality or the perplexity parameter.
#'
#' Signatures of the other modelling functions are fixed.
#' For \code{fun.extract} it is \code{function(model)} and for \code{fun.apply_model...} it is \code{function(model, input)}.
#'
#' ## Additional inputs to model-building functions
#'
#' If a clustering tool uses the original high-dimensional expression data as well as a projection (generated in the previous step by some projection method), then include the parameter \code{expression} in your function signature and set parameter \code{use_original_expression_matrix} to \code{TRUE}.
#' \code{expression} is either a single matrix or a list of matrices, much like \code{input}.
#' \code{input}, then, will be the output of the preceding projection tool in that given sub-pipeline.
#' 
#' If your tool uses a \code{k}-nearest_neighbour graph (k-NNG), you are encouraged to always use one that was computed at the beginning of your pipeline evaluation.
#' (The \code{k}-NNG will be created if \code{SingleBench} knows it will run one or more tool that need it.)
#' To do this, set \code{use_knn_graph} to \code{TRUE} and add the argument \code{knn} to the signature of your model-building functions.
#' \code{knn} will then be a list of two names matrices: \code{Indices} for indices of nearest neighbours (row-wise) and \code{Distances} for distances to those neighbours.
#' 
#' **Warning**: the entries in \code{Indices} will be \code{1}-indexed and the matrices do not contain a column for the 'zero-th' neighbour (for each point, the zero-th neighbour is itself).
#' To modify the \code{knn} object (switch to 0-indexing or include zero-th neighobur), use the convertor \code{kNNTweak} inside your model-building function.
#' For instance, to convert \code{knn} to only a matrix of indices that does include zero-th neighbours, is 1-indexed and \code{k} is lowered from its original value to \code{30}, use: \code{knn <- kNNGTweak(knn, only_indices = TRUE, zero_index = TRUE, zeroth_neighbours = TRUE, new_k = 30)}.
#'
#' ## *n*-parameters
#' 
#' Most tools can accept custom numeric parameters.
#' Any one of the arguments to a model-building function can be chosen as the *n*-parameter by the user: then, \code{SingleBench} can do parameter sweeps over different values of these parameters.
#' Dimension-reduction tools, if possible, should have a parameter \code{latent_dim} for iterating over latent-space dimensionality.
#' Clustering tools, if possible, should have a parameter \code{n_clusters} for iterating over target cluster count.
#' If there is an option to determine number of clusters automatically, it might be a good idea to use \code{n_clusters = 0} for this.
#'
#' For methods that are made to run on multiple CPU cores, set \code{prevent_parallel_execution} to \code{TRUE} (otherwise, \code{SingleBench} may attempt to run them in parallel if the user wants repeated runs for stability analysis).
#'
#' @param name string: name of tool
#' @param type string: type of tool type (either '\code{projection}' or '\code{clustering}')
#' @param r_packages string vector: names of all \code{R} packages needed by the modelling functions
#' @param python_modules optional string vector: names of \code{Python} modules needed by the modelling function (called via \code{reticulate}). Default value is \code{NULL}
#' @param fun.build_model.single_input optional function: modelling function which accepts a single coordinate matrix as input data. Minimal signature \code{function(input, latent_dim)} or \code{function(input, n_clusters)}
#' @param fun.build_model.batch_input optional function: modelling function which accepts a list of coordinate matrices as input data. Minimal signature \code{function(input, latent_dim)} or \code{function(input, n_clusters)}. Default value is \code{NULL}
#' @param fun.build_model optional function: different parameter name for \code{fun.build_model.single_input}, if \code{fun.build_model.batch_input} is left as \code{NULL}
#' @param fun.extract optional function: modelling function which accepts a model generated by \code{fun.build_model.single_input} or \code{fun.build_model.batch_input} as input. Signature \code{function(model)}. If unspecified, the \code{model} object itself is taken as result
#' @param fun.apply_model.single_input optional function: modelling function which accepts a model generated by \code{fun.build_model.single_input} or \code{fun.build_model.batch_input} and new coordinate matrix as input. Signature \code{function(model, input)}. Default value is \code{NULL}
#' @param fun.apply_model.batch_input optional function: modelling function which accepts a model generated by \code{fun.build_model.single_input} or \code{fun.build_model.batch_input} and a new list coordinate matrices as input. Signature \code{function(model, input)}. Default value is \code{NULL}
#' @param fun.apply_model optional function: different parameter name for \code{fun.apply_model.single_input}, if \code{fun.apply_model.batch_input} is left as \code{NULL}
#' @param prevent_parallel_execution logical: whether running the tool in parallel on multiple CPU cores should be prevented. Default value is \code{TRUE}
#' @param use_python logical: whether the tool uses \code{Python} via \code{reticulate}. This is automatically set to \code{TRUE} if any \code{Python} modules are required. Otherwise, default value is \code{FALSE}
#' @param use_original_expression_matrix logical: whether the tool uses original expression matrix apart from the output of the preceding dimension-reduction tool. Default value is \code{FALSE}
#' @param use_knn_graph logical: whether the tool uses a \code{k}-nearest-neighbour graph of the input data. Default value is \code{FALSE}
#'
#' @return
#' This function returns a wrapper function that can be used in constructing a benchmark pipeline using \code{Fix}, \code{Module} and \code{Subpipeline}.
#'
#' @export
WrapTool <- function(
  name,
  type,
  r_packages                      = NULL,
  python_modules                  = NULL,
  fun.build_model.single_input    = NULL,
  fun.build_model.batch_input     = NULL,
  fun.build_model                 = NULL,
  fun.extract                     = function(model) model,
  fun.apply_model.single_input    = NULL,
  fun.apply_model.batch_input     = NULL,
  fun.apply_model                 = NULL,
  prevent_parallel_execution      = TRUE,
  use_python                      = !is.null(python_modules) && length(python_modules) > 0,
  use_original_expression_matrix  = FALSE,
  use_knn_graph                   = FALSE
) {

  .WrapTool.ValidityChecks(environment())
  
  r_packages_missing <- r_packages[!r_packages %in% installed.packages()[, 1]]
  if (length(r_packages_missing) > 0) {
    message(paste0('! -> ', type, ' wrapper "', name, '" is missing R packages: ', paste(r_packages_missing, sep = ', ')))
  }
  if (length(python_modules) > 0) {
    python_modules_missing <- which(!purrr::map_lgl(python_modules, reticulate::py_module_available))
    if (length(python_modules_missing) > 0) {
      message(paste0('! -> ', type, ' wrapper "', name, '" is missing Python modules: ', paste(unlist(python_modules[python_modules_missing]), collapse = ', ')))
    }
    prevent_parallel_execution <- TRUE
  }
  
  ## Outward-facing model-training function
  train <- function(
    input,              # input matrix
    expression = NULL,  # full expression matrix (original)
    knn = NULL,         # k-nearest-neighbour graph (optional)
    params,             # additional parameters
    out.systime = NULL  # out-variable to assign elapsed time
  ) {
    
    if (is.null(fun.build_model.single_input))
      fun.build_model.single_input <- fun.build_model
    
    ## Concatenate batch input if necessary
    input_is_batch                     <- is.list(input)
    batch_input_is_handled_differently <- !is.null(fun.build_model.batch_input)
    if (input_is_batch && !batch_input_is_handled_differently)
      input <- do.call(rbind, input)
    
    ## Compile list of parameters to be passed to model-building function
    parameter_list <- list(input)
    names(parameter_list) <- 'input'
    if (use_original_expression_matrix) {
      parameter_list <- c(parameter_list, list(expression))
      names(parameter_list)[length(parameter_list)] <- 'expression'
    }
    if (use_knn_graph) {
      parameter_list <- c(parameter_list, list(knn))
      names(parameter_list)[length(parameter_list)] <- 'knn'
    }
    parameter_list <- c(parameter_list, params)
    
    ## Deploy model-building function and measure elapsed time
    systime <-
      if (!input_is_batch || (input_is_batch && !batch_input_is_handled_differently))
        system.time(
          res <- do.call(fun.build_model.single_input, parameter_list)
        )
      else # if (input_is_batched && batch_input_is_handled_differently)
        system.time(
          res <- do.call(fun.build_model.batch_input, parameter_list)
        )
    
    ## Assign elapsed time value to an out-variable (side-effect)
    if (!is.null(out.systime))
      eval.parent(substitute(out.systime <- systime))
    res
  }
  
  ## Outward-facing function to extract results from a trained model
  extract <- fun.extract
  
  ## Outward-facing function to map new data onto an existing model
  map <-
    if (is.null(fun.apply_model) & is.null(fun.apply_model.single_input))
      NULL
    else
      function(model, input, out.systime = NULL) {
        if (is.null(fun.apply_model.single_input))
          fun.apply_model.single_input <- fun.apply_model
        
        ## Concatenate batch input if necessary
        input_is_batch <- is.list(input)
        batch_input_is_handled_differently <- !is.null(fun.apply_model.batch_input)
        if (input_is_batch && !batch_input_is_handled_differently)
          input <- do.call(rbind, input)
        
        ## Deploy model-applying function and measure elapsed time
        systime <-
          if (!input_is_batch || (input_is_batch && !batch_input_is_handled_differently))
            system.time(
              res <- fun.apply_model.single_input(model, input)
            )
          else # if (input_is_batched && batch_input_is_handled_differently)
            system.time(
              res <- fun.apply_model.batch_input(model, input)
            )
        ## Assign elapsed time value to an out-variable (side-effect)
        if (!is.null(out.systime))
          eval.parent(substitute(out.systime <- systime))
        res
      }
  wrapper <- list(
    name           = name,
    type           = type,
    r_packages     = r_packages,
    uses_python    = use_python,
    python_modules = python_modules,
    
    args_to_train     = formalArgs(if (!is.null(fun.build_model.single_input)) fun.build_model.single_input else fun.build_model),
    defaults_to_train = formals(if (!is.null(fun.build_model.single_input)) fun.build_model.single_input else fun.build_model),
    
    train          = train,
    extract        = extract,
    map            = map,
    
    prevent_parallel_execution      = prevent_parallel_execution,
    uses_original_expression_matrix = use_original_expression_matrix,
    uses_knn_graph                  = use_knn_graph
  )
  class(wrapper) <-
    c(
      'BenchmarkToolWrapper',
      if (type == 'projection')
        'BenchmarkToolWrapper_Projection'
      else if (type == 'clustering')
        'BenchmarkToolWrapper_Clustering'
    )
  message(paste0('-> ', name))
  wrapper
}

CloneProjectionTool <- function(
  idx.subpipeline
) {
  wrapper <- list(name = 'ClonePrevious', idx = idx.subpipeline)
  class(wrapper) <- c('BenchmarkToolWrapper', 'BenchmarkToolWrapper_Projection')
  wrapper
}

CloneClusteringTool <- function(
  idx.subpipeline
) {
  wrapper <- list(name = 'ClonePrevious', idx = idx.subpipeline)
  class(wrapper) <- c('BenchmarkToolWrapper', 'BenchmarkToolWrapper_Clustering')
  wrapper
}
