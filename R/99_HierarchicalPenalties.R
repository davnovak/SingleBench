
#' Create a match/mismatch scoring matrix with hierarchical penalties
#'
#' Creates a matrix of scores for evaluating the prediction of cell population labels.
#' Each pair of labels gets assigned a score for match (if they are identical) or mismatch.
#' Rows of the matrix correspond to true labels and columns correspond to predicted labels.
#' 
#'  ## Hierarchical penalties
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
#' @param labels character vector: all possible manual labels
#' @param g_c numeric: constant generalisation penalty. Default value is \code{0}
#' @param g_a numeric: additive generalisation penalty. Default value is \code{0.2}
#' @param m_c numeric: constant misidentification penalty. Default value is \code{0.4}
#' @param m_a numeric: additive misidentification penalty. Default value is \code{0}
#'
#' @seealso
#' 
#' * **\code{GetPenaltyScoringMatrix}**: extracts a penalty scoring matric from a \code{Benchmark}-type object
#'
#' @export
CreatePenaltyScoringMatrix <- function(
  labels,
  g_c = 0,
  g_a = 0.2,
  m_c = 0.4,
  m_a = 0
) {
  N <- length(labels)
  m <- matrix(NA, ncol = N, nrow = N, dimnames = list(labels, labels))
  for (j in seq_len(N)) {
    for (i in seq_len(N)) {
      m[i, j] <- 1 - penalty(labels[i], labels[j], g_c, g_a, m_c, m_a)
    }
  }
  m
}

SavePenaltyScoringMatrix <- function(
  benchmark,
  obj
) {
  .h5writeNamedMatrix(benchmark, obj, '/Input/Annotation/ScoringMatrix')
}

penalty <- function(
  true,
  pred,
  g_c = 0,
  g_a = 0.2,
  m_c = 0.4,
  m_a = 0
) {
  if (true == pred)
    return(0)
  
  true <- strsplit(true, '/')[[1]][-1]
  n_true <- length(true)
  pred <- strsplit(pred, '/')[[1]][-1]
  n_pred <- length(pred)
  
  g_steps<- 0 # total steps toward root to get from pred to true
  m_steps <- 0 # total steps away from root to get from pred to true
  
  if (n_pred < n_true) {
    g_steps <- n_true - n_pred
  } else {
    m_steps <- n_true - n_pred
  }
  
  n_overlap <- min(n_true, n_pred)
  idx_clash <- which(true[1:n_overlap] != pred[1:n_overlap])[1]
  
  if (!is.na(idx_clash)) {
    extra_steps <- n_overlap - idx_clash + 1
    g_steps <- g_steps + extra_steps
    m_steps <- m_steps + extra_steps
  }
  
  g_c + g_a * g_steps + m_c + m_a * m_steps
}
