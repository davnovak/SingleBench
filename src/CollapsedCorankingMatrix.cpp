#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
# define NUM_THREADS(N) ((N) >= 0 ? (N) : omp_get_num_procs() + (N) + 1)

// [[Rcpp::export]]
NumericMatrix
  CollapsedCorankingMatrix(
    NumericMatrix idcs1,
    NumericMatrix idcs2
  ) {
    const size_t k = idcs1.ncol();
    const size_t n = idcs1.nrow();
    
    NumericMatrix M(k + 1, k + 1);
    
#pragma omp parallel for
    for (size_t i = 0; i < k; ++i) {
      int rowsum = 0;
      for (size_t j = 0; j < k; ++j) {
        int val = 0;
        for (size_t l = 0; l < n; ++l) {
          if (idcs1(l, i) == idcs2(l, j)) {
            val += 1;
          }
        } // l
        M(i, j) = val;
        rowsum += val;
      } // j
      M(i, k) = n - rowsum;
    } // i
    
#pragma omp parallel for
    for (size_t j = 0; j < k; ++j) {
      int colsum = 0;
      for (size_t i = 0; i < k; ++i) {
        colsum += M(i, j);
      } // i
      M(k, j) = n - colsum;
    } // j
    
    return M;
  }
