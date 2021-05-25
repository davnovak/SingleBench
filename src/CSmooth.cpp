#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
# define NUM_THREADS(N) ((N) >= 0 ? (N) : omp_get_num_procs() + (N) + 1)

// [[Rcpp::export]]
NumericMatrix
  pSmoothSingleIter(
    NumericMatrix    exprs,   // coordinate matrix
    NumericMatrix    idcs,    // matrix of nearest neighbour indices
    NumericMatrix    weights, // matrix of weights (based on initial distances)
    const double     phi1,    // different-segment weight coefficient
    const double     phi2,    // same-segment weight coefficient
    const double     lambda,  // learning coefficient
    std::vector<int> segs,    // segment per data point
    const size_t     k        // nearest neighbour count
  ) {
    
    const size_t n = exprs.nrow();
    const size_t d = exprs.ncol();
    
    NumericMatrix res(n, d);
    
#pragma omp parallel for
    for (size_t irow = 0; irow < n; ++irow) {
      
      // Sum change vectors
      std::vector<double> nx(d);
      double total_phi = 0;
      for (size_t irank = 0; irank < k; ++irank) {
        size_t ineighb = idcs(irow, irank) - 1; // reindex
        double phi;
        if (segs[irow] == segs[ineighb]) {
          phi = phi2;
        } else {
          phi = phi1;
        }
        total_phi += phi;
        for (size_t idim = 0; idim < d; ++idim) {
          double delta = (exprs(ineighb, idim) - exprs(irow, idim)) * phi; 
          nx[idim] += delta;
        } // idim
      } // irank
      
      // Normalise and shift
      for (size_t idim = 0; idim < d; ++idim) {
        res(irow, idim) = exprs(irow, idim) + (nx[idim] / total_phi) * lambda * weights(irow, idim);
      }
    } // irow
    
    return res;
  }
