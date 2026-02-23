#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalMatrix find_local_maxima_cpp(NumericMatrix z, int neigh_size) {
  int nx = z.nrow();
  int ny = z.ncol();
  LogicalMatrix max_bool(nx, ny);
  
  // ignore edges
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      double current_val = z(i, j);
      bool is_max = true;
      
      // Define search window, clamped to matrix boundaries
      int i_start = std::max(0, i - neigh_size);
      int i_end   = std::min(nx - 1, i + neigh_size);
      int j_start = std::max(0, j - neigh_size);
      int j_end   = std::min(ny - 1, j + neigh_size);
      
      for (int wi = i_start; wi <= i_end; ++wi) {
        for (int wj = j_start; wj <= j_end; ++wj) {
          if (z(wi, wj) > current_val) {
            is_max = false;
            break;
          }
        }
        if (!is_max) break;
      }
      max_bool(i, j) = is_max;
    }
  }
  return max_bool;
}