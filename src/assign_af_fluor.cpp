#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::IntegerVector assign_af_fluor(const arma::mat& raw_data, 
                                    const arma::mat& spectra, 
                                    const arma::mat& af_spectra, 
                                    int n_threads = 4) {
  
  // Setup dimensions
  int n_cells = raw_data.n_rows;
  int n_channels = raw_data.n_cols;
  int n_fluors = spectra.n_rows;
  int n_af = af_spectra.n_rows;
  
  // Pseudoinverse: (S*S')^-1 * S
  mat P = solve(spectra * spectra.t(), spectra);
  mat S_t = spectra.t();
  mat AF_t = af_spectra.t(); // channels x af
  
  // v.library: (n_fluors x n_af) - how much each AF looks like fluors
  mat v_library = P * AF_t;
  
  // r.library: (n_channels x n_af) - the residual AF
  mat r_library = AF_t - (S_t * v_library);
  
  // Denominator for k calculation: colSums(r.library^2)
  vec r_dots(n_af);
  for(int j = 0; j < n_af; ++j) {
    r_dots[j] = dot(r_library.col(j), r_library.col(j));
  }
  
  // Transpose raw_data to make cell access contiguous (column-major memory)
  mat Y_t = raw_data.t();
  
  // Output: Best AF index per cell
  Rcpp::IntegerVector best_indices(n_cells);
  
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  
  // Parallel Cell Loop
#pragma omp parallel
{
  // Thread-local buffers to avoid repeated allocations
  vec init_fluor(n_fluors);
  
#pragma omp for schedule(static)
  for (int i = 0; i < n_cells; ++i) {
    const double* y_ptr = Y_t.colptr(i);
    
    // Step A: Initial unmix for this cell (P * y)
    for(int f = 0; f < n_fluors; ++f) {
      double sum = 0;
      for(int c = 0; c < n_channels; ++c) {
        sum += P(f, c) * y_ptr[c];
      }
      init_fluor[f] = sum;
    }
    
    double min_l1_error = datum::inf;
    int best_af_val = 0;
    
    // Step B: Iterate through AF variants
    for (int j = 0; j < n_af; ++j) {
      // Calculate k_i: (y . r_j) / dot(r_j, r_j)
      double numerator = 0;
      const double* r_ptr = r_library.colptr(j);
      for(int c = 0; c < n_channels; ++c) {
        numerator += y_ptr[c] * r_ptr[c];
      }
      double ki = numerator / r_dots[j];
      
      // Calculate L1 error: sum(abs(init_fluor - ki * v_j))
      double current_l1 = 0;
      const double* v_ptr = v_library.colptr(j);
      for(int f = 0; f < n_fluors; ++f) {
        current_l1 += std::abs(init_fluor[f] - (ki * v_ptr[f]));
      }
      
      if (current_l1 < min_l1_error) {
        min_l1_error = current_l1;
        best_af_val = j + 1; // R-style 1-indexing
      }
    }
    best_indices[i] = best_af_val;
  }
}

return best_indices;
}