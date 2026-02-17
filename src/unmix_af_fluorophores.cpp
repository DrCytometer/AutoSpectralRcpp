#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::List unmix_af_fluorophores(const arma::mat& raw_data,
                                 const arma::mat& fluor_spectra,
                                 const arma::mat& af_spectra,
                                 int n_threads = 4) {
  
  int n_cells = raw_data.n_rows;
  int n_channels = raw_data.n_cols;
  int n_fluors = fluor_spectra.n_rows;
  int n_af = af_spectra.n_rows;
  
  // global pre-calculations and transpositions
  mat P = solve(fluor_spectra * fluor_spectra.t(), fluor_spectra);
  mat S_t = fluor_spectra.t();
  mat AF_t = af_spectra.t();
  
  // v_library: how much each AF looks like fluors (fluors x af)
  mat v_library = P * AF_t;
  
  // r_library: the residual part of AF (channels x af)
  mat r_library = AF_t - (S_t * v_library);
  
  // Denominators for k-calculation: colSums(r.library^2)
  vec r_dots(n_af);
  for(int j = 0; j < n_af; ++j) {
    r_dots[j] = dot(r_library.col(j), r_library.col(j));
  }
  
  // Transpose raw_data for contiguous cell access
  mat Y_t = raw_data.t();
  
  // Results containers
  mat unmixed_res(n_cells, n_fluors + 1);
  mat fitted_af(n_cells, n_channels);
  Rcpp::IntegerVector af_indices(n_cells);
  
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n_cells; ++i) {
    const double* y_ptr = Y_t.colptr(i);
    
    // --- Initial Unmix (Base Fluorophores) ---
    vec init_fluor(n_fluors);
    for(int f = 0; f < n_fluors; ++f) {
      double sum = 0;
      for(int c = 0; c < n_channels; ++c) {
        sum += P(f, c) * y_ptr[c];
      }
      init_fluor[f] = sum;
    }
    
    // --- Calculate Best AF based on Fluorophore Error ---
    double min_total_error = datum::inf;
    int best_idx = 0;
    double best_k = 0;
    
    for (int j = 0; j < n_af; ++j) {
      // Calculate k_i for this cell/AF variant
      // k = (raw_data %*% r_library) / r_dots
      double numerator = 0;
      const double* r_ptr = r_library.colptr(j);
      for(int c = 0; c < n_channels; ++c) {
        numerator += y_ptr[c] * r_ptr[c];
      }
      double ki = numerator / r_dots[j];
      
      // Calculate L1 Error: sum(|init_fluor - ki * v_library_j|)
      double current_error = 0;
      const double* v_ptr = v_library.colptr(j);
      for(int f = 0; f < n_fluors; ++f) {
        current_error += std::abs(init_fluor[f] - (ki * v_ptr[f]));
      }
      
      if (current_error < min_total_error) {
        min_total_error = current_error;
        best_idx = j;
        best_k = ki;
      }
    }
    
    // --- Finalize results for the best-fitted AF ---
    const double* best_v_ptr = v_library.colptr(best_idx);
    for (int f = 0; f < n_fluors; ++f) {
      unmixed_res(i, f) = init_fluor[f] - (best_k * best_v_ptr[f]);
    }
    
    unmixed_res(i, n_fluors) = best_k;
    af_indices[i] = best_idx + 1;
    
    // --- Reconstruct Fitted AF ---
    const double* af_spec_ptr = AF_t.colptr(best_idx);
    for (int c = 0; c < n_channels; ++c) {
      fitted_af(i, c) = best_k * af_spec_ptr[c];
    }
  }
  
  return Rcpp::List::create(
    Rcpp::_["unmixed"] = unmixed_res, 
    Rcpp::_["af.idx"] = af_indices,
    Rcpp::_["fitted.af"] = fitted_af
  );
}