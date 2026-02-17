#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::List unmix_af_residuals(const arma::mat& raw_data,
                              const arma::mat& fluor_spectra,
                              const arma::mat& af_spectra,
                              int n_threads = 4) {
  
  int n_cells = raw_data.n_rows;
  int n_channels = raw_data.n_cols;
  int n_fluors = fluor_spectra.n_rows;
  int n_af = af_spectra.n_rows;
  
  // global pre-calculations
  mat P = solve(fluor_spectra * fluor_spectra.t(), fluor_spectra);
  
  // Transpose AF spectra so each AF variant is a contiguous column
  mat AF_t = af_spectra.t();
  
  mat p_af_library(n_fluors, n_af);
  mat r_af_library(n_channels, n_af);
  vec dot_r_af_af(n_af);
  
  for (int j = 0; j < n_af; j++) {
    vec af_j = AF_t.col(j);
    vec p_j = P * af_j;
    p_af_library.col(j) = p_j;
    vec r_j = af_j - (fluor_spectra.t() * p_j);
    r_af_library.col(j) = r_j;
    dot_r_af_af[j] = dot(r_j, r_j);
  }
  
  // Transpose raw_data for contiguous cell access (channels x cells)
  mat Y_t = raw_data.t();
  
  // Results containers
  mat unmixed(n_cells, n_fluors + 1);
  mat fitted_af(n_cells, n_channels);
  Rcpp::IntegerVector af_indices(n_cells);
  
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n_cells; ++i) {
    const double* y_ptr = Y_t.colptr(i);
    
    // --- Initial Fluorophore-only unmix ---
    vec init_fluor(n_fluors);
    for(int f = 0; f < n_fluors; ++f) {
      double sum = 0;
      for(int c = 0; c < n_channels; ++c) {
        sum += P(f, c) * y_ptr[c];
      }
      init_fluor[f] = sum;
    }
    
    // --- Calculate Cell Residual ---
    vec res(n_channels);
    for(int c = 0; c < n_channels; ++c) {
      double fit = 0;
      for(int f = 0; f < n_fluors; ++f) {
        fit += fluor_spectra(f, c) * init_fluor[f];
      }
      res[c] = y_ptr[c] - fit;
    }
    
    // --- Score AF variants ---
    double max_score = -datum::inf;
    int best_idx = 0;
    for(int a = 0; a < n_af; ++a) {
      double score = 0;
      const double* af_ptr = AF_t.colptr(a);
      for(int c = 0; c < n_channels; ++c) {
        score += res[c] * af_ptr[c];
      }
      if(score > max_score) {
        max_score = score;
        best_idx = a;
      }
    }
    
    // --- Exact Correction ---
    double dot_res_raf = 0;
    const double* r_af_ptr = r_af_library.colptr(best_idx);
    for(int c = 0; c < n_channels; ++c) {
      dot_res_raf += res[c] * r_af_ptr[c];
    }
    
    double k_af = dot_res_raf / dot_r_af_af[best_idx];
    
    // Store unmixed fluors
    const double* p_af_ptr = p_af_library.colptr(best_idx);
    for(int f = 0; f < n_fluors; ++f) {
      unmixed(i, f) = init_fluor[f] - (k_af * p_af_ptr[f]);
    }
    
    unmixed(i, n_fluors) = k_af;
    af_indices[i] = best_idx + 1;
    
    // fitted_af = intensity * actual_spectrum
    const double* best_af_ptr = AF_t.colptr(best_idx);
    for(int c = 0; c < n_channels; ++c) {
      fitted_af(i, c) = k_af * best_af_ptr[c];
    }
  }
  
  return Rcpp::List::create(
    Rcpp::_["unmixed"] = unmixed, 
    Rcpp::_["af.idx"] = af_indices,
    Rcpp::_["fitted.af"] = fitted_af
  );
}