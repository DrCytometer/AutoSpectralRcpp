#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// Compiled counterpart of assign.af.joint.cov.l2() (assign_af_joint_cov_l2.R).
// Scores each cell/variant pair by the product of two proportional squared
// errors -- a covariance-weighted fluorophore error and a raw-space
// residual error -- both expanded as quadratics in the per-cell, per-variant
// AF abundance k so every variant can be evaluated with a fixed amount of
// per-cell work rather than a materialised (cell x variant) matrix
// subtraction per variant, as the R version does. See assign_af_joint_cov_l2.R
// for the derivation; this file only translates the arithmetic to per-cell
// scalar loops for OpenMP parallelism.
//
// Assumes `spectra` has already had any "AF" row removed by the caller,
// matching the convention of assign_af_fluor.cpp.

// [[Rcpp::export]]
Rcpp::IntegerVector assign_af_joint_cov_l2_cpp(const arma::mat& raw_data,
                                               const arma::mat& spectra,
                                               const arma::mat& af_spectra,
                                               int n_threads = 4) {

  // Setup dimensions
  int n_cells    = raw_data.n_rows;
  int n_channels = raw_data.n_cols;
  int n_fluors   = spectra.n_rows;
  int n_af       = af_spectra.n_rows;

  // ---------------------------------------------------------------------
  // Global precomputation (single-threaded, mirrors assign.af.joint.cov.l2)
  // ---------------------------------------------------------------------

  // Pseudoinverse: (S*S')^-1 * S
  mat P    = solve(spectra * spectra.t(), spectra);   // n_fluors x n_channels
  mat S_t  = spectra.t();                              // n_channels x n_fluors
  mat AF_t = af_spectra.t();                            // n_channels x n_af

  // v.library: how much each AF variant looks like each fluorophore
  mat v_library = P * AF_t;                             // n_fluors x n_af

  // r.library: residual AF signal, orthogonal to the fluorophore span
  mat r_library = AF_t - (S_t * v_library);              // n_channels x n_af

  // Raw (unfloored) self-dot of each variant's residual. Used both as the
  // true curvature term in the residual error expansion and, floored, as
  // the denominator of the k estimate (identifiability guard: a variant
  // lying almost inside the fluorophore span has a vanishing residual
  // direction, so its k would otherwise explode).
  vec r_dots(n_af);
  for (int j = 0; j < n_af; ++j) {
    r_dots[j] = dot(r_library.col(j), r_library.col(j));
  }
  double max_r_dot = 0.0;
  for (int j = 0; j < n_af; ++j) {
    if (r_dots[j] > max_r_dot) max_r_dot = r_dots[j];
  }
  double floor_val = 0.01 * std::max(max_r_dot, 1e-10);
  vec denom_floor(n_af);
  for (int j = 0; j < n_af; ++j) {
    denom_floor[j] = std::max(r_dots[j], floor_val);
  }

  // Covariance-based fluorophore error weights: propagate AF spectral
  // covariance into fluorophore space via the unmixing matrix, then take
  // the per-channel SD-scale weight from the diagonal.
  mat af_cov    = arma::cov(af_spectra);                 // n_channels x n_channels
  mat fluor_cov = P * af_cov * P.t();                     // n_fluors x n_fluors
  vec af_error_weights(n_fluors);
  for (int f = 0; f < n_fluors; ++f) {
    af_error_weights[f] = std::sqrt(std::abs(fluor_cov(f, f)));
  }

  // c.fluor[j] = sum_f w_f * v_library(f,j)^2 -- curvature term of the
  // weighted squared fluorophore error, precomputed once per variant
  vec c_fluor(n_af, fill::zeros);
  for (int j = 0; j < n_af; ++j) {
    double acc = 0.0;
    for (int f = 0; f < n_fluors; ++f) {
      double v = v_library(f, j);
      acc += af_error_weights[f] * v * v;
    }
    c_fluor[j] = acc;
  }

  // Transpose raw_data to make per-cell access contiguous (column-major)
  mat Y_t = raw_data.t();   // n_channels x n_cells

  // Output: best AF index per cell
  Rcpp::IntegerVector best_indices(n_cells);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  // ---------------------------------------------------------------------
  // Parallel per-cell loop
  // ---------------------------------------------------------------------
#pragma omp parallel
{
  // Thread-local buffers, declared `static thread_local` so each OS thread
  // allocates them once and reuses the memory across every call to this
  // function for the life of the R session, not just across cells within a
  // single call. Explicitly resized every call since panel shape can
  // differ between calls; set_size() is a no-op when the size already
  // matches.
  static thread_local vec unmixed_i;
  static thread_local vec unmixed_nonneg;
  static thread_local vec resid_initial;

  unmixed_i.set_size(n_fluors);
  unmixed_nonneg.set_size(n_fluors);
  resid_initial.set_size(n_channels);

#pragma omp for schedule(static)
  for (int i = 0; i < n_cells; ++i) {
    const double* y_ptr = Y_t.colptr(i);

    // Step A: initial (no-AF) unmix for this cell, raw and non-negative-
    // clipped copies both needed downstream (base.e.fluor uses the raw
    // values; the raw-space residual is built from the clipped ones).
    for (int f = 0; f < n_fluors; ++f) {
      double sum = 0.0;
      for (int c = 0; c < n_channels; ++c) {
        sum += P(f, c) * y_ptr[c];
      }
      unmixed_i[f]      = sum;
      unmixed_nonneg[f] = std::max(sum, 0.0);
    }

    // Step B: raw-space residual against the clipped unmix
    for (int c = 0; c < n_channels; ++c) {
      double proj = 0.0;
      for (int f = 0; f < n_fluors; ++f) {
        proj += S_t(c, f) * unmixed_nonneg[f];
      }
      resid_initial[c] = y_ptr[c] - proj;
    }

    // Step C: baseline (variant-free) weighted squared errors
    double base_e_fluor = 1e-6;
    for (int f = 0; f < n_fluors; ++f) {
      base_e_fluor += af_error_weights[f] * unmixed_i[f] * unmixed_i[f];
    }
    double base_e_resid = 1e-6;
    for (int c = 0; c < n_channels; ++c) {
      base_e_resid += resid_initial[c] * resid_initial[c];
    }

    double min_score  = datum::inf;
    int    best_af_val = 0;

    // Step D: iterate through AF variants, expanding both quadratics
    for (int j = 0; j < n_af; ++j) {

      // k_ij: estimated AF intensity for this cell/variant
      double numerator = 0.0;
      const double* r_ptr = r_library.colptr(j);
      for (int c = 0; c < n_channels; ++c) {
        numerator += y_ptr[c] * r_ptr[c];
      }
      double k_ij = numerator / denom_floor[j];
      if (k_ij < 0.0) k_ij = 0.0;

      // e.fluor_j = base.e.fluor - 2 k (w.unmixed . v_j) + k^2 c.fluor_j
      double cross_fluor = 0.0;
      const double* v_ptr = v_library.colptr(j);
      for (int f = 0; f < n_fluors; ++f) {
        cross_fluor += af_error_weights[f] * unmixed_i[f] * v_ptr[f];
      }
      double e_fluor = base_e_fluor - 2.0 * k_ij * cross_fluor + k_ij * k_ij * c_fluor[j];
      if (e_fluor < 0.0) e_fluor = 0.0;

      // e.resid_j = base.e.resid - 2 k (resid . r_j) + k^2 ||r_j||^2
      double cross_resid = 0.0;
      for (int c = 0; c < n_channels; ++c) {
        cross_resid += resid_initial[c] * r_ptr[c];
      }
      double e_resid = base_e_resid - 2.0 * k_ij * cross_resid + k_ij * k_ij * r_dots[j];
      if (e_resid < 0.0) e_resid = 0.0;

      double score = (e_fluor / base_e_fluor) * (e_resid / base_e_resid);

      if (score < min_score) {
        min_score   = score;
        best_af_val = j + 1;   // R-style 1-indexing
      }
    }
    best_indices[i] = best_af_val;
  }
}

  return best_indices;
}
