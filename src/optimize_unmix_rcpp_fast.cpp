// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// --- OLS unmix ---
inline arma::rowvec unmix_ols(const arma::rowvec& raw_row,
                              const arma::mat& spectra) {
  arma::mat M = spectra.t(); // detectors × fluorophores
  arma::mat unmixing_matrix = arma::solve(M.t() * M, M.t());
  return raw_row * unmixing_matrix.t();
}

// --- WLS unmix ---
inline arma::rowvec unmix_wls(const arma::rowvec& raw_row,
                              const arma::mat& spectra,
                              const arma::vec& weights) {
  arma::mat M = spectra.t();                       // detectors × fluorophores
  arma::mat MW = M.each_col() % weights;           // weighted M
  arma::rowvec raw_w = raw_row % weights.t();      // weighted raw_row

  arma::mat A = M.t() * MW;                        // (fluor × fluor)
  arma::rowvec b = raw_w * M;                      // (1 × fluor)

  arma::rowvec coef = arma::solve(A, b.t()).t();
  return coef;
}

// [[Rcpp::export]]
arma::mat optimize_unmix_rcpp_fast(const arma::mat& remaining_raw,
                                   const arma::mat& unmixed,
                                   const arma::mat& spectra,
                                   const arma::vec& pos_thresholds,
                                   const arma::uvec& optimize_idx_r,
                                   const std::vector<arma::mat>& variantsList,
                                   const arma::vec& weights,
                                   const bool weighted = false,
                                   const int nthreads = 1) {

  arma::uword n_cells = remaining_raw.n_rows;
  arma::uword n_detectors = remaining_raw.n_cols;
  arma::uword n_fluors = spectra.n_rows;

  arma::mat result = unmixed;

  // Validate weights
  arma::vec w = (weights.n_elem == 0) ? arma::ones<arma::vec>(n_detectors) : weights;
  if (w.n_elem != n_detectors) Rcpp::stop("weights length mismatch");

  // Determine actual number of threads
#ifdef _OPENMP
  int nthreads_actual = omp_get_max_threads();
#else
  int nthreads_actual = 1;
#endif

  // Preallocate temporary storage per thread
  std::vector<arma::rowvec> tmp_raw_row(nthreads_actual);
  std::vector<arma::rowvec> tmp_cell_unmixed(nthreads_actual);
  std::vector<arma::rowvec> tmp_fitted(nthreads_actual);
  std::vector<arma::rowvec> tmp_resid(nthreads_actual);
  std::vector<arma::rowvec> tmp_delta(nthreads_actual);
  std::vector<arma::rowvec> tmp_resid_candidate(nthreads_actual);
  std::vector<arma::mat> tmp_spectra_final(nthreads_actual, spectra); // per-thread copy

  // Select unmix function pointer once
  auto unmix_fn = [&](const arma::rowvec& raw_row,
                      const arma::mat& S) -> arma::rowvec {
                        if (weighted) return unmix_wls(raw_row, S, w);
                        else return unmix_ols(raw_row, S);
                      };

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
#endif
  for (arma::uword cell = 0; cell < n_cells; cell++) {

#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

    // Copy raw row and initial unmixed
    tmp_raw_row[tid] = remaining_raw.row(cell);
    tmp_cell_unmixed[tid] = result.row(cell);

    // Identify positive fluorophores
    std::vector<arma::uword> pos_vec;
    for (arma::uword f = 0; f < n_fluors; f++) {
      if (tmp_cell_unmixed[tid](f) >= pos_thresholds(f)) pos_vec.push_back(f);
    }
    if (pos_vec.empty()) continue;

    arma::uvec pos_idx(pos_vec.size());
    for (arma::uword i = 0; i < pos_vec.size(); i++) pos_idx(i) = pos_vec[i];

    // Precompute cell spectra for positive fluorophores
    arma::mat cell_spectra(pos_idx.n_elem, n_detectors);
    for (arma::uword i = 0; i < pos_idx.n_elem; i++)
      cell_spectra.row(i) = spectra.row(pos_idx(i));

    // Initial unmix
    arma::rowvec unmixed_cell = unmix_fn(tmp_raw_row[tid], cell_spectra);

    // Compute initial residual
    tmp_fitted[tid] = unmixed_cell * cell_spectra;
    tmp_resid[tid] = tmp_raw_row[tid] - tmp_fitted[tid];
    double error_final = arma::accu(arma::abs(tmp_resid[tid]));

    // Variant optimization
    for (arma::uword vi = 0; vi < pos_idx.n_elem; vi++) {
      arma::uword f_idx = pos_idx(vi);
      if (f_idx >= variantsList.size()) continue;
      const arma::mat& fl_variants = variantsList[f_idx];
      if (fl_variants.n_rows == 0) continue;

      for (arma::uword v = 0; v < fl_variants.n_rows; v++) {
        tmp_delta[tid] = fl_variants.row(v) - cell_spectra.row(vi);
        tmp_delta[tid] *= unmixed_cell(vi);
        tmp_resid_candidate[tid] = tmp_resid[tid] - tmp_delta[tid];

        double err = arma::accu(arma::abs(tmp_resid_candidate[tid]));
        if (err < error_final) {
          cell_spectra.row(vi) = fl_variants.row(v);
          tmp_resid[tid] = tmp_resid_candidate[tid];
          error_final = err;
        }
      }
    }

    // Copy optimized positive fluorophores into spectra_final (per-thread copy)
    tmp_spectra_final[tid] = spectra; // reset base
    for (arma::uword vi = 0; vi < pos_idx.n_elem; vi++)
      tmp_spectra_final[tid].row(pos_idx(vi)) = cell_spectra.row(vi);

    // Final unmix
    result.row(cell) = unmix_fn(tmp_raw_row[tid], tmp_spectra_final[tid]);
  }

  return result;
}
