// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Compute unmixing matrix (OLS or WLS)
inline arma::mat compute_unmixing_matrix_general(const arma::mat& spectra,
                                                 const arma::vec& weights,
                                                 bool weighted) {
  arma::mat M = spectra.t();
  arma::uword D = M.n_rows;

  if (!weighted || weights.n_elem == 0) {
    arma::mat A = M.t() * M;
    arma::mat B = M.t();
    arma::mat X;
    bool ok = arma::solve(X, A, B);
    if (!ok) X = arma::pinv(A) * B;
    return X;
  } else {
    arma::mat MtW = M.t();
    for (arma::uword j = 0; j < D; ++j)
      MtW.col(j) *= weights(j);
    arma::mat A = MtW * M;
    arma::mat B = MtW;
    arma::mat X;
    bool ok = arma::solve(X, A, B);
    if (!ok) X = arma::pinv(A) * B;
    return X;
  }
}

// Safe “fast solve” for intermediate steps
inline arma::mat compute_unmixing_matrix_fast(const arma::mat& spectra,
                                              const arma::vec& weights,
                                              bool weighted) {
  arma::mat M = spectra.t();
  arma::uword D = M.n_rows;
  arma::mat A, B, X;

  if (!weighted || weights.n_elem == 0) {
    A = M.t() * M;
    B = M.t();
  } else {
    arma::mat MtW = M.t();
    for (arma::uword j = 0; j < D; ++j)
      MtW.col(j) *= weights(j);
    A = MtW * M;
    B = MtW;
  }

  // Only use fast solve if A is square and size > 1
  if (A.n_rows == A.n_cols && A.n_rows > 1) {
    bool ok = arma::solve(X, A, B, arma::solve_opts::fast);
    if (!ok) X = arma::pinv(A) * B;
  } else {
    bool ok = arma::solve(X, A, B);
    if (!ok) X = arma::pinv(A) * B;
  }

  return X;
}

// [[Rcpp::export]]
arma::mat optimize_unmix_rcpp_exact(const arma::mat& remaining_raw,
                                    const arma::mat& unmixed,
                                    const arma::mat& spectra,
                                    const arma::vec& pos_thresholds,
                                    const arma::uvec& optimize_idx_r,
                                    const std::vector<arma::mat>& variantsList,
                                    const arma::vec& weights,
                                    const bool weighted = false,
                                    const int nthreads = 1) {

  const arma::uword n_cells = remaining_raw.n_rows;
  const arma::uword n_detectors = remaining_raw.n_cols;
  const arma::uword n_fluors = spectra.n_rows;

  arma::mat result = unmixed;

  arma::vec w = (weights.n_elem == 0)
    ? arma::ones<arma::vec>(n_detectors)
      : weights;
  if (w.n_elem != n_detectors)
    Rcpp::stop("weights length mismatch");

  const bool use_weighted = weighted && (w.n_elem == n_detectors);

#ifdef _OPENMP
  if (nthreads > 0) {
    omp_set_num_threads(nthreads);
  }
  const int nthreads_actual = omp_get_max_threads();
#else
  const int nthreads_actual = 1;
#endif

  // Per-thread preallocation
  std::vector<arma::rowvec> tmp_raw_row(nthreads_actual, arma::rowvec(n_detectors));
  std::vector<arma::rowvec> tmp_cell_unmixed(nthreads_actual, arma::rowvec(n_fluors));
  std::vector<arma::rowvec> tmp_fitted(nthreads_actual, arma::rowvec(n_detectors));
  std::vector<arma::rowvec> tmp_resid(nthreads_actual, arma::rowvec(n_detectors));
  std::vector<arma::rowvec> tmp_resid_candidate(nthreads_actual, arma::rowvec(n_detectors));
  std::vector<arma::mat> tmp_spectra_final(nthreads_actual, spectra);
  std::vector<arma::mat> tmp_pos_unmixing(nthreads_actual);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (arma::uword ci = 0; ci < n_cells; ++ci) {
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    arma::rowvec& raw_row = tmp_raw_row[tid];
    arma::rowvec& cell_unmixed_row = tmp_cell_unmixed[tid];
    arma::rowvec& fitted = tmp_fitted[tid];
    arma::rowvec& resid = tmp_resid[tid];
    arma::rowvec& resid_candidate = tmp_resid_candidate[tid];
    arma::mat& spectra_final = tmp_spectra_final[tid];

    raw_row = remaining_raw.row(ci);
    cell_unmixed_row = result.row(ci);

    // Identify positive fluorophores
    std::vector<arma::uword> pos_vec;
    pos_vec.reserve(n_fluors);
    for (arma::uword f = 0; f < n_fluors; ++f) {
      if (cell_unmixed_row(f) >= pos_thresholds(f))
        pos_vec.push_back(f);
    }
    if (pos_vec.empty()) continue;

    arma::uvec pos_idx(pos_vec.size());
    for (arma::uword i = 0; i < pos_vec.size(); ++i)
      pos_idx(i) = pos_vec[i];

    const arma::uword pos_n = pos_idx.n_elem;
    arma::mat cell_spectra(pos_n, n_detectors);
    for (arma::uword i = 0; i < pos_n; ++i)
      cell_spectra.row(i) = spectra.row(pos_idx(i));

    tmp_pos_unmixing[tid].set_size(pos_n, n_detectors);
    arma::mat unmix_pos = compute_unmixing_matrix_fast(cell_spectra, w, use_weighted);
    tmp_pos_unmixing[tid].rows(0, pos_n - 1) = unmix_pos;

    arma::rowvec unmixed_pos = raw_row * unmix_pos.t();
    fitted = unmixed_pos * cell_spectra;
    resid = raw_row - fitted;
    double error_final = arma::accu(arma::abs(resid));
    spectra_final = spectra;

    // Order fluorophores by intensity
    std::vector<std::pair<double, arma::uword>> order_vec;
    order_vec.reserve(pos_n);
    for (arma::uword i = 0; i < pos_n; ++i)
      order_vec.emplace_back(static_cast<double>(unmixed_pos(i)), i);
    std::sort(order_vec.begin(), order_vec.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    // Variant optimization
    for (const auto& pr : order_vec) {
      arma::uword local_vi = pr.second;
      arma::uword f_idx = pos_idx(local_vi);
      if (f_idx >= variantsList.size()) continue;
      const arma::mat& fl_variants = variantsList[f_idx];
      if (fl_variants.n_rows == 0) continue;

      for (arma::uword v = 0; v < fl_variants.n_rows; ++v) {
        cell_spectra.row(local_vi) = fl_variants.row(v);
        arma::mat unmix_pos_candidate =
          compute_unmixing_matrix_fast(cell_spectra, w, use_weighted);
        arma::rowvec unmixed_pos_candidate = raw_row * unmix_pos_candidate.t();

        fitted = unmixed_pos_candidate * cell_spectra;
        resid_candidate = raw_row - fitted;
        double error_curr = arma::accu(arma::abs(resid_candidate));

        if (error_curr < error_final) {
          error_final = error_curr;
          resid = resid_candidate;
          spectra_final.row(f_idx) = fl_variants.row(v);
          cell_spectra.row(local_vi) = fl_variants.row(v);
          tmp_pos_unmixing[tid].rows(0, pos_n - 1) = unmix_pos_candidate;
          unmixed_pos = unmixed_pos_candidate;
        } else {
          cell_spectra.row(local_vi) = spectra_final.row(f_idx);
        }
      }

      // Recompute after best variant
      cell_spectra.row(local_vi) = spectra_final.row(f_idx);
      arma::mat unmix_pos_after =
        compute_unmixing_matrix_fast(cell_spectra, w, use_weighted);
      tmp_pos_unmixing[tid].rows(0, pos_n - 1) = unmix_pos_after;
      unmixed_pos = raw_row * unmix_pos_after.t();
      fitted = unmixed_pos * cell_spectra;
      resid = raw_row - fitted;
      error_final = arma::accu(arma::abs(resid));
    }

    // Final full unmix (regular solve)
    arma::mat unmix_full = compute_unmixing_matrix_general(spectra_final, w, use_weighted);
    result.row(ci) = raw_row * unmix_full.t();
  }

  return result;
}
