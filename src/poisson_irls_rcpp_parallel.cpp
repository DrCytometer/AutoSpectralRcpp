// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat poisson_irls_rcpp_parallel(const arma::mat& raw_data_in,
                                     const arma::mat& spectra,
                                     const arma::mat& beta_init_in,
                                     const int maxit = 25,
                                     const double tol = 1e-6,
                                     const int n_threads = 1) {

  int n_cells = raw_data_in.n_rows;
  arma::mat X = spectra.t(); // detectors x fluorophores
  arma::mat result(n_cells, spectra.n_rows); // preallocate result

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel for
  for (int i = 0; i < n_cells; i++) {
    arma::rowvec y = clamp(raw_data_in.row(i), 1e-6, datum::inf);
    arma::colvec beta = clamp(beta_init_in.row(i).t(), 1e-6, datum::inf);

    for (int iter = 0; iter < maxit; iter++) {
      arma::colvec eta = X * beta;
      eta = clamp(eta, 1e-6, datum::inf);

      arma::colvec mu = eta;
      arma::colvec z = eta + (y.t() - mu);
      arma::colvec w = 1.0 / mu;

      arma::mat Xw = X.each_col() % sqrt(w);
      arma::colvec zw = sqrt(w) % z;
      arma::colvec beta_new;

      bool success = true;
      try {
        beta_new = solve(Xw.t() * Xw, Xw.t() * zw);
      } catch (...) {
        success = false;
      }

      if (!success) break;

      if (norm(beta_new - beta, 2) < tol) {
        beta = beta_new;
        break;
      }

      beta = beta_new;
    }

    result.row(i) = beta.t();
  }

  return result;
}
