// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// Helper function to compute Poisson deviance
inline double poisson_deviance(const arma::rowvec& y, const arma::colvec& mu) {
  double dev = 0.0;
  for (size_t i = 0; i < y.n_elem; i++) {
    if (y[i] > 0) {
      dev += 2.0 * (y[i] * std::log(y[i] / mu[i]) - (y[i] - mu[i]));
    } else {
      dev += 2.0 * mu[i];
    }
  }
  return dev;
}

//' Poisson IRLS Parallel
//'
//' @export
//' @param raw_data_in matrix of cells x detectors
//' @param spectra matrix of fluorophores x detectors
//' @param beta_init_in initial unmixed matrix of cells x fluorophores
//' @param maxit maximum number of IRLS iterations
//' @param tol convergence tolerance
//' @param n_threads number of threads to use (OpenMP)
//' @param divergence_threshold maximum ratio of a cell's deviance to its
//' initial (WLS-seeded) deviance before that iteration is treated as
//' divergent and the cell reverts towards `beta_init`
//' @param max_halving_steps maximum number of halving steps to use to approach
//' convergence
//' @param noise_floor Numeric scalar or per-detector vector, default `125`.
//' Lower clamp on the denominator of the IRLS weight, preventing
//' \eqn{w_d \to \infty} for near-dark channels. Signal units, same
//' convention as `noise.floor` elsewhere in the package.
//' @param final_validation_ratio Numeric, default `1.5`. A cell's final
//' deviance is allowed to exceed its initial (WLS-seeded) deviance by at
//' most this multiple before the cell is reverted to `beta_init`.
//' @return A list with `beta` (matrix of cells x fluorophores) and
//' `converged` (numeric vector, one entry per cell: `0` converged within
//' `tol`, `1` reverted mid-loop due to divergence, `2` reverted because the
//' final fit produced an invalid (near-zero or negative) predicted mean,
//' `3` reverted because final deviance exceeded `final_validation_ratio`
//' times initial deviance, `4` reached `maxit` without meeting `tol` but the
//' last accepted step was kept)
// [[Rcpp::export]]
Rcpp::List poisson_irls_rcpp_parallel(const arma::mat& raw_data_in,
                                      const arma::mat& spectra,
                                      const arma::mat& beta_init_in,
                                      const int maxit = 25,
                                      const double tol = 1e-6,
                                      const int n_threads = 1,
                                      const double divergence_threshold = 1e4,
                                      const int max_halving_steps = 20,
                                      Rcpp::Nullable<Rcpp::NumericVector> noise_floor = R_NilValue,
                                      const double final_validation_ratio = 1.5) {

  int n_cells = raw_data_in.n_rows;
  arma::mat X = spectra.t(); // detectors x fluorophores
  int n_detectors = X.n_rows;
  arma::mat result(n_cells, spectra.n_rows);
  arma::vec status(n_cells, arma::fill::zeros); // 0 = converged, 1 = reverted, 2 = stalled at maxit

  // Per-detector noise floor for the IRLS weight only (signal units, same
  // convention as noise_floor in unmix_autospectral_joint_pipeline.cpp).
  // Never applied to the deviance, which must reflect the true fitted mean.
  arma::vec noise_floor_vec;
  if (noise_floor.isNotNull()) {
    noise_floor_vec = Rcpp::as<arma::vec>(noise_floor.get());
  }
  if (noise_floor_vec.is_empty()) {
    noise_floor_vec.set_size(n_detectors);
    noise_floor_vec.fill(125.0);
  } else if (noise_floor_vec.n_elem == 1) {
    double scalar_floor = noise_floor_vec[0];
    noise_floor_vec.set_size(n_detectors);
    noise_floor_vec.fill(scalar_floor);
  } else if ((int)noise_floor_vec.n_elem != n_detectors) {
    Rcpp::stop("'noise_floor' must have length 1 or match the number of detectors.");
  }

  double abs_tol = tol;

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel
{
  // Thread-local buffers
  arma::colvec eta, z, w, sqrtw, beta_new, beta_trial, eta_trial, beta_old;
  arma::mat Xw;

#pragma omp for
  for (int i = 0; i < n_cells; i++) {

    arma::rowvec y = raw_data_in.row(i);
    arma::colvec beta = beta_init_in.row(i).t();
    // spectra are non-negative (normalized 0-1), so flooring beta here
    // guarantees eta = X * beta >= 0 at the starting point. Without this,
    // step-halving (which only interpolates toward beta_old) can never
    // escape a starting point whose eta is already invalid.
    beta.transform( [](double val) { return (val <= 0.0) ? 1e-6 : val; } );
    arma::colvec beta_init = beta;

    bool converged = false;
    bool diverging = false;

    // Initial deviance
    eta = X * beta;
    eta.clamp(1e-5, datum::inf);  // Stronger lower bound
    double dev_old = poisson_deviance(y, eta);
    double dev_init = dev_old;

    for (int iter = 0; iter < maxit && !converged && !diverging; iter++) {

      beta_old = beta;

      // Current linear predictor and mean
      eta = X * beta;
      eta.clamp(1e-5, datum::inf);  // Keep strictly positive for the solve below

      // Poisson variance = mean, floored per detector at noise_floor so
      // near-dark channels get a bounded weight instead of w -> 1/1e-5
      w = 1.0 / arma::max(eta, noise_floor_vec);
      z = y.t();
      sqrtw = sqrt(w);

      // Build weighted least squares system
      Xw = X.each_col() % sqrtw;

      // Solve using QR decomposition (more stable than normal equations)
      try {
        arma::mat Q, R;
        arma::qr_econ(Q, R, Xw);
        arma::colvec Qty = Q.t() * (sqrtw % z);
        beta_new = solve(trimatu(R), Qty);
      } catch (...) {
        break;
      }

      // Step halving with strict deviance monitoring
      bool step_ok = false;
      double step_size = 1.0;
      double dev_new = datum::inf;

      // First check if full step is valid
      beta_trial = beta_old + (beta_new - beta_old);
      beta_trial.clamp(-1e10, 1e10);
      eta_trial = X * beta_trial;

      // If full step produces negative eta, find the boundary
      if (eta_trial.min() <= 1e-10) {
        arma::colvec eta_old = X * beta_old;
        arma::colvec eta_proposed = eta_trial;
        arma::colvec eta_delta = eta_proposed - eta_old;

        step_size = 1.0;
        for (size_t j = 0; j < eta_old.n_elem; j++) {
          if (eta_delta[j] < -1e-10 && eta_old[j] > 1e-10) {
            double step_limit = (1e-8 - eta_old[j]) / eta_delta[j];
            if (step_limit < step_size && step_limit > 0) {
              step_size = step_limit;
            }
          }
        }

        step_size *= 0.99;  // Back off slightly from boundary
      }

      // Now do step halving from this starting point
      for (int halve = 0; halve <= max_halving_steps; halve++) {

        beta_trial = beta_old + step_size * (beta_new - beta_old);
        beta_trial.clamp(-1e10, 1e10);

        eta_trial = X * beta_trial;

        // Check if all eta values are positive
        if (eta_trial.min() > 1e-5) {  // Stronger constraint
          dev_new = poisson_deviance(y, eta_trial);

          // Accept if deviance decreased
          if (dev_new < dev_old) {
            beta = beta_trial;
            step_ok = true;
            break;
          }

          // On last attempt, accept if not too much worse
          if (halve == max_halving_steps && dev_new < dev_old * 1.001) {
            beta = beta_trial;
            step_ok = true;
            break;
          }
        }

        step_size *= 0.5;
      }

      if (!step_ok) {
        // Couldn't find improving step
        break;
      }

      // Check for divergence (deviance increasing significantly)
      if (dev_new > dev_old * 1.1 || dev_new > dev_init * divergence_threshold) {
        diverging = true;
        beta = beta_init;  // Revert to initial
        break;
      }

      // Check convergence based on deviance change BEFORE updating dev_old
      double dev_change = std::abs(dev_new - dev_old) / (std::abs(dev_old) + 0.1);

      // Update deviance for next iteration
      dev_old = dev_new;

      if (dev_change < abs_tol) {
        converged = true;
      }

      // Also check parameter convergence
      double param_change = norm(beta - beta_old, 2);
      double param_norm = norm(beta_old, 2);
      if (param_change < abs_tol * (param_norm + abs_tol)) {
        converged = true;
      }
    }

    // Final validation - revert if final solution is worse than initial.
    // status: 0 converged, 1 reverted mid-loop (divergence), 2 reverted
    // (final eta invalid), 3 reverted (final deviance too far above
    // dev_init), 4 stalled at maxit without reverting
    double cell_status;
    eta = X * beta;
    if (diverging) {
      cell_status = 1.0;
    } else if (eta.min() < 1e-5) {
      beta = beta_init;
      cell_status = 2.0;
    } else {
      eta.clamp(1e-5, datum::inf);
      double dev_final = poisson_deviance(y, eta);

      if (dev_final > dev_init * final_validation_ratio) {
        beta = beta_init;
        cell_status = 3.0;
      } else {
        cell_status = converged ? 0.0 : 4.0;
      }
    }

    result.row(i) = beta.t();
    status(i) = cell_status;
  }
}

return Rcpp::List::create(
  Rcpp::Named("beta") = result,
  Rcpp::Named("converged") = status
);
}
