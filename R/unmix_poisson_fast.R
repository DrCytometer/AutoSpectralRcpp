# unmix_poisson_fast.r

#' @title Poisson IRLS Unmixing (via C++)
#'
#' @description
#' Unmixes spectral flow cytometry data via iteratively reweighted least squares
#' (IRLS) towards a Poisson distribution, starting with weighted least squares
#' (WLS). This is a wrapper function in R around the C++ worker function.
#'
#' @param raw.data Matrix of cells x detectors
#' @param spectra Matrix of fluorophores x detectors
#' @param weights Optional numeric vector of weights, one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means.
#' @param maxit Numeric. The maximum number of iterations to be performed.
#' Default is `100`
#' @param tol Numeric. Tolerance for convergence. Default is 1e-6. Higher
#' numbers Will converge faster; smaller numbers may improve convergence.
#' @param n_threads Numeric. Number of parallel processes to be used in C++.
#' Recommended is `asp$worker.process.n` (or one less than `availableCores()` ).
#' @param divergence.threshold Numeric. Maximum allowed ratio of a cell's
#' current deviance to its initial (WLS-seeded) deviance before that cell is
#' treated as divergent mid-iteration and reverted towards its WLS estimate.
#' Default is `10`.
#' @param divergence.handling String. How to handle cells flagged as
#' non-convergent (reverted or stalled) by the C++ IRLS. Options are `NonNeg`
#' (negative components of the returned estimate are clamped to zero), `WLS`
#' (the cell is replaced entirely by its WLS initial estimate) or `Balance`
#' (a weighted average of the WLS estimate and the returned estimate).
#' Default is `Balance`.
#' @param balance.weight Numeric. Weight given to the WLS estimate when
#' averaging under `Balance`; `1 - balance.weight` is given to the returned
#' IRLS estimate. Used only when `divergence.handling = "Balance"`. Default
#' is `0.5`.
#' @param noise.floor Numeric scalar or per-detector vector, default `125`.
#' Lower clamp on the denominator of the IRLS weight, preventing the weight
#' on near-dark channels from growing without bound. Signal units, same
#' convention as `noise.floor` elsewhere in the package.
#'
#' @return Matrix of unmixed fluorophore intensities
#'
#' @export

unmix.poisson.fast <- function(
    raw.data,
    spectra,
    weights = NULL,
    maxit = 100,
    tol = 1e-6,
    n_threads = 0,
    divergence.threshold = 10,
    divergence.handling = "Balance",
    balance.weight = 0.5,
    noise.floor = 125
) {

  # safety checks on inputs so we don't crash when calling C++
  if ( !is.matrix( raw.data ) || !is.matrix( spectra ) ) {
    stop( "Check structure of 'raw.data' and 'spectra'. These should be matrices." )
  }
  if ( !is.numeric( maxit ) || length( maxit ) != 1 || maxit <= 0 ) {
    stop( "'maxit' must be a positive numeric value of length 1" )
  }
  if ( !is.numeric( tol ) || length( tol ) != 1 || tol <= 0 ) {
    stop( "'tol' must be a positive numeric value of length 1" )
  }
  if ( !is.numeric( n_threads ) || length( n_threads ) != 1 || n_threads < 0 ) {
    stop( "'n_threads' must be a non-negative numeric value of length 1" )
  }
  if ( !is.numeric( divergence.threshold) || length( divergence.threshold ) != 1 ||
      divergence.threshold <= 0 ) {
    stop( "'divergence.threshold' must be a positive numeric value of length 1" )
  }
  if ( ncol( raw.data ) != ncol( spectra ) ) {
    stop( "'raw.data' and 'spectra' must both contain detectors as the columns." )
  }


  # set any negative or zero values to a very small number
  raw.data <- pmax( raw.data, 1e-6 )
  spectra[ spectra <= 0 ] <- 1e-6

  # WLS initial unmixing
  wls.unmix <- unmix.wls.fast( raw.data, spectra, weights, noise.floor = noise.floor )
  beta.init <- wls.unmix

  # call fast C++ IRLS
  irls.result <- poisson_irls_rcpp_parallel(
    raw_data = raw.data,
    spectra = spectra,
    beta_init = beta.init,
    maxit = maxit,
    tol = tol,
    n_threads = n_threads,
    divergence_threshold = divergence.threshold,
    noise_floor = noise.floor
  )

  unmixed <- irls.result$beta

  # cells the C++ side flagged as reverted (divergence or failed final
  # validation) or stalled (hit maxit without meeting tol)
  non.convergent <- which( irls.result$converged != 0 )

  if ( divergence.handling == "Balance" ) {
    unmixed[ non.convergent, ] <- balance.weight*wls.unmix[ non.convergent, ] +
      ( 1 - balance.weight )*unmixed[ non.convergent, ]
  } else if ( divergence.handling == "WLS" ) {
    unmixed[ non.convergent, ] <- wls.unmix[ non.convergent, ]
  } else if ( divergence.handling == "NonNeg" ) {
    reverted.sub <- unmixed[ non.convergent, , drop = FALSE ]
    reverted.sub[ reverted.sub < 0 ] <- 0
    unmixed[ non.convergent, ] <- reverted.sub
  }

  colnames( unmixed ) <- rownames( spectra )

  return( unmixed )
}
