# unmix_poisson_fast.r

#' @title Poisson IRLS Unmixing (via C++)
#'
#' @description
#' Unmixes spectral flow cytometry data via iteratively reweighted least squares
#'     (IRLS) towards a Poisson distribution, starting with weighted least squares
#'     (WLS). This is a wrapper function in R around the C++ worker function.
#'
#' @param raw.data Matrix of cells x detectors
#' @param spectra Matrix of fluorophores x detectors
#' @param maxit Numeric. The maximum number of iterations to be performed.
#'     Default is 100.
#' @param tol Numeric. Tolerance for convergence. Default is 1e-6. Higher numbers
#'     Will converge faster; smaller numbers may improve convergence.
#' @param n_threads Numeric. Number of parallel processes to be used in C++.
#'     Recommended is asp$worker.process.n (or one less than availableCores()).
#' @param divergence.threshold Numeric. Threshold for reversion to the initial
#'     WLS unmixing for any point that changes dramatically during Poisson IRLS.
#'     Default is 1e4.
#' @return Matrix of unmixed fluorophore intensities
#' @export


unmix.poisson.fast <- function( raw.data, spectra, maxit = 100,  tol = 1e-6,
                                n_threads = 0, divergence.threshold = 1e4 ) {

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
  wls.unmix <- unmix.wls( raw.data, spectra )
  beta.init <- wls.unmix
  beta.init[ beta.init <= 0 ] <- 1e-6

  # call fast C++ IRLS
  unmixed <- poisson_irls_rcpp_parallel(
    raw_data = raw.data,
    spectra = spectra,
    beta_init = beta.init,
    maxit = maxit,
    tol = tol,
    n_threads = n_threads,
    divergence_threshold = divergence.threshold
  )

  # check for nonconvergent points from C++ unmixing and revert towards WLS
  non.convergent <- which( unmixed == 1e-6, arr.ind = TRUE )
  unmixed[ non.convergent, ] <- 0.1*wls.unmix[ non.convergent, ] + 0.9*unmixed[ non.convergent, ]

  colnames( unmixed ) <- rownames( spectra )

  return( unmixed )
}
