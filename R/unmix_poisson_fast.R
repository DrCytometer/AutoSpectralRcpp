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
#' @return Matrix of unmixed fluorophore intensities
#' @export


unmix.poisson.fast <- function( raw.data, spectra, maxit = 100,  tol = 1e-6,
                                n_threads = 0 ) {

  # set any negative or zero values to a very small number
  raw.data <- pmax( raw.data, 1e-6 )
  # spectra <- abs( spectra )
  spectra.t <- t( spectra )
  spectra.t[ spectra.t <= 0 ] <- 1e-6

  # WLS initiation
  beta.init <- unmix.wls( raw.data, spectra )
  beta.init[ beta.init <= 0 ] <- 1e-6

  # call fast C++ IRLS
  unmixed <- poisson_irls_rcpp_parallel(
    raw_data = raw.data,
    spectra = spectra,
    beta_init = beta.init,
    maxit = maxit,
    tol = tol,
    n_threads = n_threads
  )

  colnames( unmixed ) <- rownames( spectra )

  return( unmixed )
}
