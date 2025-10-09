# optimize_unmix.r

#' @title Optimize Unmix
#' Unmix using the AUtoSpectral method to adapt fluorophore spectral signatures
#' on a per-cell basis.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra. Matrix.
#' @param unmixed Initial unmixing to be optimized. Expectation is that the data
#' will derive from AutoSpectral and will have the autofluorescence component
#' removed. Format: cells (rows) x fluorophores (columns). Matrix.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param pos.thresholds Named numeric. The thresholds to be considered "positive"
#' for a fluorophore in the unnmixing. Expectation: create this by calling
#' `get.spectral.variants`, which generates the thresholds based on the 99.5th
#' percentile of the unstained sample.
#' @param optimize.fluors Vector of fluorophores to be optimized. Must match the
#' rownames in `spectra` and the colnames of `unmixed`.
#' @param variants List of spectral signature matrices for the variants of each
#' fluorophore. Format: Named list with elements fluorophores, each element being
#' a matrix of spectra in format variants (rows) by detectors (columns). Prepare
#' using `get.spectral.variants`.
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `numeric(0)`, in which case OLS unmixing will be used.
#' @param weighted Logical, whether to use ordinary or weighted least squares
#' unmixing as the base algorithm. Default is `FALSE` and will use OLS.
#' @param nthreads Numeric, number of threads to use for parallel processing.
#' Default is `0`, which will allow OpenMP to use all available cores.
#' @param speed Selector for the precision-speed trade-off. Options are the
#' default `fast`, which selects the best spectral fit per cell by updating the
#' predicted values for each fluorophore independently without repeating the
#' unnmixing, `medium` which uses a Woodbury-Sherman-Morrison rank-one updating
#' of the unnmixing matrix for better results and a moderate slow-down, or `slow`,
#' which explicitly recomputes the unmixing matrix for each variant for maximum
#' precision.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.

optimize.unmix <- function( raw.data,
                            unmixed,
                            spectra,
                            pos.thresholds,
                            optimize.fluors,
                            variants,
                            weights = numeric( 0 ),
                            weighted = FALSE,
                            nthreads = 0,
                            speed = "fast" ) {

  # If optimize.fluors are names, try to match to rownames(spectra)
  if ( is.character( optimize.fluors ) ) {

    if ( is.null( rownames( spectra ) ) )
      stop( "spectra must have rownames when optimize.fluors are names" )

    optimize.idx <- match( optimize.fluors, rownames( spectra ) )

    if ( any( is.na( optimize.idx ) ) )
      stop( "Some optimize.fluors not found in spectra rownames" )

  } else {
    optimize.idx <- as.integer( optimize.fluors )
  }

  # convert weights to arma vector: length must equal ncol(raw.data) or length 0
  if ( length( weights ) > 0 && length( weights ) != ncol( raw.data ) ) {
    stop( "weights, if provided, must have length = number of detectors (ncol(raw.data))" )
  }

  variants.list <- lapply( variants, function( x ) x )

  if ( speed == "fast" ) {
    res <- optimize_unmix_rcpp_fast( remaining_raw = as.matrix( raw.data ),
                                     unmixed = as.matrix( unmixed ),
                                     spectra = as.matrix( spectra ),
                                     pos_thresholds = as.numeric( pos.thresholds ),
                                     optimize_idx_r = as.integer( optimize.idx ),
                                     variantsList = variants.list,
                                     weights = as.numeric( weights ),
                                     weighted = weighted,
                                     nthreads = as.integer( nthreads ) )
  } else if ( speed == "medium" ) {
    res <- optimize_unmix_rcpp_woodbury( remaining_raw = as.matrix( raw.data ),
                                         unmixed = as.matrix( unmixed ),
                                         spectra = as.matrix( spectra ),
                                         pos_thresholds = as.numeric( pos.thresholds ),
                                         optimize_idx_r = as.integer( optimize.idx ),
                                         variantsList = variants.list,
                                         weights = as.numeric( weights ),
                                         weighted = weighted,
                                         nthreads = as.integer( nthreads ) )
  } else if ( speed == "slow" ) {
    res <- optimize_unmix_rcpp_exact( remaining_raw = as.matrix( raw.data ),
                                      unmixed = as.matrix( unmixed ),
                                      spectra = as.matrix( spectra ),
                                      pos_thresholds = as.numeric( pos.thresholds ),
                                      optimize_idx_r = as.integer( optimize.idx ),
                                      variantsList = variants.list,
                                      weights = as.numeric( weights ),
                                      weighted = weighted,
                                      nthreads = as.integer( nthreads ) )
  } else {
    warning( "Incorrect method supplied to `speed`. Select from `fast`, `medium` or `slow`." )
    return( unmixed )
  }

  return( as.matrix( res ) )
}
