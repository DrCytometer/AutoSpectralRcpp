# unmix_autospectral_rcpp.r

#' @title Unmix AutoSpectral Rcpp
#'
#' @description
#' Unmix using the AutoSpectral method to extract autofluorescence and optimize
#' fluorophore signatures at the single cell level.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants`. Default is `NULL`.
#' @param weighted Logical, whether to use ordinary or weighted least squares
#' unmixing as the base algorithm. Default is `FALSE` and will use OLS.
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `NULL`, in which case weighting will be done by
#' channel means (Poisson variance). Only used if `weighted`.
#' @param calculate.error Logical, whether to calculate the RMSE unmixing model
#' accuracy and include it as an output. Default is `FALSE`.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#' signature for each cell is determined by which unmixing brings the cell
#' closest to 0 (`use.dist0` = `TRUE`) or by which unmixing minimizes the
#' per-cell residual (`use.dist0` = `FALSE`). Default is `TRUE`.
#' @param verbose Logical, whether to send messages to the console.
#' Default is `TRUE`.
#' @param parallel Logical, whether to use parallel processing for the per-cell
#' unmixing for fluorophore optimization. Default is `FALSE`.
#' @param threads Numeric, number of threads to use for parallel processing.
#' Default is `0`, which will use all available cores. Recommended: use
#' `asp$worker.process.n`.
#' @param speed Selector for the precision-speed trade-off. Options are the
#' default `fast`, which selects the best spectral fit per cell by updating the
#' predicted values for each fluorophore independently without repeating the
#' unnmixing, `medium` which uses a Woodbury-Sherman-Morrison rank-one updating
#' of the unnmixing matrix for better results and a moderate slow-down, or `slow`,
#' which explicitly recomputes the unmixing matrix for each variant for maximum
#' precision.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.autospectral.rcpp <- function( raw.data, spectra, af.spectra,
                                     spectra.variants = NULL,
                                     weighted = FALSE, weights = NULL,
                                     calculate.error = FALSE,
                                     use.dist0 = TRUE,
                                     verbose = TRUE,
                                     parallel = TRUE, threads = 0,
                                     speed = "fast" ) {

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  if ( is.null( af.spectra ) )
    stop( "Multiple AF spectra must be provided." )
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided." )

  # select unmixing algorithm
  if ( weighted )
    unmix <- unmix.wls
  else
    unmix <- unmix.ols

  if ( weighted & is.null( weights ) ) {
    weights <- colMeans( raw.data )
    weights <- 1 / ( weights + 1e-6 )
  }

  fluorophores <- rownames( spectra )
  af.n <- nrow( af.spectra )
  fluorophore.n <- nrow( spectra )
  detector.n <- ncol( spectra )
  combined.spectra <- matrix( NA_real_, nrow = fluorophore.n + 1, ncol = detector.n )
  colnames( combined.spectra ) <- colnames( spectra )
  fluors.af <- c( fluorophores, "AF" )
  rownames( combined.spectra ) <- fluors.af
  combined.spectra[ 1:fluorophore.n, ] <- spectra
  af.only <- is.null( spectra.variants )

  # initial unmixing without any AF
  if ( verbose ) message( "Initializing unmix" )
  unmixed <- unmix( raw.data, spectra, weights )

  # calculate error
  if ( use.dist0 ) {
    error <- rowSums( abs( unmixed[ , fluorophores, drop = FALSE ] ) )
    initial.error <- sum( error )
  } else {
    error <- rowSums( abs( raw.data - ( unmixed %*% spectra ) ) )
    initial.error <- sum( error )
  }

  if ( calculate.error ) {
    # track full residuals for RMSE calculation
    residual <- raw.data - ( unmixed %*% spectra )
    RMSE <- sqrt( mean( residual^2 ) )
  }

  # add AF column, set AF Index
  initial.af <- matrix( 0, nrow = nrow( raw.data ), ncol = 2 )
  colnames( initial.af ) <- c( "AF", "AF Index" )
  unmixed <- cbind( unmixed, initial.af )
  fitted.af <- matrix( 0, nrow = nrow( raw.data ), ncol = ncol( spectra ) )

  # unmix for each af.spectrum
  if ( verbose ) message( "Extracting AF cell-by-cell" )

  for ( af in seq_len( af.n ) ) {
    # set this AF as the spectrum to use
    combined.spectra[ fluorophore.n + 1, ] <- af.spectra[ af, , drop = FALSE ]
    # unmix with this AF
    unmixed.af <- unmix( raw.data, combined.spectra, weights )

    if ( use.dist0 ) {
      # determine which cells have less dist0 error with this AF
      error.af <- rowSums( abs( unmixed.af[ , fluorophores, drop = FALSE ] ) )
      improved <- which( error.af < error )
    } else {
      # determine which cells have less residual error with this AF
      error.af <- rowSums( abs( raw.data - ( unmixed.af %*% combined.spectra ) ) )
      improved <- which( error.af < error )
    }

    error[ improved ] <- error.af[ improved ]
    unmixed[ improved, fluors.af ] <- unmixed.af[ improved, ]
    unmixed[ improved, "AF Index" ] <- af

    # update full residuals if tracking error by RMSE
    if ( calculate.error )
      residual[ improved, ] <- raw.data[ improved, ] -
      ( unmixed.af[ improved, ] %*% combined.spectra )

    # update AF fitted values and residuals with improved cells
    if ( !af.only )
      fitted.af[ improved, ] <- unmixed.af[ improved, "AF", drop = FALSE ] %*%
        af.spectra[ af, , drop = FALSE ]
  }

  if ( calculate.error )
    RMSE <- sqrt( mean( residual^2 ) )

  if ( verbose & calculate.error ) {
    af.error.reduction <- ( initial.error - sum( error ) ) / initial.error * 100
    af.error.reduction <- round( af.error.reduction, 2 )
    message( paste( "Unmixing error reduced by", af.error.reduction, "percent with per-cell AF extraction." ) )
  }

  if ( af.only ) {
    if ( calculate.error )
      return( list(
        RMSE = RMSE,
        unmixed.data = unmixed
      ) )
    else
      return( unmixed )
  }

  ### per cell fluorophore optimization
  # unpack spectral variants
  pos.thresholds <- spectra.variants$thresholds
  variants <- spectra.variants$variants

  if ( is.null( pos.thresholds ) )
    stop( "Check that spectral variants have been calculated using get.spectra.variants" )

  if ( is.null( variants ) )
    stop( "Multiple fluorophore spectral variants must be provided." )
  if ( !( length( variants ) > 1 ) )
    stop( "Multiple fluorophore spectral variants must be provided." )

  # restrict optimization to fluors present in names( variants )
  # if no match, return unmixed with warning
  optimize.fluors <- fluorophores[ fluorophores %in% names( variants ) ]
  if ( !( length( optimize.fluors ) > 0 ) ) {
    warning( "No matching fluorophores between supplied spectra and spectral variants.
             No spectral optimization performed." )
    if ( calculate.error )
      return( list(
        RMSE = RMSE,
        unmixed.data = unmixed
      ) )
    else
      return( unmixed )
  }

  if ( verbose ) message( "Optimizing fluorophore unmixing cell-by-cell" )

  # for fluor optimization, use raw.data with AF fitted component subtracted
  remaining.raw <- raw.data - fitted.af

  unmixed[ , fluorophores ] <- optimize.unmix( raw.data = remaining.raw,
                                               unmixed = unmixed[ , fluorophores ],
                                               spectra = spectra,
                                               pos.thresholds = pos.thresholds,
                                               optimize.fluors = optimize.fluors,
                                               variants = variants,
                                               weights = weights,
                                               weighted = weighted,
                                               nthreads = threads,
                                               speed = speed )

  if ( calculate.error )
    return( list(
      RMSE = RMSE,
      unmixed.data = unmixed
    ) )
  else
    return( unmixed )
}
