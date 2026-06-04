# unmix_autospectral_rcpp.r

#' @title Unmix AutoSpectral Rcpp
#'
#' @description
#' Unmix using the AutoSpectral method to extract autofluorescence and optimize
#' fluorophore signatures at the single cell level.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#'   detectors in columns. Columns should be fluorescent data only and must
#'   match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#'   and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#'   between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#'   using `get.af.spectra`.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#'   of spectral signature variations for each fluorophore. Prepare using
#'   `get.spectral.variants`. Default is `NULL`.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#'   signature for each cell is determined by which unmixing brings the fluorophore
#'   signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing minimizes the
#'   per-cell residual (`use.dist0` = `FALSE`). Default is `TRUE`.
#' @param verbose Logical, default `TRUE`. Whether to send messages to the console.
#' @param speed Default is `fast`. Selector for the precision-speed trade-off in
#'   AutoSpectral per-cell fluorophore optimization. Options are `slow`, `medium`
#'   and `fast`. From v1.0.0, this controls the number of variants tested per cell
#'   (and per fluorophore). More variants takes longer, but gives better resolution
#'   in some unmixed data. When `speed = fast`, as single variant will be tested;
#'   for `medium`, three will be tested and for `slow`, 10 variants will be tested.
#' @param parallel Logical, default is `TRUE`. Always use parallel processing in
#'   the C++ version for faster results.
#' @param threads Numeric, default is `1`.
#' @param n.variants Number of variants to test per cell. Allows explicit control
#'   over the number used, as opposed to `speed`, which selects from pre-defined
#'   choices. Providing a numeric value to `n.variants` will override `speed`,
#'   allowing up to `n.variants` (or the max available) variants to be tested. The
#'   default is `NULL`, in which case `n.variants` will be ignored.
#' @param pipeline Character, one of `"joint"` (default) or `"legacy"`.
#'   `"joint"` uses the new covariance-weighted joint AF + variant optimisation
#'   pipeline (`unmix_autospectral_joint_cpp` or its pure-R equivalent).
#'   `"legacy"` uses the previous sequential pipeline
#'   (`unmix_autospectral_pipeline_cpp`), which selects the AF spectrum first
#'   and then optimises fluorophore variants independently.
#' @param n.passes Integer, default `2L`. Number of joint optimisation passes
#'   per cell. Only used when `pipeline = "joint"`. Higher values allow more
#'   swaps to be committed per cell at the cost of additional computation.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.autospectral.rcpp <- function(
    raw.data,
    spectra,
    af.spectra,
    spectra.variants = NULL,
    use.dist0 = TRUE,
    verbose = TRUE,
    speed = c("fast", "medium", "slow"),
    parallel = TRUE,
    threads = 1L,
    n.variants = NULL,
    pipeline = c( "joint", "legacy" ),
    n.passes = 2L
) {

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  if ( is.null( af.spectra ) )
    stop( "Multiple AF spectra must be provided.",
          call. = FALSE )
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided.",
          call. = FALSE )

  # check for data/spectra column matching
  raw.data.cols <- colnames( raw.data )
  spectra.cols <- colnames( spectra )

  if ( !identical( raw.data.cols, spectra.cols ) ) {

    # ensure both actually have the same columns before reordering
    if ( all( spectra.cols %in% raw.data.cols ) &&
         length( spectra.cols ) == length( raw.data.cols ) ) {
      # reorder raw.data to match the order of spectra
      raw.data <- raw.data[ , spectra.cols ]
      message( "Columns reordered to match spectra." )
    } else {
      stop( "Column names in spectra and raw.data do not match perfectly;
            cannot reorder by name alone.",
            call. = FALSE)
    }
  }

  # prepare fluorophore names and inputs
  fluorophores <- rownames( spectra )
  detectors <- colnames( spectra )

  # handle optimization logic (sanitize for C++)
  optimize <- FALSE
  optimize.fluors <- character(0)

  # check if optimization is requested and possible
  if ( !is.null( spectra.variants ) && length( spectra.variants ) > 0 ) {

    if ( is.null( spectra.variants$thresholds ) )
      stop( "Check that spectral variants have been calculated using `get.spectra.variants()`",
            call. = FALSE )
    if ( is.null( spectra.variants$variants ) )
      stop( "Multiple fluorophore spectral variants must be provided.",
            call. = FALSE )
    if ( !( length( spectra.variants$variants ) > 0 ) )
      stop( "Multiple fluorophore spectral variants must be provided.",
            call. = FALSE )

    # set positivity thresholds vector
    pos.thresholds <- rep( Inf, length( fluorophores ) )
    names( pos.thresholds ) <- fluorophores
    # fill with data
    pos.thresholds[ names( spectra.variants$thresholds ) ] <- spectra.variants$thresholds

    # unpack spectral variants
    variants <- spectra.variants$variants

    # restrict optimization to fluors present in names( variants )
    optimize.fluors <- fluorophores[ fluorophores %in% names( variants ) ]

    if ( length( optimize.fluors ) > 0 ) {

      if ( is.null( spectra.variants$delta.list ) ) {
        message(
          paste(
            "For best results, re-calculate `spectra.variants` using AutoSpectral",
            "version 1.0.0 or greater. See `?get.spectra.variants`."
          )
        )
        # calculate deltas for each fluorophore's variants
        delta.list <- lapply( optimize.fluors, function( fl ) {
          variants[[ fl ]] - matrix(
            spectra[ fl, ],
            nrow = nrow( variants[[ fl ]] ),
            ncol = ncol( spectra ),
            byrow = TRUE
          )
        } )
        names( delta.list ) <- optimize.fluors

        # precompute delta norms
        delta.norms <- lapply( delta.list, function( d ) {
          sqrt( rowSums( d^2 ) )
        } )
      } else {
        delta.list <- spectra.variants$delta.list
        delta.norms <- spectra.variants$delta.norms
      }

      # remove fluors where delta.norms are all zero/empty
      optimize.fluors <- sanitize.optimization.inputs(
        spectra,
        optimize.fluors,
        variants,
        delta.norms
      )

      # filter using pre-computed necessity scores when available
      if ( !is.null( spectra.variants$optimize.recommended ) ) {
        rec     <- spectra.variants$optimize.recommended
        to.skip <- names( rec )[ !rec ]
        to.skip <- to.skip[ to.skip %in% optimize.fluors ]

        if ( length( to.skip ) > 0 ) {
          if ( verbose ) {
            message( paste0(
              "\033[34m",
              "Skipping variant optimization for (low necessity score): ",
              paste( to.skip, collapse = ", " ),
              "\033[0m"
            ) )
          }
          optimize.fluors <- optimize.fluors[ !optimize.fluors %in% to.skip ]
        }
      }

      # anything left?
      if ( length( optimize.fluors ) > 0 ) {
        optimize <- TRUE
      }
    }
  }

  # tell the user what's been selected
  if ( verbose ) {
    if ( optimize )
      message( "Per-cell autofluorescence extraction and fluorophore optimization will be performed." )
    else
      message( "Only per-cell autofluorescence extraction will be performed." )
  }

  # set variant count based on speed
  if ( is.null( n.variants ) ) {
    n.variants <- switch( match.arg( speed ), "fast" = 1L, "medium" = 3L, "slow" = 10L )
  }

  pipeline <- match.arg( pipeline )

  # handle NULLs safely
  pass.variants   <- if ( optimize ) variants      else list()
  pass.deltas     <- if ( optimize ) delta.list    else list()
  pass.norms      <- if ( optimize ) delta.norms   else list()
  pass.thresholds <- if ( optimize ) pos.thresholds else numeric(0)

  if ( verbose ) {
    if ( pipeline == "joint" )
      message( "Running joint AutoSpectral C++ pipeline..." )
    else
      message( "Running legacy AutoSpectral C++ pipeline..." )
  }

  # call C++ pipeline
  # returns matrix [Fluors | AF_val | AF_idx]
  if ( pipeline == "joint" ) {
    unmixed <- unmix_autospectral_joint_cpp(
      raw_data_in   = raw.data,
      spectra       = as.matrix( spectra ),
      af_spectra    = as.matrix( af.spectra ),
      fluor_names   = fluorophores,
      pos_thresholds = as.numeric( pass.thresholds ),
      variants_list  = pass.variants,
      delta_list     = pass.deltas,
      n_passes       = n.passes,
      n_threads      = threads
    )
  } else {
    unmixed <- unmix_autospectral_pipeline_cpp(
      raw_data_in = raw.data,
      spectra = as.matrix( spectra ),
      af_spectra = as.matrix( af.spectra ),
      fluor_names = fluorophores,
      optimize_fluors = optimize.fluors,
      pos_thresholds = as.numeric( pass.thresholds ),
      variants = pass.variants,
      delta_list = pass.deltas,
      delta_norms = pass.norms,
      k_opt = n.variants,
      n_threads = threads,
      optimize = optimize,
      use_dist0 = use.dist0
    )
  }

  # add names for output
  colnames( unmixed ) <- c( fluorophores, "AF", "AF Index" )

  if ( verbose ) message( "Unmixing complete." )
  return( unmixed )
}
