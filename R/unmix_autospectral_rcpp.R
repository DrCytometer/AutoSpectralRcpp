# unmix_autospectral_rcpp.r

# Internal helper: validates and unpacks a spectra.variants object.
# Returns a list with elements:
#   $optimize          logical
#   $optimize.fluors   character vector (possibly empty)
#   $variants          list of matrices, named by fluorophore
#   $delta.list        list of delta matrices, named by fluorophore
#   $delta.norms       list of delta norm vectors, named by fluorophore
#   $pos.thresholds    named numeric vector, length == nrow(spectra)
#
# Side-effects: emits messages for skipped fluorophores (consistent with
# existing unmix.autospectral.rcpp behaviour).

.unpack.spectra.variants <- function(
    spectra.variants,
    spectra,
    verbose = TRUE
) {

  fluorophores <- rownames( spectra )

  # default no-optimization return
  empty <- list(
    optimize        = FALSE,
    optimize.fluors = character(0),
    variants        = list(),
    delta.list      = list(),
    delta.norms     = list(),
    pos.thresholds  = numeric(0)
  )

  if ( is.null( spectra.variants ) || length( spectra.variants ) == 0 )
    return( empty )

  if ( is.null( spectra.variants$thresholds ) )
    stop( "Check that spectral variants have been calculated using `get.spectra.variants()`",
          call. = FALSE )
  if ( is.null( spectra.variants$variants ) || !( length( spectra.variants$variants ) > 0 ) )
    stop( "Multiple fluorophore spectral variants must be provided.",
          call. = FALSE )

  # positivity thresholds — full-length vector aligned to spectra row order
  pos.thresholds <- rep( Inf, length( fluorophores ) )
  names( pos.thresholds ) <- fluorophores
  pos.thresholds[ names( spectra.variants$thresholds ) ] <- spectra.variants$thresholds

  variants       <- spectra.variants$variants
  optimize.fluors <- fluorophores[ fluorophores %in% names( variants ) ]

  if ( length( optimize.fluors ) == 0 )
    return( utils::modifyList( empty, list( pos.thresholds = pos.thresholds ) ) )

  # delta.list / delta.norms: use pre-computed if available, else derive
  if ( is.null( spectra.variants$delta.list ) ) {
    message(
      paste(
        "For best results, re-calculate `spectra.variants` using AutoSpectral",
        "version 1.0.0 or greater. See `?get.spectra.variants`."
      )
    )
    delta.list <- lapply( optimize.fluors, function( fl ) {
      variants[[ fl ]] - matrix(
        spectra[ fl, ],
        nrow = nrow( variants[[ fl ]] ),
        ncol = ncol( spectra ),
        byrow = TRUE
      )
    } )
    names( delta.list ) <- optimize.fluors

    delta.norms <- lapply( delta.list, function( d ) sqrt( rowSums( d^2 ) ) )
  } else {
    delta.list  <- spectra.variants$delta.list
    delta.norms <- spectra.variants$delta.norms
  }

  # remove fluorophores with no meaningful variation
  optimize.fluors <- sanitize.optimization.inputs(
    spectra,
    optimize.fluors,
    variants,
    delta.norms
  )

  # apply necessity filter when available
  if ( !is.null( spectra.variants$optimize.recommended ) ) {
    rec     <- spectra.variants$optimize.recommended
    to.skip <- names( rec )[ !rec ]
    to.skip <- to.skip[ to.skip %in% optimize.fluors ]

    if ( length( to.skip ) > 0 ) {
      if ( verbose )
        message( paste0(
          "\033[34m",
          "Skipping variant optimization for (low necessity score): ",
          paste( to.skip, collapse = ", " ),
          "\033[0m"
        ) )
      optimize.fluors <- optimize.fluors[ !optimize.fluors %in% to.skip ]
    }
  }

  list(
    optimize        = length( optimize.fluors ) > 0,
    optimize.fluors = optimize.fluors,
    variants        = variants,
    delta.list      = delta.list,
    delta.norms     = delta.norms,
    pos.thresholds  = pos.thresholds
  )
}

#' @title Unmix AutoSpectral Rcpp
#'
#' @description
#' Unmix using the AutoSpectral method to extract autofluorescence and optimize
#' fluorophore signatures at the single cell level.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#'   detectors in columns. Columns should be fluorescent data only and must
#'   match those of \code{spectra}.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#'   and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#'   between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#'   using \code{get.af.spectra}.
#' @param spectra.variants Named list produced by \code{get.spectra.variants},
#'   or \code{NULL} (default) for AF-only mode.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#'   signature for each cell is determined by which unmixing brings the fluorophore
#'   signals closest to 0 (\code{use.dist0} = \code{TRUE}) or by which unmixing
#'   minimizes the per-cell residual (\code{use.dist0} = \code{FALSE}). Default
#'   is \code{TRUE}. Used by the "legacy" pipeline.
#' @param verbose Logical, default \code{TRUE}. Whether to send messages to the
#'   console.
#' @param speed Default is \code{fast}. Selector for the precision-speed trade-off
#'   in AutoSpectral per-cell fluorophore optimization. Options are \code{slow},
#'   \code{medium} and \code{fast}. From v1.0.0, this controls the number of
#'   variants tested per cell (and per fluorophore). More variants takes longer,
#'   but gives better resolution in some unmixed data. When \code{speed = fast},
#'   as single variant will be tested; for \code{medium}, three will be tested and
#'   for \code{slow}, 10 variants will be tested.
#' @param parallel Logical, default is \code{TRUE}. Always use parallel processing in
#'   the C++ version for faster results.
#' @param threads Integer, number of OpenMP threads. Default \code{1L}.
#' @param n.variants Number of variants to test per cell. Allows explicit control
#'   over the number used, as opposed to \code{speed}, which selects from pre-defined
#'   choices. Providing a numeric value to \code{n.variants} will override \code{speed},
#'   allowing up to \code{n.variants} (or the max available) variants to be tested.
#'   The default is `NULL`, in which case \code{n.variants} will be ignored.
#' @param pipeline Character, one of \code{"legacy"} (default) or \code{"joint"} .
#'   \code{"joint"} uses the new covariance-weighted joint AF + variant optimisation
#'   pipeline (\code{unmix_autospectral_joint_cpp} or its pure-R equivalent).
#'   \code{"legacy"} uses the previous sequential pipeline
#'   (\code{unmix_autospectral_pipeline_cpp}), which selects the AF spectrum first
#'   and then optimises fluorophore variants independently.
#' @param n.af.passes Integer, number of autofluorescence extraction passes per cell.
#'   Default \code{1L}. Only used when \code{pipeline = "joint"}. Higher values allow
#'   a bit more AF to be extracted and may be useful for tissue samples. Neglible
#'   difference is expected for PBMC data.
#' @param n.passes Integer, number of joint fluorophore optimisation passes per cell.
#'   Default \code{1L}. Only used when \code{pipeline = "joint"}. Higher values allow
#'   more swaps to be committed per cell at the cost of additional computation.
#' @param cell.weight Logical. If \code{TRUE}, applies per-cell detector
#'   weighting to the unmixing solve. Weights are \eqn{w_d = 1 /
#'   \max(|\hat{y}_d|, \text{noise.floor})}, where \eqn{\hat{y}_d} is the
#'   fitted signal in detector \eqn{d} from an initial unweighted solve.
#'   Default \code{FALSE}. Useful for ID7000 files. Used in the new
#'   \code{unmix_autospectral_joint_cpp}  pipeline.
#' @param noise.floor Numeric. Lower clamp on the denominator of the per-cell
#'   detector weights, preventing \eqn{w_d \to \infty} for near-dark channels.
#'   Only used when \code{cell.weight = TRUE}. Default \code{1}.
#' @param alpha Numeric, default \code{0.5}. Weighting for balancing residual and
#'   covariance spillover minimization. Used in the new \code{unmix_autospectral_joint_cpp}
#'   pipeline.
#' @param collinear.threshold Numeric, default \code{0.5}. Cosine similarity value to
#'   trigger conflict assessment in assignment of per-cell fluorophore spectral
#'   variants. Higher values will allow committing variants that minimize
#'   spillover from the brighter (higher abundance) fluorophore on a cell, even
#'   if that might incorrectly reduce the calculated abundance of the lower
#'   expression collinear fluorophore being spilled into. The lower the value,
#'   the more conflict resolution will be required and the longer the calculation
#'   will take. The default \code{0.5} is about right.
#' @param joint.pair.resolution Logical, default \code{TRUE}. Whether to perform a
#'   round of "conflict resolution" for collinear fluorophores that compete for
#'   the same spillover/residual space. When \code{TRUE}, a small combinatorial test
#'   will be performed exclusively for the competing pair(s) to attempt to find
#'   an optimal solution.
#' @param refine.af.quantile Numeric, default \code{0.5}. Fraction of cells to
#'   take forward for any additional AF passes (determined by \code{n.af.passes}).
#'   Set to \code{1} to assess all cells in subsequent round; otherwise the highest
#'   abundance AF cells will be reassessed.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.autospectral.rcpp <- function(
    raw.data,
    spectra,
    af.spectra,
    spectra.variants      = NULL,
    use.dist0             = TRUE,
    verbose               = TRUE,
    speed                 = c("fast", "medium", "slow"),
    parallel              = TRUE,
    threads               = 1L,
    n.variants            = NULL,
    pipeline              = c("legacy", "joint"),
    n.af.passes           = 1L,
    n.passes              = 1L,
    cell.weight           = FALSE,
    noise.floor           = 125L,
    alpha                 = 0.5,
    collinear.threshold   = 0.5,
    joint.pair.resolution = TRUE,
    refine.af.quantile    = 0.5
) {

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  if ( is.null( af.spectra ) || nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided.", call. = FALSE )

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

  fluorophores <- rownames( spectra )

  # unpack and validate spectra.variants
  sv <- .unpack.spectra.variants( spectra.variants, spectra, verbose )

  # tell the user what's been selected
  if ( verbose ) {
    if ( sv$optimize )
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
  pass.variants   <- if ( sv$optimize ) sv$variants[ sv$optimize.fluors ]    else list()
  pass.deltas     <- if ( sv$optimize ) sv$delta.list[ sv$optimize.fluors ]  else list()
  pass.thresholds <- if ( sv$optimize ) sv$pos.thresholds                    else numeric(0)
  pass.norms      <- if ( sv$optimize ) sv$delta.norms[ sv$optimize.fluors ] else list()

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
      raw_data_in           = raw.data,
      spectra               = as.matrix( spectra ),
      af_spectra            = as.matrix( af.spectra ),
      fluor_names           = fluorophores,
      pos_thresholds        = as.numeric( pass.thresholds ),
      variants_list         = pass.variants,
      delta_list            = pass.deltas,
      n_passes              = as.integer( n.passes ),
      n_threads             = as.integer( threads ),
      cell_weight           = isTRUE( cell.weight ),
      noise_floor           = as.double( noise.floor ),
      alpha                 = as.numeric( alpha ),
      collinear_thresh      = as.numeric( collinear.threshold ),
      joint_pair_resolution = isTRUE( joint.pair.resolution ),
      n_af_passes           = as.integer( n.af.passes ),
      refine_af_quantile    = as.numeric( refine.af.quantile )
    )
  } else {
    unmixed <- unmix_autospectral_pipeline_cpp(
      raw_data_in     = raw.data,
      spectra         = as.matrix( spectra ),
      af_spectra      = as.matrix( af.spectra ),
      fluor_names     = fluorophores,
      optimize_fluors = sv$optimize.fluors,
      pos_thresholds  = as.numeric( pass.thresholds ),
      variants        = pass.variants,
      delta_list      = pass.deltas,
      delta_norms     = pass.norms,
      k_opt           = as.integer( n.variants ),
      n_threads       = as.integer( threads ),
      optimize        = isTRUE( sv$optimize ),
      use_dist0       = isTRUE( use.dist0 )
    )
  }

  # add names for output
  colnames( unmixed ) <- c( fluorophores, "AF", "AF Index" )

  if ( verbose ) message( "Unmixing complete." )
  return( unmixed )
}
