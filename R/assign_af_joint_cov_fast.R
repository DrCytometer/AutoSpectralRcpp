# assign_af_joint_cov_fast.r

#' @title Assign AF Spectrum By Joint Covariance-Weighted Error, Fast
#'
#' @description
#' Assigns each cell to the best-fitting autofluorescence spectral variant
#' using a joint scoring criterion that multiplies two proportional error
#' terms: a covariance-weighted L1 (abs) fluorophore error and an L2
#' (Euclidean norm) raw-space residual error. Compiled counterpart of
#' `assign.af.joint.cov`, which uses OpenMP for parallel processing. See
#' `assign.af.joint.cov` for the full derivation.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with AF variants in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param threads Numeric, default is `1`.
#'
#' @return The indices of the best-fitting AF spectrum per cell
#'
#' @export

assign.af.joint.cov.fast <- function(raw.data, spectra, af.spectra, threads = 1) {
  # remove AF if present
  if ("AF" %in% rownames(spectra))
    spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE]

  # Call C++ function
  # Ensure inputs are matrices
  af.idx <- assign_af_joint_cov_cpp(
    as.matrix(raw.data),
    as.matrix(spectra),
    as.matrix(af.spectra),
    n_threads = threads
  )

  return(af.idx)
}
