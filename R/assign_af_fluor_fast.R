# assign_af_fluor_fast.r

#' @title Assign AF Spectrum By Fluorophore Projection Fast
#'
#' @description
#' Projects the autofluorescence spectral variants into fluorophore unmixed
#' space to determine which best fits each cell (event). A fast approximation for
#' brute force sequential unmixing method in early versions of AutoSpectral.
#' Provides essentially identical results to minimization of fluorophore signal
#' (dist0 method). Substantially faster. Use L1 (absolute value) minimization,
#' which works better than the standard L2 (squared error). Wrapper for the C++
#' version, which uses OpenMP for parallel processing.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param threads Numeric, default is `1`.
#' 
#' @return The indices of the best-fitting AF spectrum per cell
#' 
#' @export

assign.af.fluor.fast <- function(raw.data, spectra, af.spectra, threads = 1) {
  # remove AF if present
  if ("AF" %in% rownames(spectra)) 
    spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE]
  
  # Call C++ function
  # Ensure inputs are matrices
  af.idx <- assign_af_fluor(
    as.matrix(raw.data), 
    as.matrix(spectra), 
    as.matrix(af.spectra), 
    n_threads = threads
  )
  
  return(af.idx)
}