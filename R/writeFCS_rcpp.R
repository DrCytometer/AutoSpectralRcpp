# writeFCS.R  (AutoSpectralRcpp)
#
# Drop-in replacement for AutoSpectral::writeFCS using fcs_rcpp_write_data().
#
# The R side handles everything the pure-R version does: keyword enforcement,
# TEXT segment construction, and the offset self-resolution loop.  The C++
# function takes over for the binary write, eliminating the transient
# allocations that R's writeBin path requires.
#
# Memory comparison for 500k events x 183 channels (~732 MB matrix):
#   R writeBin path:   t(mat) ~732 MB + as.vector ~732 MB + float buf ~366 MB
#                      = ~1.8 GB transient on top of the input matrix
#   C++ path:          one CHUNK_ROWS * n_col * 4-byte float buffer (~4-96 MB)
#                      = essentially zero transient relative to the matrix

#' @title Write FCS File (Rcpp-accelerated)
#'
#' @description
#' A drop-in replacement for \code{AutoSpectral::writeFCS} that uses a compiled
#' C++ writer for the binary data segment.  Eliminates the large transient
#' allocations made by R's \code{writeBin(as.vector(t(mat)), ...)} path,
#' substantially reducing peak memory for large high-channel-count files.
#' Automatically replaces \code{AutoSpectral::writeFCS} in the
#' \code{AutoSpectral} namespace when \code{AutoSpectralRcpp} is loaded
#' (see \code{AutoSpectralRcpp::::.onLoad}).
#'
#' @param mat The expression data matrix (events x channels).
#' @param keys Named list of FCS TEXT segment keywords.
#' @param file.name The output filename (not including directory).
#' @param output.dir Directory in which to write the file.
#'
#' @export
#'
#' @seealso \code{\link[AutoSpectral:writeFCS]{AutoSpectral::writeFCS}}

writeFCS <- function(mat, keys, file.name, output.dir) {

  delim <- "|"

  # --- 1. Enforce mandatory keywords ----------------------------------------
  keys[["$TOT"]]      <- as.character(nrow(mat))
  keys[["$PAR"]]      <- as.character(ncol(mat))
  keys[["$DATATYPE"]] <- "F"
  keys[["$BYTEORD"]]  <- "1,2,3,4"
  keys[["$NEXTDATA"]] <- "0"

  keys[["$BEGINDATA"]] <- "0"   # placeholder, resolved below
  keys[["$ENDDATA"]]   <- "0"

  # --- 2. Build initial TEXT segment ----------------------------------------
  text.segment <- paste0(
    delim,
    paste0(names(keys), delim, unlist(keys), delim, collapse = "")
  )

  TEXT.start       <- 58L
  TEXT.end         <- nchar(text.segment, "bytes") + TEXT.start - 1L
  data.stream.bytes <- nrow(mat) * ncol(mat) * 4L   # float32

  # --- 3. Resolve self-referential $BEGINDATA / $ENDDATA --------------------
  # TEXT.end grows as the offset strings get longer; iterate to convergence.
  kw.len.old <- 2L
  repeat {
    DATA.start  <- TEXT.end + 1L
    DATA.end    <- DATA.start + data.stream.bytes - 1L
    kw.len.new  <- nchar(DATA.start) + nchar(DATA.end)
    if (kw.len.new > kw.len.old) {
      TEXT.end   <- TEXT.end + kw.len.new - kw.len.old
      kw.len.old <- kw.len.new
    } else break
  }

  # Stamp resolved offsets into TEXT segment
  text.segment <- sub(
    "\\|\\$BEGINDATA\\|0\\|",
    paste0("|$BEGINDATA|", DATA.start, "|"),
    text.segment
  )
  text.segment <- sub(
    "\\|\\$ENDDATA\\|0\\|",
    paste0("|$ENDDATA|", DATA.end, "|"),
    text.segment
  )

  # --- 4. Build 58-byte header ----------------------------------------------
  # Per FCS 3.1 spec: if any offset exceeds 8 characters, write 0 in the
  # header field and rely on the TEXT segment keywords instead.
  h_t <- if (nchar(TEXT.start) > 8 || nchar(TEXT.end) > 8) c(0L, 0L) else c(TEXT.start, TEXT.end)
  h_d <- if (nchar(DATA.start) > 8 || nchar(DATA.end) > 8) c(0L, 0L) else c(DATA.start, DATA.end)

  header <- paste0(
    sprintf("%-10s", "FCS3.1"),
    sprintf("%8d%8d", h_t[1], h_t[2]),
    sprintf("%8d%8d", h_d[1], h_d[2]),
    sprintf("%8d%8d", 0L, 0L)   # analysis segment offsets
  )

  # Validate: header must be exactly 58 bytes (ASCII only, so chars == bytes)
  if (nchar(header, "bytes") != 58L)
    stop("Internal error: FCS header is ", nchar(header, "bytes"),
         " bytes, expected 58. Please report this as a bug.")

  # --- 5. Hand off binary write to C++ -------------------------------------
  # fcs_rcpp_write_data() writes header + text_segment + data (transposed,
  # cast to float32, chunked) + footer.  swap=FALSE since $BYTEORD is
  # always "1,2,3,4" (little-endian) for files written by this package.
  fcs_rcpp_write_data(
    file_path    = file.path(output.dir, file.name),
    header       = header,
    text_segment = text.segment,
    data_mat     = mat,
    swap         = FALSE
  )

  invisible(NULL)
}
