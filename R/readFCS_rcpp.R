# readFCS.R  (AutoSpectralRcpp)
#
# Drop-in replacement for AutoSpectral::readFCS using fcs_rcpp_read_data().
#
# Preserves the full public interface of AutoSpectral::readFCS:
#   - fcs.path, return.keywords, start.row, end.row
#   - same return value: matrix, or list(data, keywords) when return.keywords=TRUE
#
# Advantages over the pure-R version:
#   - reads float32 directly into the output matrix with no intermediate
#     double vector, so peak memory is ~1x the matrix size rather than ~2x
#   - the C++ read loop is substantially faster for large channel counts
#     (e.g. ID7000 with 183 channels)
#
# Keyword parsing and row-subsetting remain in R

#' @title Read FCS File (Rcpp-accelerated)
#'
#' @description
#' A drop-in replacement for \code{AutoSpectral::readFCS} that uses a compiled
#' C++ reader for the binary data segment.  Lower peak memory usage and faster
#' than the pure-R version, particularly for high-channel-count instruments
#' (e.g. ID7000, Xenith).  Automatically replaces \code{AutoSpectral::readFCS}
#' in the \code{AutoSpectral} namespace when \code{AutoSpectralRcpp} is loaded
#' (see \code{AutoSpectralRcpp::::.onLoad}).
#'
#' @param fcs.path A character string specifying the file path for the .fcs
#'   file to be read.
#' @param return.keywords Logical, default \code{FALSE}. When \code{TRUE},
#'   returns a list with elements \code{data} (the expression matrix) and
#'   \code{keywords} (named list of TEXT segment keywords).
#' @param start.row Optional numeric. First event (row) to read. Default
#'   \code{NULL} reads from the first event.
#' @param end.row Optional numeric. Last event (row) to read. Default
#'   \code{NULL} reads to the last event.
#' @param columns Optional character vector of \code{$PnN} channel names to
#'   read. Default \code{NULL} reads and returns every channel (unchanged
#'   behaviour). Supplying a subset avoids decoding and materializing columns
#'   that will be discarded immediately afterward -- useful for wide panels
#'   (e.g. ID7000, Xenith) where only spectral + a few passthrough channels
#'   are actually needed.
#'
#' @return If \code{return.keywords = FALSE} (default): a numeric matrix with
#'   events in rows and channels in columns.  If \code{return.keywords = TRUE}:
#'   a list with elements \code{$data} and \code{$keywords}.
#'
#' @export
#'
#' @seealso \code{AutoSpectral::readFCS}
#'
#' @references
#' Granjeud, Samuel. \emph{flowCoreUtils}.
#' \url{https://github.com/i-cyto/flowCoreUtils}
#'
#' Laniewski, Nathan. \emph{flowstate}.
#' \url{https://github.com/nlaniewski/flowstate}
#'
#' Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J, Jiang M,
#' Finak G (2025). \emph{flowCore: Basic structures for flow cytometry data}.
#' \doi{10.18129/B9.bioc.flowCore}

readFCS <- function(
    fcs.path,
    return.keywords = FALSE,
    start.row       = NULL,
    end.row         = NULL,
    columns         = NULL
) {

  fcs.path <- path.expand(fcs.path)

  # --- 1. Parse TEXT segment keywords ---------------------------------------
  con <- file(fcs.path, open = "rb")
  on.exit(close(con))

  header <- readChar(con, 58)
  txt.st <- as.numeric(trimws(substr(header, 11, 18)))
  txt.en <- as.numeric(trimws(substr(header, 19, 26)))

  seek(con, txt.st)
  txt.raw  <- readBin(con, "raw", n = txt.en - txt.st + 1)
  txt      <- rawToChar(txt.raw[txt.raw != as.raw(0)])

  delim <- substr(txt, 1, 1)
  kv    <- strsplit(txt, delim, fixed = TRUE)[[1]]
  if (length(kv) > 0 && kv[1] == "") kv <- kv[-1]
  if (length(kv) %% 2 != 0)          kv <- c(kv, "")

  keywords <- stats::setNames(
    as.list(kv[seq(2, length(kv), 2)]),
    kv[seq(1, length(kv), 2)]
  )

  # Must close before fcs_rcpp_read_data opens its own handle
  close(con)
  on.exit(NULL)

  # --- 2. Resolve dimensions and byte offset --------------------------------
  total.events <- as.numeric(keywords[["$TOT"]])
  n.par        <- as.numeric(keywords[["$PAR"]])
  data.st      <- as.numeric(keywords[["$BEGINDATA"]])

  byte.ord <- keywords[["$BYTEORD"]]
  swap     <- !identical(byte.ord, "1,2,3,4")

  # --- 2a. Sanity-check DATA segment size against $TOT/$PAR -----------------
  # Catches header/TEXT corruption or a mismatched $BEGINDATA/$ENDDATA pair
  # before it ever reaches the C++ seek/read
  expected.bytes <- total.events * n.par * 4
  if (!is.null(keywords[["$ENDDATA"]])) {
    data.en <- as.numeric(keywords[["$ENDDATA"]])
    actual.bytes <- data.en - data.st + 1
    if (abs(actual.bytes - expected.bytes) > n.par * 4) {
      stop(
        "DATA segment size implied by $TOT (", total.events, ") x $PAR (",
        n.par, ") x 4 bytes = ", expected.bytes,
        " bytes does not match $BEGINDATA/$ENDDATA (", actual.bytes,
        " bytes). File may be corrupt or truncated: ", fcs.path,
        call. = FALSE
      )
    }
  }

  read.start       <- if (is.null(start.row)) 1 else as.numeric(start.row)
  read.end         <- if (is.null(end.row))   total.events else as.numeric(end.row)
  num.rows.to.read <- read.end - read.start + 1

  # Offset into the data segment for the requested row range
  byte.offset <- data.st + ((read.start - 1) * n.par * 4)

  # --- 2b. Validate DATATYPE/PnB match the reader's float32 assumption ------
  # Both the C++ and pure-R readers hard-code 4-byte IEEE-754 float reads.
  data.type <- keywords[["$DATATYPE"]]
  if (!is.null(data.type) && !identical(data.type, "F")) {
    stop(
      "Unsupported $DATATYPE '", data.type, "' in ", fcs.path,
      ": this reader only supports $DATATYPE/F/ (32-bit float). ",
      "Files with $DATATYPE/D/ or $DATATYPE/I/ require a different code path.",
      call. = FALSE
    )
  }

  pnb.vals <- vapply(seq_len(n.par), function(i) {
    val <- keywords[[paste0("$P", i, "B")]]
    if (is.null(val)) NA_real_ else as.numeric(val)
  }, numeric(1L))

  if (any(!is.na(pnb.vals) & pnb.vals != 32)) {
    bad <- which(!is.na(pnb.vals) & pnb.vals != 32)
    stop(
      "Parameter(s) ", paste(bad, collapse = ", "),
      " in ", fcs.path, " declare $PnB != 32 (mixed bit depths are not ",
      "supported by this reader, which assumes a uniform 32-bit float ",
      "record layout).",
      call. = FALSE
    )
  }

  # --- 3. Resolve the full parameter name list, then an optional subset -----
  all.col.names <- unname(vapply(seq_len(n.par), function(i) {
    val <- keywords[[paste0("$P", i, "N")]]
    if (is.null(val)) paste0("Channel_", i) else val
  }, character(1L)))

  if (is.null(columns)) {
    selected.idx <- NULL
    col.names    <- all.col.names
  } else {
    selected.idx <- match(columns, all.col.names)
    if (anyNA(selected.idx)) {
      stop(
        "Requested column(s) not found in ", fcs.path, ": ",
        paste(columns[is.na(selected.idx)], collapse = ", "),
        call. = FALSE
      )
    }
    col.names <- all.col.names[selected.idx]
  }

  # --- 4. Read binary data via C++ -------------------------------------------
  data.mat <- fcs_rcpp_read_data_test(
    file_path     = fcs.path,
    byte_offset   = byte.offset,
    n_row         = num.rows.to.read,
    n_par         = n.par,
    swap          = swap,
    selected_cols = if (is.null(selected.idx)) NULL else selected.idx - 1L
  )
  colnames(data.mat) <- col.names

  # --- 5. Return ------------------------------------------------------------
  if (return.keywords) list(data = data.mat, keywords = keywords)
  else data.mat
}
