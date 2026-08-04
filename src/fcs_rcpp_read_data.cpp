//

#include <Rcpp.h>
#include <fstream>
#include <stdexcept>
#include <cstdint>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fcs_rcpp_read_data(
    const std::string& file_path,
    double byte_offset,
    double n_row,
    double n_par,
    bool swap
) {
  if (n_row <= 0) stop("n_row must be positive");
  if (n_par <= 0) stop("n_par must be positive");

  // R doubles are exact integers to 2^53 — plenty of headroom. Cast once,
  // here, to a guaranteed-64-bit type. Never declare these as `long`: it's
  // only 32 bits under the Windows (LLP64) toolchain R ships with, unlike
  // Linux/macOS (LP64), which is what caused the original crash.
  const int64_t byte_offset_i = static_cast<int64_t>(byte_offset);
  const int64_t n_row_i       = static_cast<int64_t>(n_row);
  const int64_t n_par_i       = static_cast<int64_t>(n_par);
  const int64_t n_vals        = n_row_i * n_par_i;

  std::ifstream con(file_path, std::ios::binary | std::ios::ate);
  if (!con.is_open()) stop("Cannot open file: " + file_path);

  // Bounds check against actual file size before seeking/reading — same
  // safeguard as test1.cpp, carried over so we don't lose it by switching
  // to the streaming implementation.
  const int64_t file_size = static_cast<int64_t>(con.tellg());
  const int64_t bytes_needed = static_cast<int64_t>(sizeof(float)) * n_vals;

  if (byte_offset_i < 0 || byte_offset_i > file_size)
    stop("byte_offset (" + std::to_string(byte_offset_i) +
      ") is outside the file (" + std::to_string(file_size) +
      " bytes) — header/TEXT offsets may be corrupt.");

  if (byte_offset_i + bytes_needed > file_size)
    stop("Requested read extends " +
      std::to_string(byte_offset_i + bytes_needed - file_size) +
      " bytes past end of file — $TOT/$PAR/$BEGINDATA may be corrupt "
      "or inconsistent.");

  con.seekg(byte_offset_i, std::ios::beg);
  if (con.fail())
    stop("Failed to seek to offset " + std::to_string(byte_offset_i));

  // Allocate only the final double matrix; stream one row of floats at a
  // time through a small reusable scratch buffer instead of materializing
  // the whole float copy of the dataset alongside it. Peak extra memory is
  // now one row (n_par floats) rather than the full dataset.
  NumericMatrix data_mat(static_cast<int>(n_row_i), static_cast<int>(n_par_i));
  double* out = data_mat.begin();

  std::vector<float> row_buf(static_cast<size_t>(n_par_i));
  const std::streamsize row_bytes =
    static_cast<std::streamsize>(n_par_i) * sizeof(float);

  for (int64_t row = 0; row < n_row_i; row++) {
    con.read(reinterpret_cast<char*>(row_buf.data()), row_bytes);

    if (con.gcount() != row_bytes)
      stop("Short read at row " + std::to_string(row) +
        " — file may be truncated or byte_offset/n_row/n_par are "
        "incorrect. Expected " + std::to_string(row_bytes) +
          " bytes, got " + std::to_string(con.gcount()));

    if (swap) {
      for (int64_t col = 0; col < n_par_i; col++) {
        char* p = reinterpret_cast<char*>(&row_buf[static_cast<size_t>(col)]);
        std::swap(p[0], p[3]);
        std::swap(p[1], p[2]);
      }
    }

    for (int64_t col = 0; col < n_par_i; col++) {
      out[col * n_row_i + row] = static_cast<double>(row_buf[static_cast<size_t>(col)]);
    }
  }

  return data_mat;
}
