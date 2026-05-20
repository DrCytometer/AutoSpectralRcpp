// fcs_rcpp_write_data.cpp
//
// Writes a pre-built FCS header + TEXT segment + data matrix to disk.
//
// The R side (writeFCS in AutoSpectralRcpp) handles all keyword construction,
// offset resolution, and header formatting — the same logic as the pure-R
// writeFCS.  This function handles only the binary write, which is where the
// memory and performance cost lies.
//
// Key advantage over the R writeBin path:
//   R's writeBin(as.vector(t(mat)), ...) allocates:
//     (1) a transposed double matrix  (~n_row * n_col * 8 bytes)
//     (2) a flat double vector        (~n_row * n_col * 8 bytes)
//     (3) an internal float32 buffer  (~n_row * n_col * 4 bytes)
//   Total transient: ~20 bytes/value on top of the input matrix.
//
//   This function transposes and casts to float32 in a small stack-allocated
//   chunk buffer, writing directly to disk with fwrite().  Peak extra
//   allocation: one chunk (~CHUNK_ROWS * n_col * 4 bytes, default 4 MB).

#include <Rcpp.h>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

using namespace Rcpp;

// Chunk size for the transpose+cast+write loop.
// 131072 rows * 4 bytes/float * n_col: for 183 channels this is ~96 MB per
// chunk, which is well within budget.  Adjust if needed.
static const int CHUNK_ROWS = 131072;

// [[Rcpp::export]]
void fcs_rcpp_write_data(
        const std::string& file_path,
        const std::string& header,        // exactly 58 bytes
        const std::string& text_segment,  // TEXT segment as built by R
        NumericMatrix      data_mat,       // events (rows) x channels (cols)
        bool               swap           // byte-swap for big-endian output
) {
    const long n_row = data_mat.nrow();
    const long n_col = data_mat.ncol();

    // validate header length
    if ((long)header.size() != 58)
        stop("header must be exactly 58 bytes, got " +
             std::to_string(header.size()));

    // open in binary write mode
    FILE* f = std::fopen(file_path.c_str(), "wb");
    if (!f)
        stop("Cannot open file for writing: " + file_path);

    // --- write header (58 bytes) ---
    if (std::fwrite(header.data(), 1, header.size(), f) != header.size()) {
        std::fclose(f);
        stop("Failed to write FCS header to: " + file_path);
    }

    // --- write TEXT segment ---
    if (std::fwrite(text_segment.data(), 1, text_segment.size(), f)
            != text_segment.size()) {
        std::fclose(f);
        stop("Failed to write TEXT segment to: " + file_path);
    }

    // --- write DATA segment: transpose + cast double->float32 + optional swap ---
    // Rcpp NumericMatrix is column-major: data_mat(row, col) = col*n_row + row.
    // FCS requires row-major float32: event0[ch0,ch1,...], event1[ch0,ch1,...].
    // We process CHUNK_ROWS events at a time to bound peak memory.

    const long chunk_rows = std::min((long)CHUNK_ROWS, n_row);
    std::vector<float> buf(chunk_rows * n_col);

    long rows_remaining = n_row;
    long row_offset     = 0;

    while (rows_remaining > 0) {
        long rows_this_chunk = std::min((long)CHUNK_ROWS, rows_remaining);

        // transpose chunk from column-major double into row-major float32
        for (long col = 0; col < n_col; col++) {
            const double* col_ptr = &data_mat(row_offset, col);  // col-major ptr
            for (long r = 0; r < rows_this_chunk; r++) {
                buf[r * n_col + col] = static_cast<float>(col_ptr[r]);
            }
        }

        // byte-swap if writing big-endian (rare but spec-legal)
        if (swap) {
            for (long k = 0; k < rows_this_chunk * n_col; k++) {
                char* p = reinterpret_cast<char*>(&buf[k]);
                std::swap(p[0], p[3]);
                std::swap(p[1], p[2]);
            }
        }

        long n_vals  = rows_this_chunk * n_col;
        long written = static_cast<long>(
            std::fwrite(buf.data(), sizeof(float), n_vals, f)
        );
        if (written != n_vals) {
            std::fclose(f);
            stop("Failed to write data chunk at row " +
                 std::to_string(row_offset) + ": wrote " +
                 std::to_string(written) + " of " +
                 std::to_string(n_vals) + " values");
        }

        row_offset     += rows_this_chunk;
        rows_remaining -= rows_this_chunk;
    }

    // --- write footer (FCS standard: 8 ASCII zeros) ---
    const char footer[8] = {'0','0','0','0','0','0','0','0'};
    if (std::fwrite(footer, 1, 8, f) != 8) {
        std::fclose(f);
        stop("Failed to write FCS footer to: " + file_path);
    }

    std::fclose(f);
}
