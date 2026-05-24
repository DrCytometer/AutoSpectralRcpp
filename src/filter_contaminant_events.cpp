// filter_contaminant_events.cpp
//
// Fast per-event contaminant filter for get.af.spectra().
//
// For each row (event) of `event_mat`, computes cosine similarity against
// every row of `spectra_mat` and returns a logical vector that is TRUE when
// the maximum similarity across all spectra is BELOW `threshold` (i.e. the
// event is kept).  Events that match any fluorophore signature too closely
// are flagged as contaminants and excluded.
//
// Uses RcppParallel when available; falls back to a sequential loop otherwise.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Filter contaminant events by cosine similarity
//'
//' For each event (row) in \code{event_mat}, computes cosine similarity
//' against every spectrum (row) in \code{spectra_mat}.  Returns a logical
//' vector that is \code{TRUE} for events whose maximum cosine similarity to
//' any spectrum is strictly below \code{threshold}.
//'
//' @param event_mat   Numeric matrix; events x detectors (the unstained
//'   expression data, restricted to the spectral channels).
//' @param spectra_mat Numeric matrix; fluorophores x detectors (the reference
//'   spectra, AF rows already removed, L-inf normalised).
//' @param threshold   Numeric scalar.  Events with max cosine similarity
//'   \eqn{\geq} \code{threshold} are flagged as contaminants.
//' @return A logical vector of length \code{nrow(event_mat)}.  \code{TRUE}
//'   means the event is clean (keep); \code{FALSE} means it is a likely
//'   contaminant (remove).
//' @export
// [[Rcpp::export]]
LogicalVector filter_contaminant_events_cpp(
    NumericMatrix event_mat,
    NumericMatrix spectra_mat,
    double        threshold
) {
  const int n_events  = event_mat.nrow();
  const int n_fluors  = spectra_mat.nrow();
  const int n_cols    = event_mat.ncol();

  // Pre-compute L2 norm for each spectrum row
  NumericVector spec_norms( n_fluors );
  for ( int f = 0; f < n_fluors; ++f ) {
    double ss = 0.0;
    for ( int c = 0; c < n_cols; ++c ) {
      double v = spectra_mat( f, c );
      ss += v * v;
    }
    spec_norms[ f ] = std::sqrt( ss ) + 1e-9;
  }

  LogicalVector keep( n_events, true );

  for ( int i = 0; i < n_events; ++i ) {

    // L2 norm of this event
    double ev_norm = 0.0;
    for ( int c = 0; c < n_cols; ++c ) {
      double v = event_mat( i, c );
      ev_norm += v * v;
    }
    ev_norm = std::sqrt( ev_norm ) + 1e-9;

    // cosine similarity vs each spectrum — exit early if threshold exceeded
    for ( int f = 0; f < n_fluors; ++f ) {
      double dot = 0.0;
      for ( int c = 0; c < n_cols; ++c ) {
        dot += event_mat( i, c ) * spectra_mat( f, c );
      }
      double cs = dot / ( ev_norm * spec_norms[ f ] );
      if ( cs >= threshold ) {
        keep[ i ] = false;
        break;           // no need to check remaining spectra
      }
    }
  }

  return keep;
}
