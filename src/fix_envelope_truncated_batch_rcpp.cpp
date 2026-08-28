// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Identical to the helper in fix_huber_slope_rcpp.cpp: exact selection-based
// median via std::nth_element, writing into a caller-owned scratch buffer so
// repeated calls (every IRLS iteration, every mask pass, every target) reuse
// one heap allocation instead of paying malloc/free each time.
static double median_into( std::vector<double>& scratch, const double* src, int n ) {

  std::copy( src, src + n, scratch.begin() );

  std::size_t half = static_cast<std::size_t>( n ) / 2;
  std::nth_element( scratch.begin(), scratch.begin() + half, scratch.begin() + n );

  if ( n % 2 == 1 ) return scratch[ half ];

  double upper = scratch[ half ];
  double lower = *std::max_element( scratch.begin(), scratch.begin() + half );
  return ( lower + upper ) / 2.0;
}

// Huber-weighted IRLS on a fixed, already-selected (xs, ys) pair, into
// caller-owned scratch buffers sized to the caller's maximum working-set
// width. Same closed-form update as fix_huber_slope_rcpp(), extracted so it
// can be called once per mask pass per target without re-deriving it.
static double huber_slope_on(
    const double* xs, const double* ys, int ni, double start_intercept,
    double start_slope, double k, int max_iter, double tol,
    std::vector<double>& w, std::vector<double>& w_next,
    std::vector<double>& resid, std::vector<double>& absdev,
    std::vector<double>& scratch, bool* converged_out ) {

  const double meps = std::numeric_limits<double>::epsilon();

  double coef0 = start_intercept, coef1 = start_slope;

  for ( int i = 0; i < ni; ++i )
    resid[ i ] = ys[ i ] - ( coef0 + coef1 * xs[ i ] );

  double med0 = median_into( scratch, resid.data(), ni );
  for ( int i = 0; i < ni; ++i ) absdev[ i ] = std::fabs( resid[ i ] - med0 );
  double scale0 = std::max( median_into( scratch, absdev.data(), ni ) / 0.6745, meps );

  for ( int i = 0; i < ni; ++i ) {
    double u = std::fabs( resid[ i ] / scale0 );
    w[ i ] = std::min( 1.0, k / std::max( u, meps ) );
  }

  bool converged = false;

  for ( int iter = 0; iter < max_iter; ++iter ) {

    double sw = 0.0, swx = 0.0, swy = 0.0;
    for ( int i = 0; i < ni; ++i ) {
      sw  += w[ i ];
      swx += w[ i ] * xs[ i ];
      swy += w[ i ] * ys[ i ];
    }
    if ( sw <= 0.0 ) { *converged_out = false; return 0.0; }

    double x_mean = swx / sw, y_mean = swy / sw;

    double sxx = 0.0, sxy = 0.0;
    for ( int i = 0; i < ni; ++i ) {
      double dx = xs[ i ] - x_mean;
      sxx += w[ i ] * dx * dx;
      sxy += w[ i ] * dx * ( ys[ i ] - y_mean );
    }
    if ( sxx <= 0.0 ) { *converged_out = false; return 0.0; }

    double slope = sxy / sxx;
    double intercept = y_mean - slope * x_mean;

    for ( int i = 0; i < ni; ++i )
      resid[ i ] = ys[ i ] - ( intercept + slope * xs[ i ] );

    double med = median_into( scratch, resid.data(), ni );
    for ( int i = 0; i < ni; ++i ) absdev[ i ] = std::fabs( resid[ i ] - med );
    double scale = std::max( median_into( scratch, absdev.data(), ni ) / 0.6745, meps );

    for ( int i = 0; i < ni; ++i ) {
      double u = std::fabs( resid[ i ] / scale );
      w_next[ i ] = std::min( 1.0, k / std::max( u, meps ) );
    }

    double d0 = std::fabs( intercept - coef0 );
    double d1 = std::fabs( slope     - coef1 );
    double moved_scale = std::max( std::max( std::fabs( intercept ), std::fabs( slope ) ), 1.0 );
    bool moved = std::max( d0, d1 ) < tol * moved_scale;

    coef0 = intercept;
    coef1 = slope;
    std::swap( w, w_next );

    if ( moved ) { converged = true; break; }
  }

  *converged_out = converged;
  return converged ? coef1 : 0.0;
}

//' Batched Truncated Spillover Slope, One Source Against Every Target (C++)
//'
//' @description
//' Port of the truncated-estimator half of `.fix.envelope.slope()` --
//' `select.negative()` plus the `max.mask.passes` refinement loop plus
//' `.fix.huber.slope()` -- run for every target sharing one source's
//' abundance in a single call, instead of once per (source, target) pair.
//'
//' Every working buffer (`w`, `resid`, the selection-index scratch space) is
//' allocated once, sized to `n`, and reused across every target and every
//' mask pass; only the `ni <= n` elements actually selected for a given
//' target's current pass are touched. This is the piece that a batched pure-R
//' implementation cannot get right no matter how it is vectorised: R's
//' `sweep()` / `outer()` return fresh objects every call, so a version built
//' from them reallocates matrices sized to the full working set on every one
//' of up to `targets x mask.passes x irls.iterations` steps, several times
//' over. Batching only pays off once that allocation, and the two
//' `stats::median()` calls' R-level dispatch per iteration, are removed
//' entirely -- which needs compiled code, not a different R vectorisation.
//'
//' Reproduces `select.negative()`'s bulk cap: every bright (source-positive)
//' event is kept regardless of file size, and the origin-hugging bulk below
//' the source threshold is subsampled down to `max_truncated_events` events
//' total once the negative-selected population exceeds it. The bulk sits at
//' the origin and carries no leverage on the fitted slope, but left uncapped
//' it sets the Huber M-estimator's iterative MAD-based scale almost
//' entirely once it outnumbers the bright tail by orders of magnitude --
//' diluting the tail's relative influence and biasing the fitted slope
//' toward zero. Each target gets its own `std::mt19937` stream, seeded from
//' `seed` and reused across passes (not reseeded per pass), so the sampling
//' is reproducible under R's `set.seed()` and identical regardless of
//' `n_threads` -- a target's subsample depends only on its own seed, never
//' on which thread happened to process it. Not bit-identical to
//' `.fix.envelope.slope()`, which draws its bulk subsample from a single
//' shared RNG stream advanced pair by pair, but statistically equivalent to
//' it.
//'
//' @param x Numeric vector, length `n`, the shared source abundance.
//' @param Y Numeric matrix, `n x m`, one column per target.
//' @param Threshold_target Numeric matrix, `n x m`, per-event per-target
//'   positivity boundary.
//' @param threshold_source Numeric vector, length `n`, the source's own
//'   per-event positivity boundary. Recycle to length `n` on the R side
//'   before calling.
//' @param start_slope Numeric vector, length `m`, per-target warm starts.
//' @param seed Integer vector, length `m`, one RNG seed per target, drawn
//'   from R's ambient RNG stream so a caller's `set.seed()` covers this
//'   call.
//' @param max_coefficient double, largest accepted residual coefficient.
//'   Default `0.2`.
//' @param max_mask_passes int, refinement passes after the initial fit.
//'   Default `3`.
//' @param mask_tolerance double, relative change below which a pass is
//'   treated as settled and refinement stops early. Default `0.05`.
//' @param min_events int, minimum selected events for a fit. Default `200`.
//' @param max_truncated_events int, cap on the events used for the robust
//'   fit. Everything above the source threshold is kept and the negative
//'   bulk below it is subsampled. Default `20000`.
//' @param k double, Huber tuning constant. Default `1.345`.
//' @param max_iter int, maximum IRLS iterations per fit. Default `100`.
//' @param tol double, IRLS relative coefficient-change tolerance. Default
//'   `1e-4`.
//' @param n_threads int, OpenMP threads to split targets across. Default
//'   `1` (serial). Every target's fit reads only its own column of `Y` and
//'   `Threshold_target` and writes only its own output slot, so this is
//'   embarrassingly parallel across targets -- but only raise it when this
//'   call is not itself already running inside another parallel context
//'   (`mclapply()` over samples, say). Oversubscribing cores that way is
//'   the same failure mode `RhpcBLASctl::blas_set_num_threads(1)` already
//'   guards against for BLAS inside forked workers elsewhere in this
//'   package -- one source of truth for how many cores are in play at a
//'   time, not two competing ones.
//'
//' @return A list with `slope` (numeric, length `m`, `NA` where the selected
//'   population never reached `min_events`) and `n` (integer, length `m`,
//'   the final pass's selected event count, for `span`/diagnostic use back
//'   in R).
//'
//' @export
// [[Rcpp::export]]
List fix_envelope_truncated_batch_rcpp(
    NumericVector x, NumericMatrix Y, NumericMatrix Threshold_target,
    NumericVector threshold_source, NumericVector start_slope,
    IntegerVector seed,
    double max_coefficient      = 0.2,
    int max_mask_passes         = 3,
    double mask_tolerance       = 0.05,
    int min_events              = 200,
    int max_truncated_events    = 20000,
    double k                    = 1.345,
    int max_iter                = 100,
    double tol                  = 1e-4,
    int n_threads                = 1 ) {

  int n = x.size();
  int m = Y.ncol();

  const double* xp  = x.begin();
  const double* tsp = threshold_source.begin();

  NumericVector slope_out( m, NA_REAL );
  IntegerVector n_out( m, 0 );

  // Raw pointers into already-allocated REALSXP/INTSXP storage. Writing
  // through these from worker threads touches no R API and triggers no
  // allocation or garbage collection, which is what makes it safe -- R's
  // own API is not thread-safe, but a plain write to a distinct element of
  // memory that already exists is just a memory write.
  double* slope_ptr = slope_out.begin();
  int*    n_ptr      = n_out.begin();

#ifdef _OPENMP
  omp_set_num_threads( n_threads );
#endif

  // Every target's fit depends only on its own column of Y and
  // Threshold_target and the shared x, and writes only its own output
  // slot, so this parallelises directly across targets with no
  // cross-target dependency to manage. schedule(dynamic) rather than the
  // implicit static split, because how many mask passes and IRLS
  // iterations a target needs before settling varies target to target --
  // a static up-front division can leave threads idle waiting on whichever
  // one drew the slow targets.
#pragma omp parallel
{
  // Thread-local: each thread allocates these once, on entry to the
  // parallel region, and reuses them for every target and every pass it
  // handles -- the same "allocate once, reuse" rule the serial version
  // followed, now applied per thread rather than once globally.
  std::vector<int>    idx( n ), bright_idx( n ), bulk_idx( n );
  std::vector<double> xs( n ), ys( n );
  std::vector<double> w( n ), w_next( n ), resid( n ), absdev( n ), scratch( n );

#pragma omp for schedule(dynamic)
  for ( int j = 0; j < m; ++j ) {

    double slope_curr = start_slope[ j ];
    if ( !R_finite( slope_curr ) || std::fabs( slope_curr ) > max_coefficient )
      slope_curr = 0.0;

    double slope_result = NA_REAL;
    int    n_result      = 0;

    // One generator per target, seeded from R and reused across every pass
    // so successive passes' bulk subsamples are independent draws rather
    // than repeats of the same one -- the same behaviour repeated calls to
    // R's sample() give the scalar estimator across its own mask passes.
    std::mt19937 gen( static_cast<unsigned int>( seed[ j ] ) );

    for ( int pass = 0; pass <= max_mask_passes; ++pass ) {

      // Every target-negative event is bucketed by source positivity as
      // the scan runs, so the bright/bulk split select.negative() needs is
      // available with no second pass over the selection.
      int n_bright = 0, n_bulk = 0;
      for ( int i = 0; i < n; ++i ) {
        if ( Y( i, j ) - slope_curr * xp[ i ] < Threshold_target( i, j ) ) {
          if ( xp[ i ] > tsp[ i ] ) bright_idx[ n_bright++ ] = i;
          else                      bulk_idx[ n_bulk++ ]    = i;
        }
      }

      int ni_total = n_bright + n_bulk;
      int ni;

      // Mirrors select.negative(): every bright (source-positive) event is
      // kept regardless of file size, and the origin-hugging bulk is
      // subsampled once the total exceeds max_truncated_events.
      if ( ni_total <= max_truncated_events ) {

        ni = 0;
        for ( int t = 0; t < n_bright; ++t ) idx[ ni++ ] = bright_idx[ t ];
        for ( int t = 0; t < n_bulk;    ++t ) idx[ ni++ ] = bulk_idx[ t ];

      } else {

        int n_bulk_keep = std::max( max_truncated_events - n_bright, min_events );
        n_bulk_keep = std::min( n_bulk_keep, n_bulk );

        // Partial Fisher-Yates over the first n_bulk_keep slots only, so
        // the cost is O(n_bulk_keep), not O(n_bulk).
        for ( int t = 0; t < n_bulk_keep; ++t ) {
          std::uniform_int_distribution<int> pick( t, n_bulk - 1 );
          int r = pick( gen );
          std::swap( bulk_idx[ t ], bulk_idx[ r ] );
        }

        ni = 0;
        for ( int t = 0; t < n_bright; ++t )    idx[ ni++ ] = bright_idx[ t ];
        for ( int t = 0; t < n_bulk_keep; ++t ) idx[ ni++ ] = bulk_idx[ t ];
      }

      if ( ni < min_events ) {
        if ( pass == 0 ) { slope_result = NA_REAL; n_result = ni; }
        break;
      }

      for ( int t = 0; t < ni; ++t ) {
        xs[ t ] = xp[ idx[ t ] ];
        ys[ t ] = Y( idx[ t ], j );
      }

      bool converged = false;
      double slope_next = huber_slope_on(
        xs.data(), ys.data(), ni, 0.0, slope_curr, k, max_iter, tol,
        w, w_next, resid, absdev, scratch, &converged );

      // huber_slope_on() returns 0 on non-convergence, matching
      // .fix.huber.slope()'s "non-convergence is a zero slope" contract --
      // but the caller (.fix.envelope.slope()) treats that as a genuine
      // fitted value, not as "no fit", so it is not remapped to NA here.
      n_result     = ni;
      slope_result = slope_next;

      bool valid = std::fabs( slope_next ) <= max_coefficient;
      bool settled = std::fabs( slope_next - slope_curr ) <=
        mask_tolerance * std::fabs( slope_next );

      slope_curr = slope_next;

      if ( !valid || settled ) break;
    }

    slope_ptr[ j ] = slope_result;
    n_ptr[ j ]     = n_result;
  }
} // omp parallel

  return List::create( _[ "slope" ] = slope_out, _[ "n" ] = n_out );
}
