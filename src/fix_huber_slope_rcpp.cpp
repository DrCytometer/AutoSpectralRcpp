#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// Exact selection-based median, matching stats::median.default() on a
// vector with no missing values. Partial-sorts (std::nth_element, O(n)
// average) rather than fully sorting, since only the middle value(s) are
// ever read - the same shortcut R's own median() takes internally.
//
// Writes into `scratch` (copying `n` elements from `src` first, since
// nth_element must reorder its input and the caller still needs `src` in
// its original event order afterwards) rather than allocating a fresh
// buffer per call, so repeated calls across IRLS iterations reuse one
// heap allocation instead of paying malloc/free every iteration.
static double median_into( std::vector<double>& scratch,
                           const double* src, int n ) {

  std::copy( src, src + n, scratch.begin() );

  std::size_t half = static_cast<std::size_t>( n ) / 2;
  std::nth_element( scratch.begin(), scratch.begin() + half,
                    scratch.begin() + n );

  if ( n % 2 == 1 ) return scratch[ half ];

  double upper = scratch[ half ];
  double lower = *std::max_element( scratch.begin(), scratch.begin() + half );
  return ( lower + upper ) / 2.0;
}

//' Huber-Weighted IRLS Slope (C++)
//'
//' @description
//' Line-for-line port of AutoSpectral's `.fix.huber.slope()`
//' (fix_my_unmix.R): a Huber-weighted, iteratively reweighted least squares
//' fit of `y` on `x`, one predictor with an intercept. Every step is the
//' same closed-form weighted mean/variance update the R version uses - no
//' BLAS call here either, only the iteration loop and the two `median()`
//' calls per iteration moving out of the R interpreter, which is what
//' `.fix.envelope.slope()`'s call volume (every panel pair, every mask
//' pass, every spillover-estimator iteration, every inner-loop pass) makes
//' worth paying for. All scratch buffers are allocated once, outside the
//' iteration loop, and reused every iteration.
//'
//' On non-convergence this returns a zero slope rather than falling back to
//' OLS, matching the R version: a coefficient this function cannot pin down
//' is exactly the case the caller needs to see as untrustworthy.
//'
//' @param x,y Numeric vectors, the predictor and response. Must be the same
//'   length.
//' @param k Numeric, the Huber tuning constant. Default `1.345`.
//' @param max_iter Integer, maximum IRLS iterations. Default `100`.
//' @param tol Numeric, relative coefficient-change convergence tolerance.
//'   Default `1e-4`.
//' @param start Optional numeric vector `c(intercept, slope)`, a warm
//'   start. When supplied, the first iteration's weights come from its
//'   residuals instead of OLS weights. Default `NULL`.
//'
//' @return Numeric vector `c(intercept, slope)`.
//'
//' @export
// [[Rcpp::export]]
NumericVector fix_huber_slope_rcpp(
    NumericVector x, NumericVector y,
    double k               = 1.345,
    int max_iter            = 100,
    double tol              = 1e-4,
    Rcpp::Nullable<NumericVector> start = R_NilValue ) {

  int n = x.size();
  if ( n < 3 ) return NumericVector::create( 0.0, 0.0 );

  const double meps = std::numeric_limits<double>::epsilon();
  const double* xp = x.begin();
  const double* yp = y.begin();

  std::vector<double> w( n, 1.0 ), w_next( n );
  std::vector<double> resid( n ), absdev( n ), scratch( n );
  double coef0 = 0.0, coef1 = 0.0;

  if ( start.isNotNull() ) {

    NumericVector s( start );
    coef0 = s[ 0 ];
    coef1 = s[ 1 ];

    for ( int i = 0; i < n; ++i )
      resid[ i ] = yp[ i ] - ( coef0 + coef1 * xp[ i ] );

    double med0 = median_into( scratch, resid.data(), n );
    for ( int i = 0; i < n; ++i )
      absdev[ i ] = std::fabs( resid[ i ] - med0 );

    double scale = median_into( scratch, absdev.data(), n ) / 0.6745;
    scale = std::max( scale, meps );

    for ( int i = 0; i < n; ++i ) {
      double u = std::fabs( resid[ i ] / scale );
      w[ i ] = std::min( 1.0, k / std::max( u, meps ) );
    }
  }

  bool converged = false;
  int iters_used = max_iter;

  for ( int iter = 0; iter < max_iter; ++iter ) {

    double sw = 0.0, swx = 0.0, swy = 0.0;
    for ( int i = 0; i < n; ++i ) {
      sw  += w[ i ];
      swx += w[ i ] * xp[ i ];
      swy += w[ i ] * yp[ i ];
    }

    if ( sw <= 0.0 ) return NumericVector::create( 0.0, 0.0 );

    double x_mean = swx / sw;
    double y_mean = swy / sw;

    double sxx = 0.0, sxy = 0.0;
    for ( int i = 0; i < n; ++i ) {
      double dx = xp[ i ] - x_mean;
      sxx += w[ i ] * dx * dx;
      sxy += w[ i ] * dx * ( yp[ i ] - y_mean );
    }

    if ( sxx <= 0.0 ) return NumericVector::create( 0.0, 0.0 );

    double slope     = sxy / sxx;
    double intercept = y_mean - slope * x_mean;

    for ( int i = 0; i < n; ++i )
      resid[ i ] = yp[ i ] - ( intercept + slope * xp[ i ] );

    double med0 = median_into( scratch, resid.data(), n );
    for ( int i = 0; i < n; ++i )
      absdev[ i ] = std::fabs( resid[ i ] - med0 );

    double scale = median_into( scratch, absdev.data(), n ) / 0.6745;
    scale = std::max( scale, meps );

    for ( int i = 0; i < n; ++i ) {
      double u = std::fabs( resid[ i ] / scale );
      w_next[ i ] = std::min( 1.0, k / std::max( u, meps ) );
    }

    double d0 = std::fabs( intercept - coef0 );
    double d1 = std::fabs( slope     - coef1 );
    double moved_scale = std::max(
      std::max( std::fabs( intercept ), std::fabs( slope ) ), 1.0 );
    bool moved = std::max( d0, d1 ) < tol * moved_scale;

    coef0 = intercept;
    coef1 = slope;
    std::swap( w, w_next );   // pointer swap, not an O(n) copy

    if ( moved ) { converged = true; iters_used = iter + 1; break; }
  }

  NumericVector out = converged ?
  NumericVector::create( coef0, coef1 ) :
    NumericVector::create( 0.0, 0.0 );

  // Diagnostic only - not read by .fix.huber.slope(), does not change what
  // the function returns to any existing caller. Check with
  // attr(fix_huber_slope_rcpp(x, y, ...), "iters").
  out.attr( "iters" )     = iters_used;
  out.attr( "converged" ) = converged;

  return out;
}
