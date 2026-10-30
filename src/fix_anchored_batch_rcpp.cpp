// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Exact selection-based median (stats::median() convention for even n),
// through a caller-owned scratch buffer.
static double anchored_median( std::vector<double>& scratch, const double* src, int n ) {

  std::copy( src, src + n, scratch.begin() );

  std::size_t half = static_cast<std::size_t>( n ) / 2;
  std::nth_element( scratch.begin(), scratch.begin() + half, scratch.begin() + n );

  if ( n % 2 == 1 ) return scratch[ half ];

  double upper = scratch[ half ];
  double lower = *std::max_element( scratch.begin(), scratch.begin() + half );
  return ( lower + upper ) / 2.0;
}

static double anchored_scale( const std::vector<double>& r, int n,
                              std::vector<double>& absdev, std::vector<double>& scratch ) {

  const double meps = std::numeric_limits<double>::epsilon();
  double med = anchored_median( scratch, r.data(), n );
  for ( int i = 0; i < n; ++i ) absdev[ i ] = std::fabs( r[ i ] - med );
  return std::max( anchored_median( scratch, absdev.data(), n ) / 0.6745, meps );
}

// Huber IRLS slope of a line through the origin, warm-started at `start`.
// Returns false (slope 0, se NA) when it does not converge.
static bool huber_origin( const double* xs, const double* ys, int n, double start,
                          double k, int max_iter, double tol, double& slope, double& se,
                          std::vector<double>& r, std::vector<double>& absdev,
                          std::vector<double>& scratch ) {

  const double meps = std::numeric_limits<double>::epsilon();
  double b = start;

  for ( int it = 0; it < max_iter; ++it ) {

    for ( int i = 0; i < n; ++i ) r[ i ] = ys[ i ] - b * xs[ i ];
    double scale = anchored_scale( r, n, absdev, scratch );

    double sxx = 0.0, sxy = 0.0;
    for ( int i = 0; i < n; ++i ) {
      double w = std::min( 1.0, k / std::max( std::fabs( r[ i ] / scale ), meps ) );
      sxx += w * xs[ i ] * xs[ i ];
      sxy += w * xs[ i ] * ys[ i ];
    }
    if ( !( sxx > 0.0 ) ) break;

    double b_new = sxy / sxx;
    bool settled = std::fabs( b_new - b ) < tol * std::max( std::fabs( b_new ), 1.0 );
    b = b_new;

    if ( settled ) {
      for ( int i = 0; i < n; ++i ) r[ i ] = ys[ i ] - b * xs[ i ];
      double s2 = anchored_scale( r, n, absdev, scratch );
      double wxx = 0.0;
      for ( int i = 0; i < n; ++i ) {
        double w = std::min( 1.0, k / std::max( std::fabs( r[ i ] / s2 ), meps ) );
        wxx += w * xs[ i ] * xs[ i ];
      }
      slope = b;
      se = s2 / std::sqrt( std::max( wxx, meps ) );
      return true;
    }
  }

  slope = 0.0;
  se = NA_REAL;
  return false;
}

//' Anchored Spillover Slopes, One Source Against Every Target (C++)
//'
//' @description
//' Compiled form of AutoSpectral's `.fix.anchored.r()` (the anchored pair
//' estimator of `fix.my.unmix()`), with the same results. For each target:
//' the source-negative bulk events that are negative for the target, with
//' the current slope removed, give the anchor (their median source and
//' target values); a Huber line through the anchor is fitted to the clearly
//' positive source events that are negative for the target (upper
//' boundary), and after the first pass also above the target's lower
//' boundary (band); negativity is re-decided with the new slope up to
//' `max_mask_passes` times, stopping once the slope moves by less than
//' `mask_tolerance` of itself or leaves `max_coefficient`.
//'
//' The caller selects and subsamples the positive (`xb`, `Yb`, ...) and bulk
//' (`xk`, `Yk`, ...) events, so this draws no random numbers and gives the
//' same result for any `n_threads`.
//'
//' @param xb Numeric vector, the clearly positive events' source abundance.
//' @param Yb,Tb Numeric matrices (positive events x targets), their target
//'   abundance and the target's upper boundary.
//' @param Lb Numeric matrix of the same shape, the target's lower boundary,
//'   or `NULL` for a one-sided cut.
//' @param xk Numeric vector, the bulk events' source abundance.
//' @param Yk,Tk Numeric matrices (bulk events x targets), as `Yb`, `Tb`.
//' @param start_slope Numeric vector (targets), warm starts.
//' @param max_coefficient double, largest accepted coefficient.
//' @param max_mask_passes int, negativity re-decisions after the first fit.
//' @param mask_tolerance double, relative slope change that ends them.
//' @param min_count int, fewest target-negative positive and bulk events.
//' @param k,max_iter,tol Huber tuning constant, IRLS iterations and
//'   tolerance.
//' @param n_threads int, OpenMP threads across targets. Default `1`.
//'
//' @return A list with `slope`, `n` and `se` per target; `slope` and `se`
//'   are `NA` where fewer than `min_count` events were target-negative on
//'   the first pass.
//'
//' @export
// [[Rcpp::export]]
List fix_anchored_batch_rcpp(
    NumericVector xb, NumericMatrix Yb, NumericMatrix Tb, Nullable<NumericMatrix> Lb,
    NumericVector xk, NumericMatrix Yk, NumericMatrix Tk,
    NumericVector start_slope,
    double max_coefficient = 0.5,
    int max_mask_passes    = 3,
    double mask_tolerance  = 0.05,
    int min_count          = 25,
    double k               = 1.345,
    int max_iter           = 100,
    double tol             = 1e-4,
    int n_threads          = 1 ) {

  const int nb = xb.size();
  const int nk = xk.size();
  const int m  = Yb.ncol();

  const bool band = Lb.isNotNull();
  NumericMatrix Lmat;
  if ( band ) Lmat = as<NumericMatrix>( Lb.get() );

  const double* xbp = xb.begin();
  const double* xkp = xk.begin();
  const double* ybp = Yb.begin();
  const double* tbp = Tb.begin();
  const double* lbp = band ? Lmat.begin() : nullptr;
  const double* ykp = Yk.begin();
  const double* tkp = Tk.begin();
  const double* stp = start_slope.begin();

  NumericVector slope_out( m, NA_REAL );
  NumericVector se_out( m, NA_REAL );
  IntegerVector n_out( m, 0 );
  double* slope_ptr = slope_out.begin();
  double* se_ptr    = se_out.begin();
  int*    n_ptr     = n_out.begin();

#ifdef _OPENMP
  omp_set_num_threads( std::max( n_threads, 1 ) );
#endif

#pragma omp parallel
{
  const int cap = std::max( nb, nk );
  std::vector<double> xs( cap ), ys( cap ), r( cap ), absdev( cap ), scratch( cap );

#pragma omp for schedule(dynamic)
  for ( int j = 0; j < m; ++j ) {

    // R matrices are column-major: each target's column is contiguous.
    const double* yb = ybp + static_cast<std::size_t>( j ) * nb;
    const double* tb = tbp + static_cast<std::size_t>( j ) * nb;
    const double* lb = band ? lbp + static_cast<std::size_t>( j ) * nb : nullptr;
    const double* yk = ykp + static_cast<std::size_t>( j ) * nk;
    const double* tk = tkp + static_cast<std::size_t>( j ) * nk;

    double slope_curr = stp[ j ];
    if ( !std::isfinite( slope_curr ) || std::fabs( slope_curr ) > max_coefficient ) slope_curr = 0.0;

    double slope_res = NA_REAL, se_res = NA_REAL;
    int n_res = 0;

    for ( int pass = 0; pass <= max_mask_passes; ++pass ) {

      int nkn = 0;
      for ( int i = 0; i < nk; ++i ) {
        if ( yk[ i ] - slope_curr * xkp[ i ] < tk[ i ] ) { xs[ nkn ] = xkp[ i ]; ys[ nkn ] = yk[ i ]; ++nkn; }
      }

      int nbn = 0;
      for ( int i = 0; i < nb; ++i ) {
        double resid = yb[ i ] - slope_curr * xbp[ i ];
        bool keep = resid < tb[ i ];
        if ( keep && band && pass > 0 ) keep = resid > lb[ i ];
        if ( keep ) ++nbn;
      }

      if ( nkn < min_count || nbn < min_count ) break;

      double x0 = anchored_median( scratch, xs.data(), nkn );
      double y0 = anchored_median( scratch, ys.data(), nkn );

      int t = 0;
      for ( int i = 0; i < nb; ++i ) {
        double resid = yb[ i ] - slope_curr * xbp[ i ];
        bool keep = resid < tb[ i ];
        if ( keep && band && pass > 0 ) keep = resid > lb[ i ];
        if ( keep ) { xs[ t ] = xbp[ i ] - x0; ys[ t ] = yb[ i ] - y0; ++t; }
      }

      double slope_next = 0.0, se = NA_REAL;
      huber_origin( xs.data(), ys.data(), t, slope_curr, k, max_iter, tol, slope_next, se,
                    r, absdev, scratch );

      n_res = t;
      slope_res = slope_next;
      se_res = se;

      bool valid = std::fabs( slope_next ) <= max_coefficient;
      bool settled = std::fabs( slope_next - slope_curr ) <= mask_tolerance * std::fabs( slope_next );
      slope_curr = slope_next;

      if ( !valid || settled ) break;
    }

    slope_ptr[ j ] = slope_res;
    se_ptr[ j ]    = se_res;
    n_ptr[ j ]     = n_res;
  }
} // omp parallel

  return List::create( _[ "slope" ] = slope_out, _[ "n" ] = n_out, _[ "se" ] = se_out );
}
