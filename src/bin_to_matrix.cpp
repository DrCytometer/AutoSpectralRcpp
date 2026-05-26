#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' Bin a spectral expression matrix into a 2D count matrix
//'
//' Takes an n_events x n_channels matrix of already biexp-transformed values
//' and bins them into an n_y x n_x count matrix. Loops column-major to match
//' R's memory layout and avoid cache misses.
//'
//' @param data NumericMatrix of biexp-transformed expression values,
//'   n_events rows x n_channels (n_x) columns.
//' @param y_breaks NumericVector of length n_y + 1, uniformly spaced break
//'   points spanning the y-axis range. Produced by seq(y_min, y_max,
//'   length.out = n_y + 1).
//' @param n_y Integer number of y bins.
//'
//' @return IntegerMatrix of counts, n_y rows x n_x columns.
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix bin_matrix_cpp(
     NumericMatrix data,
     NumericVector y_breaks,
     int n_y
 ) {
   int n_events = data.nrow();
   int n_x      = data.ncol();
   IntegerMatrix mat( n_y, n_x );

   double y_min   = y_breaks[0];
   double y_max   = y_breaks[y_breaks.size() - 1];
   double y_range = y_max - y_min;

   // loop column-major: R matrices are column-major so this is cache-friendly
   for ( int col = 0; col < n_x; col++ ) {
     for ( int i = 0; i < n_events; i++ ) {
       double v = data( i, col );
       if ( v < y_min || v > y_max ) continue;
       int yi = (int)( ( v - y_min ) / y_range * ( n_y - 1 ) );
       if ( yi < 0 )     yi = 0;
       if ( yi >= n_y )  yi = n_y - 1;
       mat( yi, col )++;
     }
   }

   return mat;
 }
