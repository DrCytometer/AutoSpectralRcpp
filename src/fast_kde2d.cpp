#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' Fast Kernel Density Estimation in 2D
//'
//' @param x Numeric vector of x-values (typically FSC-A)
//' @param y Numeric vector of y-values (typically SSC-A)
//' @param n Number of grid points for both x and y (recommend 100)
//' @param h Vector of bandwidths for x and y (recommed 0.1*bandwidth.nrd)
//' @param x_limits Limits of the data in the x coordinates
//' @param y_limits Limits of the data in the y coordinates
//' @export
// [[Rcpp::export]]
List fast_kde2d_cpp(NumericVector x, NumericVector y, int n, NumericVector h, NumericVector x_limits, NumericVector y_limits) {
  int nx = x.size();
  if (nx < 2) return List::create(_["x"]=NumericVector(n), _["y"]=NumericVector(n), _["z"]=NumericMatrix(n,n));

  double x_min = x_limits[0], x_max = x_limits[1];
  double y_min = y_limits[0], y_max = y_limits[1];
  double dx = (x_max - x_min) / (n - 1);
  double dy = (y_max - y_min) / (n - 1);

  // Create a slightly larger internal grid to prevent edge spikes
  NumericMatrix z(n, n);

  // 1. Fast Linear Binning (O(N))
  for (int i = 0; i < nx; ++i) {
    double x_pos = (x[i] - x_min) / dx;
    double y_pos = (y[i] - y_min) / dy;
    int ix = (int)std::floor(x_pos);
    int iy = (int)std::floor(y_pos);

    if (ix >= 0 && ix < n - 1 && iy >= 0 && iy < n - 1) {
      double wx = x_pos - ix;
      double wy = y_pos - iy;
      z(ix, iy)         += (1.0 - wx) * (1.0 - wy);
      z(ix + 1, iy)     += wx * (1.0 - wy);
      z(ix, iy + 1)     += (1.0 - wx) * wy;
      z(ix + 1, iy + 1) += wx * wy;
    }
  }

  // 2. Separable Smoothing (O(n^2)) - The fast way
  auto smooth_pass = [&](NumericVector vec, double sigma, double step) {
    int sz = vec.size();
    NumericVector out(sz);
    double s2 = 2.0 * sigma * sigma;
    int radius = (int)std::ceil(3.0 * sigma / step);

    for (int i = 0; i < sz; ++i) {
      double sum = 0, weights = 0;
      for (int j = std::max(0, i - radius); j <= std::min(sz - 1, i + radius); ++j) {
        double dist = (i - j) * step;
        double w = std::exp(-(dist * dist) / s2);
        sum += vec[j] * w;
        weights += w;
      }
      // Ensure edges taper to zero if mass is low
      out[i] = (weights > 0) ? (sum / weights) : 0;
    }
    return out;
  };

  for (int j = 0; j < n; ++j) z(_, j) = smooth_pass(z(_, j), h[0], dx);
  for (int i = 0; i < n; ++i) z(i, _) = smooth_pass(z(i, _), h[1], dy);

  // 3. Normalization: Matches the MASS::kde2d formula
  // Scale = 1 / (N * hx * hy)
  double mass_scale = 1.0 / (nx * dx * dy);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      z(i, j) *= mass_scale;
    }
  }

  NumericVector x_out(n), y_out(n);
  for(int k = 0; k < n; k++) {
    x_out[k] = x_min + k * dx;
    y_out[k] = y_min + k * dy;
  }

  return List::create(_["x"] = x_out, _["y"] = y_out, _["z"] = z);
}
