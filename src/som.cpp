// som.cpp
//
// OpenMP-accelerated Self-Organizing Map (SOM) training and mapping for
// AutoSpectral, replacing FlowSOM::SOM() / EmbedSOM::SOM() as the clustering
// engine behind get.af.spectra() and get.fluor.variants().
//
// som_train_cpp() is a faithful port of the online (sequential) Kohonen SOM
// used by kohonen / FlowSOM::SOM() -- one randomly-drawn training event per
// iteration, annealed learning rate and neighbourhood radius -- kept
// numerically comparable to FlowSOM::SOM() for validation purposes.
//
// Parallelisation strategy:
//   - som_train_batch_cpp() replaces the retired online (per-event)
//     som_train_cpp(). Batch SOM has no per-event state dependency: each
//     epoch (1) assigns every event to its nearest code using the codebook
//     as it stood at the END of the previous epoch (read-only for the whole
//     epoch), then (2) recomputes every code as a neighbourhood-weighted
//     average of the per-node event sums. Synchronisation is needed once
//     per epoch (rlen times total), not once per training event (rlen * n
//     times) -- that was the root problem with parallelising the online
//     trainer (see git history / NEWS for the earlier attempt and why it
//     got slower with more threads).
//   - Ported from EmbedSOM's bsom() (see embedsom_som.cpp in project
//     history), replacing its std::thread chunk-per-thread partitioning
//     with OpenMP work-sharing, consistent with the rest of this file.
//   - Step 1 (nearest-code assignment) is parallelised over events using
//     one accumulator buffer (ncodes x px sums, ncodes counts) PER THREAD,
//     indexed by omp_get_thread_num() -- heap vectors, NOT thread_local
//     (thread_local inside an omp for loop body is unsafe per project
//     convention). Each thread only ever writes its own buffer, so there is
//     no cross-thread false sharing during the assignment loop itself.
//   - Per-thread buffers are reduced into one set of global sums serially
//     after the assignment loop: O(threads * ncodes * px), cheap relative
//     to the O(n * ncodes * px) assignment step it follows.
//   - Step 2 (neighbourhood diffusion update) is restructured from
//     EmbedSOM's source-node-major loop order (for si { for di { ... } }),
//     which accumulates into the same di across si iterations and is NOT
//     safe to parallelise directly, into a destination-node-major loop
//     (for di { for si { ... } }) -- mathematically identical, but each di
//     iteration now only reads the (already-reduced) global sums and
//     writes its own distinct row of the codebook, so it parallelises over
//     di the same way the retired trainer's per-code loops did.
//   - map_data_to_codes_cpp() is unchanged -- still fully data-parallel
//     across events with no barriers in its hot loop.

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cfloat>   // DBL_MAX -- missing from the flowsomlite reference; fixed here
#include <cmath>
#include <vector>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Distance functions (same four options as flowsomlite / kohonen).
// Column-major indexing: element (row, col) is ptr[row + col*nrow].
// ---------------------------------------------------------------------------

inline double dist_eucl( const double* p1, const double* p2, int px, int n, int ncodes ) {
  double xdist = 0.0;
  for ( int j = 0; j < px; ++j ) {
    double tmp = p1[ j * n ] - p2[ j * ncodes ];
    xdist += tmp * tmp;
  }
  return std::sqrt( xdist );
}

inline double dist_manh( const double* p1, const double* p2, int px, int n, int ncodes ) {
  double xdist = 0.0;
  for ( int j = 0; j < px; ++j ) {
    xdist += std::fabs( p1[ j * n ] - p2[ j * ncodes ] );
  }
  return xdist;
}

inline double dist_chebyshev( const double* p1, const double* p2, int px, int n, int ncodes ) {
  double xdist = 0.0;
  for ( int j = 0; j < px; ++j ) {
    double tmp = std::fabs( p1[ j * n ] - p2[ j * ncodes ] );
    if ( tmp > xdist ) xdist = tmp;
  }
  return xdist;
}

inline double dist_cosine( const double* p1, const double* p2, int px, int n, int ncodes ) {
  double nom = 0.0, denom1 = 0.0, denom2 = 0.0;
  for ( int j = 0; j < px; ++j ) {
    double a = p1[ j * n ];
    double b = p2[ j * ncodes ];
    nom    += a * b;
    denom1 += a * a;
    denom2 += b * b;
  }
  return ( -nom / ( std::sqrt( denom1 ) * std::sqrt( denom2 ) ) ) + 1.0;
}

typedef double (*dist_fun)( const double*, const double*, int, int, int );

inline dist_fun get_dist_fun( int dist ) {
  switch ( dist ) {
  case 1:  return &dist_manh;
  case 3:  return &dist_chebyshev;
  case 4:  return &dist_cosine;
  case 2:
  default: return &dist_eucl;
  }
}

// ---------------------------------------------------------------------------
// som_train_cpp(): online Kohonen SOM training
// ---------------------------------------------------------------------------

//' Train a Self-Organizing Map with batch updates (OpenMP-accelerated)
//'
//' Port of EmbedSOM's batch SOM (\code{EmbedSOM::SOM(..., batch = TRUE)}),
//' intended as the clustering engine behind \code{get.af.spectra()} and
//' \code{get.fluor.variants()} without requiring EmbedSOM as a dependency.
//'
//' This is a different algorithm from the retired online \code{som_train_cpp()}
//' (one weighted-average update per epoch over the whole dataset, rather than
//' one small update per randomly-drawn event) and does not reproduce
//' \code{FlowSOM::SOM()}'s per-event codebook trajectory bit-for-bit. It is
//' validated instead against \code{EmbedSOM::SOM(..., batch = TRUE)} output
//' for a fixed seed.
//'
//' @param data Numeric matrix, training events x features (column-major).
//' @param init_codes Numeric matrix, initial codebook (rows = SOM nodes),
//'   same columns as \code{data}. Chosen in R (e.g. a random sample of
//'   \code{data} rows) so that \code{set.seed()} controls it.
//' @param nhbrdist Numeric matrix (ncodes x ncodes), grid neighbourhood
//'   distances.
//' @param radii Numeric vector, one neighbourhood radius per epoch --
//'   \code{length(radii)} is the number of epochs (there is no separate
//'   \code{rlen}). Typically a decreasing schedule built in R, e.g.
//'   \code{seq(radius_start, radius_end, length.out = rlen)}. Must stay
//'   strictly positive: the neighbourhood kernel is
//'   \code{exp(-d^2 / radius^2)}, and a radius of exactly 0 would starve
//'   every node outside the exact best-matching unit (guarded internally
//'   via \code{min_radius}, but avoid relying on that guard).
//' @param dist Integer 1:4, distance function (1 manhattan, 2 euclidean,
//'   3 chebyshev, 4 cosine). Default 2.
//' @param n_threads Integer, OpenMP threads. Default 0 (all available
//'   cores, via \code{omp_get_max_threads()}).
//' @return Numeric matrix, the trained codebook (ncodes x features).
//' @export
// [[Rcpp::export]]
NumericMatrix som_train_batch_cpp(
    NumericMatrix data,
    NumericMatrix init_codes,
    NumericMatrix nhbrdist,
    NumericVector radii,
    int dist      = 2,
    int n_threads = 0
) {
  const int n      = data.nrow();
  const int px     = data.ncol();
  const int ncodes = init_codes.nrow();
  const int epochs = radii.size();

  if ( epochs < 1 )
    stop( "`radii` must have length >= 1 (one entry per epoch)." );

  NumericMatrix codes = clone( init_codes );

  const double* data_ptr = data.begin();
  double*       code_ptr = codes.begin();
  const double* nhbr_ptr = nhbrdist.begin();

  dist_fun distf = get_dist_fun( dist );

#ifdef _OPENMP
  if ( n_threads < 1 ) n_threads = omp_get_max_threads();
  omp_set_num_threads( n_threads );
  const int actual_threads = n_threads;
#else
  const int actual_threads = 1;
#endif

  // Per-thread accumulators, indexed by thread id -- heap vectors, NOT
  // thread_local, so safe to touch from inside #pragma omp for.
  std::vector<double> thread_sums(   (size_t) actual_threads * ncodes * px, 0.0 );
  std::vector<double> thread_counts( (size_t) actual_threads * ncodes,      0.0 );

  std::vector<double> global_sums(   (size_t) ncodes * px );
  std::vector<double> global_counts( ncodes );
  std::vector<double> prev_codes(    (size_t) ncodes * px );

  const double min_radius = 1e-10;

  for ( int epoch = 0; epoch < epochs; ++epoch ) {

    R_CheckUserInterrupt(); // safe here -- outside any parallel region

    std::fill( thread_sums.begin(),   thread_sums.end(),   0.0 );
    std::fill( thread_counts.begin(), thread_counts.end(), 0.0 );

    // --- Step 1: assign every event to its nearest code (parallel over events) ---
#pragma omp parallel default(shared)
{
  int tid = 0;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#endif
  double* my_sums   = &thread_sums[   (size_t) tid * ncodes * px ];
  double* my_counts = &thread_counts[ (size_t) tid * ncodes ];

#pragma omp for schedule(static)
  for ( int ev = 0; ev < n; ++ev ) {
    int    nearest = 0;
    double neard   = distf( &data_ptr[ ev ], &code_ptr[ 0 ], px, n, ncodes );
    for ( int cd = 1; cd < ncodes; ++cd ) {
      double tmp = distf( &data_ptr[ ev ], &code_ptr[ cd ], px, n, ncodes );
      if ( tmp < neard ) { neard = tmp; nearest = cd; }
    }
    my_counts[ nearest ] += 1.0;
    for ( int j = 0; j < px; ++j ) {
      my_sums[ nearest + j * ncodes ] += data_ptr[ ev + j * n ];
    }
  }
} // end omp parallel -- barrier guarantees every thread's buffer is complete

    // --- reduce per-thread buffers (serial, O(threads * ncodes * px)) ---
    std::fill( global_sums.begin(),   global_sums.end(),   0.0 );
    std::fill( global_counts.begin(), global_counts.end(), 0.0 );
    for ( int t = 0; t < actual_threads; ++t ) {
      const double* ts = &thread_sums[   (size_t) t * ncodes * px ];
      const double* tc = &thread_counts[ (size_t) t * ncodes ];
      for ( size_t idx = 0; idx < (size_t) ncodes * px; ++idx ) global_sums[ idx ] += ts[ idx ];
      for ( int cd = 0; cd < ncodes; ++cd )                     global_counts[ cd ] += tc[ cd ];
    }

    std::copy( code_ptr, code_ptr + (size_t) ncodes * px, prev_codes.begin() );

    // --- Step 2: neighbourhood-weighted diffusion update (parallel over destination node) ---
    const double radius        = std::max( min_radius, (double) radii[ epoch ] );
    const double inv_sq_radius = -1.0 / ( radius * radius );

#pragma omp parallel for schedule(static)
    for ( int di = 0; di < ncodes; ++di ) {
      std::vector<double> row_sum( px, 0.0 );
      double row_weight = 0.0;
      for ( int si = 0; si < ncodes; ++si ) {
        double d = nhbr_ptr[ si + ncodes * di ];
        double w = std::exp( d * d * inv_sq_radius );
        row_weight += w * global_counts[ si ];
        for ( int j = 0; j < px; ++j )
          row_sum[ j ] += global_sums[ si + j * ncodes ] * w;
      }
      if ( row_weight > 0.0 ) {
        for ( int j = 0; j < px; ++j )
          code_ptr[ di + j * ncodes ] = row_sum[ j ] / row_weight;
      } else {
        // No node contributed within kernel range -- hold previous position
        // rather than collapsing to zero (matches EmbedSOM's oldkoho fallback).
        for ( int j = 0; j < px; ++j )
          code_ptr[ di + j * ncodes ] = prev_codes[ di + j * ncodes ];
      }
    }
  }

  return codes;
}

// ---------------------------------------------------------------------------
// map_data_to_codes_cpp(): map new data onto a trained codebook
// Fully data-parallel across events -- no shared mutable state.
// ---------------------------------------------------------------------------

//' Map data onto a trained SOM codebook (OpenMP-accelerated)
//'
//' For each row of \code{data}, finds the nearest code in \code{codes}.
//' Fully data-parallel across events.
//'
//' @param data Numeric matrix, events x features (same column order as
//'   \code{codes}).
//' @param codes Numeric matrix, SOM codebook (nodes x features).
//' @param dist Integer 1:4, distance function (see \code{som_train_cpp}).
//' @param n_threads Integer, OpenMP threads. Default 1.
//' @return Numeric matrix, 2 columns: nearest code id (1-based) and distance
//'   to that code, one row per event.
//' @export
// [[Rcpp::export]]
NumericMatrix map_data_to_codes_cpp(
     NumericMatrix data,
     NumericMatrix codes,
     int dist      = 2,
     int n_threads = 1
) {
   const int nd     = data.nrow();
   const int p      = data.ncol();
   const int ncodes = codes.nrow();

   const double* data_ptr = data.begin();
   const double* code_ptr = codes.begin();

   dist_fun distf = get_dist_fun( dist );

   NumericMatrix result( nd, 2 );

#ifdef _OPENMP
   if ( n_threads < 1 ) n_threads = 1;
   omp_set_num_threads( n_threads );
#endif

   // No R_CheckUserInterrupt() here: it is not safe to call from inside a
   // parallel region.
#pragma omp parallel for schedule(static)
   for ( int i = 0; i < nd; ++i ) {
     int    minid   = 0;
     double mindist = DBL_MAX;
     for ( int cd = 0; cd < ncodes; ++cd ) {
       double tmp = distf( &data_ptr[ i ], &code_ptr[ cd ], p, nd, ncodes );
       if ( tmp < mindist ) { mindist = tmp; minid = cd; }
     }
     result( i, 0 ) = minid + 1;
     result( i, 1 ) = mindist;
   }

   return result;
}
