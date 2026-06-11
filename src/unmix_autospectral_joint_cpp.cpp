#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// Precomputed data for a single optimisable fluorophore
// ---------------------------------------------------------------------------
struct FluorPrecomp {
  bool    active     = false;
  int     master_idx = -1;

  mat  v_mats;      // n_variants x D  – candidate spectra
  mat  delta;       // n_variants x D  – (v_mats - master_row), pre-centred

  // Leakage prediction via the other-fluorophore pseudoinverse
  mat  v_lib;       // (F-1) x n_variants – how each variant leaks into others
  vec  w_leakage;   // F-1               – per-other-fluor cov-propagated weights

  // Rank-1 residual update helpers
  mat  r_lib;       // D x n_variants  – residual component of each variant delta
  vec  r_dots;      // n_variants      – dot(r_lib(:,v), r_lib(:,v))

  uword n_variants = 0;
};

// [[Rcpp::export]]
arma::mat unmix_autospectral_joint_cpp(
    arma::mat              raw_data_in,
    const arma::mat&       spectra,
    const arma::mat&       af_spectra,
    const CharacterVector& fluor_names,
    const arma::vec&       pos_thresholds,
    const List&            variants_list,
    const List&            delta_list,
    int                    n_passes    = 2,
    int                    n_threads   = 1,
    bool                   cell_weight = false,
    double                 noise_floor = 1.0
) {
  mat raw_data = raw_data_in.t();   // D x N
  const uword N   = raw_data.n_cols;
  const uword F   = spectra.n_rows;
  const uword D   = spectra.n_cols;
  const uword nAF = af_spectra.n_rows;

  // =========================================================================
  // SECTION 1 – Global pre-computations
  // =========================================================================

  // Global pseudoinverse P = (S S^T)^{-1} S,  shape F x D
  const mat P = solve(spectra * spectra.t(), spectra);

  // AF helpers
  const mat v_lib_af  = P * af_spectra.t();                        // F   x nAF
  const mat r_lib_af  = af_spectra.t() - spectra.t() * v_lib_af;  // D   x nAF

  vec r_dots_af(nAF);
  for (uword j = 0; j < nAF; ++j)
    r_dots_af[j] = std::max(dot(r_lib_af.col(j), r_lib_af.col(j)), 1e-10);

  // Covariance-propagated AF weights (force evaluation with mat() to avoid
  // Armadillo lazy-expression issues on the triple product).
  const mat af_cov_mat = mat(P * cov(af_spectra) * P.t());
  const vec w_af = sqrt(abs(af_cov_mat.diag())) + 1e-8;

  // Determine whether fluorophore variant optimisation is requested.
  // An empty variants_list signals AF-only mode: skip Section 2 setup and
  // the multi-pass candidate loop entirely.
  const bool af_only = (variants_list.size() == 0);

  // =========================================================================
  // SECTION 2 – Per-fluorophore pre-computations  (skipped in AF-only mode)
  // =========================================================================

  std::vector<std::string> cpp_names = as<std::vector<std::string>>(fluor_names);
  std::map<std::string, int> name_to_idx;
  for (size_t i = 0; i < cpp_names.size(); ++i)
    name_to_idx[cpp_names[i]] = (int)i;

  const int n_opt = af_only ? 0 : (int)variants_list.size();
  std::vector<FluorPrecomp> precomp(n_opt);
  std::vector<int> active_indices;

  if (!af_only) {
    CharacterVector vlist_names = variants_list.names();
    for (int i = 0; i < n_opt; ++i) {
      std::string fname = as<std::string>(vlist_names[i]);
      if (!name_to_idx.count(fname)) continue;

      FluorPrecomp& pc = precomp[i];
      pc.active     = true;
      pc.master_idx = name_to_idx[fname];
      pc.v_mats     = as<mat>(variants_list[i]);
      pc.n_variants = pc.v_mats.n_rows;
      if (pc.n_variants == 0) { pc.active = false; continue; }

      const rowvec master_row = spectra.row(pc.master_idx);

      pc.delta.set_size(pc.n_variants, D);
      for (uword v = 0; v < pc.n_variants; ++v)
        pc.delta.row(v) = pc.v_mats.row(v) - master_row;

      // Other-fluorophore pseudoinverse
      uvec keep(F - 1); uword ri = 0;
      for (uword r = 0; r < F; ++r) {
        if ((int)r != pc.master_idx) keep[ri++] = r;
      }

      const mat S_nof = spectra.rows(keep);             // (F-1) x D
      const mat U_nof = solve(S_nof * S_nof.t(), S_nof); // (F-1) x D

      pc.v_lib = U_nof * pc.delta.t();   // (F-1) x n_variants

      // Covariance-propagated leakage weights.
      // add a ridge (1e-4 * I) to stabilise small delta samples.
      {
        mat delta_obs   = as<mat>(delta_list[i]);        // n_obs x D
        mat delta_cov   = cov(delta_obs);
        delta_cov.diag() += 1e-4;                        // ridge regularisation
        const mat leakage_cov = mat(U_nof * delta_cov * U_nof.t()); // force eval
        pc.w_leakage = sqrt(abs(leakage_cov.diag())) + 1e-8;
      }

      // Rank-1 residual update helpers
      // r_lib(:,v) = delta(v)^T - S^T * (P * delta(v)^T)
      pc.r_lib  = pc.delta.t() - spectra.t() * (P * pc.delta.t());  // D x n_variants
      pc.r_dots.set_size(pc.n_variants);
      for (uword v = 0; v < pc.n_variants; ++v)
        pc.r_dots[v] = dot(pc.r_lib.col(v), pc.r_lib.col(v));

      active_indices.push_back(i);
    }
  } // end if (!af_only)

  // =========================================================================
  // SECTION 3 – Per-cell parallel loop
  // =========================================================================
  mat result(N, F + 2, fill::zeros);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel
{
  // Thread-local buffers – declared inside the parallel region so each
  // thread gets its own copy.  Using concrete mat/vec types (not auto or
  // expression templates) avoids lazy-evaluation pitfalls.
  vec   resid(D);
  vec   fluor_unmixed(F);
  mat   cell_S(F + 1, D);
  mat   cell_S_w(F + 1, D);   // sqrt_w-scaled copy used for WLS
  vec   unmixed_full(F + 1);
  vec   sqrt_w(D);             // per-cell detector weights (or ones)
  std::vector<int> best_v(n_opt, -1);

  struct Candidate { double score; int f_opt; uword v; };
  std::vector<Candidate> candidates;

#pragma omp for schedule(dynamic, 64)
  for (uword i = 0; i < N; ++i) {

    const vec cell_raw = raw_data.col(i);

    // =====================================================================
    // A. AF SELECTION
    // =====================================================================
    const vec init_f        = P * cell_raw;
    const vec base_resid_af = cell_raw - spectra.t() * init_f;
    const double base_resid_norm = std::max(norm(base_resid_af), 1e-8);
    const double base_fluor_l1   = std::max(dot(w_af, abs(init_f)),  1e-8);

    double best_score_af = datum::inf;
    uword  best_j_af     = 0;

    for (uword j = 0; j < nAF; ++j) {
      const double k_j = std::max(
        0.0, dot(cell_raw, r_lib_af.col(j)) / r_dots_af[j]);

      const vec r_j = base_resid_af - k_j * r_lib_af.col(j);
      // Multiplicative joint score: p_resid * p_fluor
      // Rewards large improvements on either axis without a mixing parameter.
      const double p_resid = norm(r_j) / base_resid_norm;
      const double p_fluor =
        dot(w_af, abs(init_f - k_j * v_lib_af.col(j))) / base_fluor_l1;

      const double score = p_resid * p_fluor;
      if (score < best_score_af) {
        best_score_af = score;
        best_j_af     = j;
      }
    }

    // =====================================================================
    // B. INITIAL JOINT SOLVE
    // =====================================================================
    cell_S.rows(0, F - 1) = spectra;
    cell_S.row(F)          = af_spectra.row(best_j_af);

    // Per-cell detector weights derived from the unweighted reconstruction.
    // Computed once here and held fixed for all subsequent solves this cell.
    if (cell_weight) {
      const vec coeff_init = solve(cell_S.t(), cell_raw, solve_opts::fast);  // F+1
      const vec y_hat      = cell_S.t() * coeff_init;                        // D
      for (uword d = 0; d < D; ++d)
        sqrt_w[d] = 1.0 / std::sqrt(std::max(std::abs(y_hat[d]), noise_floor));
    } else {
      sqrt_w.ones();
    }

    // WLS via sqrt_w row-scaling: solve (cell_S % sqrt_w.t())^T x = cell_raw % sqrt_w
    cell_S_w  = cell_S.each_row() % sqrt_w.t();
    unmixed_full  = solve(cell_S_w.t(), cell_raw % sqrt_w, solve_opts::fast);
    fluor_unmixed = unmixed_full.head(F);
    resid         = (cell_raw - cell_S.t() * unmixed_full) % sqrt_w;  // weighted residual

    // =====================================================================
    // B2. AF-ONLY EARLY RETURN
    // If no fluorophore variants were supplied, we are done: write output
    // and move to the next cell without entering the multi-pass loop.
    // =====================================================================
    if (af_only) {
      result(i, span(0, F - 1)) = fluor_unmixed.t();
      result(i, F)               = unmixed_full[F];
      result(i, F + 1)           = (double)best_j_af + 1.0;
      continue;
    }

    std::fill(best_v.begin(), best_v.end(), -1);

    // =====================================================================
    // C. JOINT VARIANT SELECTION
    // =====================================================================
    const int n_active = (int)active_indices.size();

    for (int pass = 0; pass < n_passes; ++pass) {

      const double rss_curr      = std::max(dot(resid, resid), 1e-12);
      const double rss_curr_sqrt = std::sqrt(rss_curr);

      candidates.clear();

      for (int ai = 0; ai < n_active; ++ai) {
        const int opt_i      = active_indices[ai];
        const FluorPrecomp& pc = precomp[opt_i];
        const double abund   = fluor_unmixed[pc.master_idx];
        if (abund < pos_thresholds[pc.master_idx]) continue;

        // Other-fluor unmixed values for leakage scoring
        vec other_unmixed(F - 1); uword oi = 0;
        for (uword r = 0; r < F; ++r) {
          if ((int)r != pc.master_idx) other_unmixed[oi++] = fluor_unmixed[r];
        }
        const double base_leakage =
          std::max(dot(pc.w_leakage, abs(other_unmixed)), 1e-8);

        const int cur_v = best_v[opt_i];   // -1 = master

        for (uword v = 0; v < pc.n_variants; ++v) {
          vec delta_r(D), delta_v_leakage(F - 1);

          if (cur_v < 0) {
            delta_r         = vec(pc.r_lib.col(v));
            delta_v_leakage = vec(pc.v_lib.col(v));
          } else {
            delta_r         = vec(pc.r_lib.col(v) - pc.r_lib.col(cur_v));
            delta_v_leakage = vec(pc.v_lib.col(v) - pc.v_lib.col(cur_v));
          }

          // Exact rank-1 RSS update
          const double cross   = dot(resid, delta_r);
          const double dr_sq   = dot(delta_r, delta_r);
          const double new_rss = rss_curr
          - 2.0 * abund * cross
          + abund * abund * dr_sq;
          const double resid_ratio =
          std::sqrt(std::max(new_rss, 0.0)) / rss_curr_sqrt;

          const vec    new_other    = other_unmixed - abund * delta_v_leakage;
          const double leakage_ratio =
            dot(pc.w_leakage, abs(new_other)) / base_leakage;

          const double composite = resid_ratio * leakage_ratio;
          if (composite < 1.0)
            candidates.push_back({composite, ai, v});
        }
      }

      if (candidates.empty()) break;

      // Sort best-first
      std::sort(candidates.begin(), candidates.end(),
                [](const Candidate& a, const Candidate& b){
                  return a.score < b.score;
                });

      // Conflict resolution: commit non-overlapping swaps
      std::vector<bool> committed(n_active, false);
      std::vector<vec>  committed_deltas;
      committed_deltas.reserve(n_active);

      std::vector<std::pair<int,uword>> commits;
      commits.reserve(n_active);

      for (const auto& cand : candidates) {
        if (committed[cand.f_opt]) continue;

        const int opt_i      = active_indices[cand.f_opt];
        const FluorPrecomp& pc = precomp[opt_i];
        const double abund   = fluor_unmixed[pc.master_idx];
        const int cur_v      = best_v[opt_i];

        // Concrete vec for the signed residual perturbation
        vec dr(D);
        if (cur_v < 0) {
          dr = vec(pc.r_lib.col(cand.v)) * abund;
        } else {
          dr = vec(pc.r_lib.col(cand.v) - pc.r_lib.col(cur_v)) * abund;
        }
        const double dr_norm = std::max(norm(dr), 1e-12);

        bool conflict = false;
        for (const vec& cd : committed_deltas) {
          const double cosine = std::abs(dot(dr, cd)) /
            (dr_norm * std::max(norm(cd), 1e-12));
          if (cosine > 0.5) { conflict = true; break; }
        }
        if (conflict) continue;

        committed[cand.f_opt] = true;
        committed_deltas.push_back(dr);
        commits.push_back({opt_i, cand.v});
      }

      if (commits.empty()) break;

      // Apply committed swaps and do one solve per pass
      for (auto& [opt_i, v] : commits) {
        best_v[opt_i] = (int)v;
        cell_S.row(precomp[opt_i].master_idx) = precomp[opt_i].v_mats.row(v);
      }

      cell_S_w      = cell_S.each_row() % sqrt_w.t();
      unmixed_full  = solve(cell_S_w.t(), cell_raw % sqrt_w, solve_opts::fast);
      fluor_unmixed = unmixed_full.head(F);
      resid         = (cell_raw - cell_S.t() * unmixed_full) % sqrt_w;

      if (norm(resid) < 1e-8 * norm(cell_raw)) break;
    }

    // =====================================================================
    // D. Write output
    // =====================================================================
    result(i, span(0, F - 1)) = fluor_unmixed.t();
    result(i, F)              = unmixed_full[F];
    result(i, F + 1)          = (double)best_j_af + 1.0;
  }
}
return result;
}
