#include <RcppArmadillo.h>
#include <string>
#include <vector>
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
// Integrated AutoSpectral pipeline using covariance-propagated error scoring
// for both AF variant selection and fluorophore variant selection.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat unmix_autospectral_cov_cpp(
    arma::mat              raw_data_in,
    const arma::mat&       spectra,
    const arma::mat&       af_spectra,
    const CharacterVector& fluor_names,
    const CharacterVector& optimize_fluors,
    const arma::vec&       pos_thresholds,
    const List&            variants,
    const List&            delta_list,
    const List&            delta_norms,
    int                    k_opt       = 1,
    int                    n_threads   = 1,
    bool                   cell_weight = false,
    double                 noise_floor = 1.0
) {
  mat raw_data = raw_data_in.t();          // D x N
  const uword n_cells  = raw_data.n_cols;
  const uword n_fluors = spectra.n_rows;
  const uword n_af     = af_spectra.n_rows;

  // -------------------------------------------------------------------------
  // 1.  AF pre-computations
  // -------------------------------------------------------------------------
  mat P    = solve(spectra * spectra.t(), spectra, solve_opts::fast);  // F x D
  mat S_t  = spectra.t();                                               // D x F
  mat AF_t = af_spectra.t();                                            // D x n_af

  mat v_lib_af = P * AF_t;                           // F x n_af
  mat r_lib_af = AF_t - (S_t * v_lib_af);            // D x n_af

  vec r_dots_af(n_af);
  for (uword j = 0; j < n_af; ++j) {
    double d = dot(r_lib_af.col(j), r_lib_af.col(j));
    r_dots_af[j] = (d < 1e-10) ? 1e-10 : d;
  }

  // Covariance-propagated AF error weights
  mat af_cov       = cov(af_spectra);                // D x D
  mat fluor_af_cov = P * af_cov * P.t();            // F x F
  vec w_af         = sqrt(abs(fluor_af_cov.diag())) + 1e-8;  // F

  // -------------------------------------------------------------------------
  // 2.  Fluorophore variant pre-computations
  // -------------------------------------------------------------------------
  std::vector<std::string> cpp_names = as<std::vector<std::string>>(fluor_names);
  std::vector<std::string> cpp_opt   = as<std::vector<std::string>>(optimize_fluors);
  std::vector<std::string> v_names   = as<std::vector<std::string>>(variants.names());
  const size_t n_var_lists = v_names.size();

  std::map<std::string, int> name_to_idx;
  for (size_t i = 0; i < cpp_names.size(); ++i)
    name_to_idx[cpp_names[i]] = (int)i;

  struct FluorPrecomp {
    bool  active     = false;
    int   master_idx = -1;
    mat   v_mats;          // n_variants x D
    mat   v_lib;           // (F-1) x n_variants
    mat   r_lib;           // D x n_variants
    vec   r_dots;          // n_variants
    vec   w;               // F-1
    uword n_variants = 0;
  };

  std::vector<FluorPrecomp> precomp(n_var_lists);
  std::vector<size_t> active_f;

  for (size_t f = 0; f < n_var_lists; ++f) {
    FluorPrecomp& pc = precomp[f];

    bool requested = false;
    for (const auto& s : cpp_opt) {
      if (s == v_names[f]) { requested = true; break; }
    }

    auto it = name_to_idx.find(v_names[f]);
    if (!requested || it == name_to_idx.end()) continue;
    pc.master_idx = it->second;

    mat d_mat = as<mat>(delta_list[f]);
    vec d_nrm = as<vec>(delta_norms[f]);

    if (d_mat.is_empty() || d_nrm.is_empty() || max(d_nrm) < 1e-12) continue;

    pc.v_mats     = as<mat>(variants[f]);
    pc.n_variants = pc.v_mats.n_rows;

    uvec keep(n_fluors - 1);
    uword ri = 0;
    for (uword r = 0; r < n_fluors; ++r) {
      if ((int)r != pc.master_idx) keep[ri++] = r;
    }

    mat S_nof = spectra.rows(keep);                               // (F-1) x D
    mat U_nof = solve(S_nof * S_nof.t(), S_nof, solve_opts::fast); // (F-1) x D

    mat delta_cov = cov(d_mat);
    mat fluor_cov = U_nof * delta_cov * U_nof.t();
    pc.w = sqrt(abs(fluor_cov.diag())) + 1e-8;

    pc.v_lib  = U_nof * pc.v_mats.t();
    pc.r_lib  = pc.v_mats.t() - (S_nof.t() * pc.v_lib);

    pc.r_dots.set_size(pc.n_variants);
    for (uword j = 0; j < pc.n_variants; ++j) {
      double d = dot(pc.r_lib.col(j), pc.r_lib.col(j));
      pc.r_dots[j] = (d < 1e-10) ? 1e-10 : d;
    }

    pc.active = true;
    active_f.push_back(f);
  }

  // -------------------------------------------------------------------------
  // 3.  Output
  // -------------------------------------------------------------------------
  mat res(n_cells, n_fluors + 2);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  // -------------------------------------------------------------------------
  // 4.  Per-cell pipeline
  // -------------------------------------------------------------------------
#pragma omp parallel
{
  // Thread-local buffers — each thread gets its own copy
  const uword D = spectra.n_cols;

  vec    cell_raw_vec(D);
  rowvec cell_raw(D);
  vec    init_f(n_fluors);
  rowvec base_resid_vec(D);
  vec    adj_f(n_fluors);
  vec    adj_resid(D);
  rowvec af_row(D);
  mat    cell_spectra_final(n_fluors + 1, D);
  mat    cell_spectra_w(n_fluors + 1, D);     // sqrt_w-scaled copy for WLS
  rowvec cell_raw_w(D);                        // sqrt_w-scaled raw data
  vec    sqrt_w(D);                            // per-cell detector weights
  rowvec unmixed_full(n_fluors + 1);
  rowvec unmixed_curr(n_fluors);
  rowvec resid(D);
  vec    cell_corrected_vec(D);
  rowvec curr_spectrum(D);
  rowvec d_j_raw(D);
  rowvec unmixed_other;
  vec    adj_other;
  vec    scores;
  rowvec t_full(n_fluors + 1);
  rowvec t_resid(D);
  rowvec t_other;
  mat    trial_spectra(n_fluors + 1, D);
  mat    trial_spectra_w(n_fluors + 1, D);
  rowvec final_full(n_fluors + 1);

  std::vector<std::pair<double, size_t>> fluor_order;
  std::vector<uword>                     top_k;

#pragma omp for schedule(dynamic, 64)
  for (uword i = 0; i < n_cells; ++i) {

    fluor_order.clear();

    // -- Load cell --
    cell_raw_vec = raw_data.col(i);            // D
    cell_raw     = cell_raw_vec.t();           // 1 x D

    // -- A. AF selection (unweighted — weights not yet available) --
    init_f = P * cell_raw_vec;                 // F

    base_resid_vec = cell_raw - (init_f.t() * spectra);  // 1 x D

    double base_e_fluor_af = dot(w_af, abs(init_f)) + 1e-6;
    double base_e_resid_af = std::sqrt(dot(base_resid_vec, base_resid_vec)) + 1e-6;

    double best_score_af = datum::inf;
    uword  best_idx_af   = 0;

    adj_f.set_size(n_fluors);
    adj_resid.set_size(D);

    for (uword j = 0; j < n_af; ++j) {
      double k_j = dot(cell_raw_vec, r_lib_af.col(j)) / r_dots_af[j];
      if (k_j < 0.0) k_j = 0.0;

      adj_f     = init_f - k_j * v_lib_af.col(j);
      adj_resid = base_resid_vec.t() - k_j * r_lib_af.col(j);

      double e_fluor_j = dot(w_af, abs(adj_f));
      double e_resid_j = std::sqrt(dot(adj_resid, adj_resid));

      double score_j = (e_resid_j / base_e_resid_af) * (e_fluor_j / base_e_fluor_af);

      if (score_j < best_score_af) {
        best_score_af = score_j;
        best_idx_af   = j;
      }
    }

    // -- B. Build AF-augmented spectra matrix [F+1 x D] --
    af_row = af_spectra.row(best_idx_af);
    cell_spectra_final.rows(0, n_fluors - 1) = spectra;
    cell_spectra_final.row(n_fluors)         = af_row;

    // Per-cell detector weights derived from the unweighted reconstruction.
    // Ŷ = cell_spectra_final^T · (P_aug · cell_raw_vec), approximated cheaply
    // as cell_spectra_final^T · init_f_aug where the AF coefficient comes from
    // a rank-1 projection.  We use the full initial OLS solve for correctness.
    if (cell_weight) {
      const vec coeff_init = solve(cell_spectra_final.t(), cell_raw_vec,
                                   solve_opts::fast);               // F+1
      const vec y_hat      = cell_spectra_final.t() * coeff_init;  // D
      for (uword d = 0; d < D; ++d)
        sqrt_w[d] = 1.0 / std::sqrt(std::max(std::abs(y_hat[d]), noise_floor));
    } else {
      sqrt_w.ones();
    }

    cell_spectra_w = cell_spectra_final.each_row() % sqrt_w.t();
    cell_raw_w     = cell_raw % sqrt_w.t();

    // Initial WLS solve, AF included
    unmixed_full = solve(cell_spectra_w.t(), cell_raw_w.t(),
                         solve_opts::fast).t();                    // 1 x (F+1)
    unmixed_curr = unmixed_full.head(n_fluors);
    double af_intensity = unmixed_full(n_fluors);
    // Weighted residual for scoring
    resid = (cell_raw - (unmixed_full * cell_spectra_final)) % sqrt_w.t(); // 1 x D

    if (!active_f.empty()) {

      // Build sorted fluorophore order (brightest first), skipping sub-threshold
      for (size_t f : active_f) {
        int midx = precomp[f].master_idx;
        if (pos_thresholds.n_elem > 0 &&
            (uword)midx < pos_thresholds.n_elem &&
            unmixed_curr[(uword)midx] < pos_thresholds[(uword)midx]) continue;
        fluor_order.push_back({ unmixed_curr[(uword)midx], f });
      }

      std::sort(fluor_order.begin(), fluor_order.end(),
                [](const std::pair<double,size_t>& a,
                   const std::pair<double,size_t>& b) {
                  return a.first > b.first;
                });

      double resid_norm = std::sqrt(dot(resid, resid));

      for (auto const& [abundance, f] : fluor_order) {
        const FluorPrecomp& pc = precomp[f];
        const uword nv = pc.n_variants;
        const uword F1 = pc.w.n_elem;  // F-1

        if (resid_norm < 1e-12) continue;

        unmixed_other.set_size(F1);
        {
          uword oi = 0;
          for (uword p = 0; p < n_fluors; ++p)
            if ((int)p != pc.master_idx) unmixed_other[oi++] = unmixed_curr[p];
        }

        double base_e_fluor     = dot(pc.w, abs(unmixed_other.t())) + 1e-6;
        double base_e_resid     = resid_norm + 1e-6;
        double align_scale      = abundance / base_e_resid;
        double base_e_fluor_inv = 1.0 / base_e_fluor;

        curr_spectrum = cell_spectra_final.row((uword)pc.master_idx);  // 1 x D

        // AF-corrected cell signal for rank-1 leakage projection
        cell_corrected_vec = (cell_raw - af_intensity * af_row).t();   // D

        scores.set_size(nv);
        d_j_raw.set_size(D);
        adj_other.set_size(F1);

        for (uword j = 0; j < nv; ++j) {
          d_j_raw = pc.v_mats.row(j) - curr_spectrum;
          double d_j_norm = std::sqrt(dot(d_j_raw, d_j_raw));

          // Pre-screen: alignment with weighted residual (consistent with WLS)
          double p_resid_j = (d_j_norm > 1e-12)
            ? 1.0 - dot(d_j_raw / d_j_norm, resid) * align_scale
          : 1.0;

          double k_j = dot(cell_corrected_vec, pc.r_lib.col(j)) / pc.r_dots[j];
          if (k_j < 0.0) k_j = 0.0;

          adj_other = unmixed_other.t() - k_j * pc.v_lib.col(j);
          double p_fluor_j = dot(pc.w, abs(adj_other)) * base_e_fluor_inv;

          scores[j] = p_resid_j * p_fluor_j;
        }

        // Top-K (ascending = better score)
        int k_eff = std::min(k_opt, (int)nv);
        top_k.resize(nv);
        std::iota(top_k.begin(), top_k.end(), 0);
        std::partial_sort(top_k.begin(), top_k.begin() + k_eff, top_k.end(),
                          [&](uword a, uword b) { return scores[a] < scores[b]; });
        top_k.resize(k_eff);

        double best_score = 1.0;

        t_other.set_size(F1);

        for (int vi = 0; vi < k_eff; ++vi) {
          uword var_idx = top_k[vi];

          // Trial spectra: substitute variant row, AF row unchanged
          trial_spectra = cell_spectra_final;
          trial_spectra.row((uword)pc.master_idx) = pc.v_mats.row(var_idx);

          // WLS solve for trial variant
          trial_spectra_w = trial_spectra.each_row() % sqrt_w.t();
          t_full  = solve(trial_spectra_w.t(), cell_raw_w.t(),
                          solve_opts::fast).t();                   // 1 x (F+1)
          // Weighted residual
          t_resid = (cell_raw - (t_full * trial_spectra)) % sqrt_w.t(); // 1 x D

          {
            uword oi = 0;
            for (uword p = 0; p < n_fluors; ++p)
              if ((int)p != pc.master_idx) t_other[oi++] = t_full[p];
          }

          double e_fluor_t = dot(pc.w, abs(t_other.t())) * base_e_fluor_inv;
          double e_resid_t = std::sqrt(dot(t_resid, t_resid));

          double composite = (e_resid_t / base_e_resid) * e_fluor_t;

          if (composite < best_score) {
            best_score            = composite;
            unmixed_curr          = t_full.head(n_fluors);
            af_intensity          = t_full(n_fluors);
            resid                 = t_resid;
            resid_norm            = e_resid_t;
            cell_spectra_final    = trial_spectra;
            cell_spectra_w        = trial_spectra_w;
          }
        }
      }

      // Final WLS solve with all accepted variants
      final_full   = solve(cell_spectra_w.t(), cell_raw_w.t(),
                           solve_opts::fast).t();
      unmixed_curr = final_full.head(n_fluors);
      af_intensity = final_full(n_fluors);
    }

    res(i, span(0, n_fluors - 1)) = unmixed_curr;
    res(i, n_fluors)              = af_intensity;
    res(i, n_fluors + 1)          = (double)best_idx_af + 1.0;
  }
} // end omp parallel

return res;
}
