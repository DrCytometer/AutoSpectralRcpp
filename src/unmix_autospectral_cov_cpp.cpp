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
    int                    k_opt     = 1,
    int                    n_threads = 1
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
#pragma omp parallel for schedule(dynamic, 64)
  for (uword i = 0; i < n_cells; ++i) {

    // -- Thread-local scratch buffers --
    // Declared static thread_local so they are allocated once per thread and
    // reused across cells (resize/set_size is cheap; no heap traffic per cell).
    alignas(64) static thread_local vec    tl_cell_raw_vec;
    alignas(64) static thread_local rowvec tl_cell_raw;
    alignas(64) static thread_local vec    tl_init_f;
    alignas(64) static thread_local rowvec tl_base_resid_vec;
    alignas(64) static thread_local vec    tl_adj_f;
    alignas(64) static thread_local vec    tl_adj_resid;
    alignas(64) static thread_local rowvec tl_af_row;
    alignas(64) static thread_local mat    tl_cell_spectra_final;
    alignas(64) static thread_local rowvec tl_unmixed_full;
    alignas(64) static thread_local rowvec tl_unmixed_curr;
    alignas(64) static thread_local rowvec tl_resid;
    alignas(64) static thread_local vec    tl_cell_corrected_vec;
    alignas(64) static thread_local rowvec tl_curr_spectrum;
    alignas(64) static thread_local rowvec tl_d_j_raw;
    alignas(64) static thread_local rowvec tl_unmixed_other;
    alignas(64) static thread_local vec    tl_adj_other;
    alignas(64) static thread_local vec    tl_scores;
    alignas(64) static thread_local rowvec tl_t_full;
    alignas(64) static thread_local rowvec tl_t_resid;
    alignas(64) static thread_local rowvec tl_t_other;
    alignas(64) static thread_local mat    tl_trial_spectra;
    alignas(64) static thread_local rowvec tl_final_full;

    // thread_local index / order containers
    static thread_local std::vector<std::pair<double, size_t>> tl_fluor_order;
    static thread_local std::vector<uword>                     tl_top_k;

    tl_fluor_order.clear();

    // -- Load cell --
    tl_cell_raw_vec = raw_data.col(i);            // D
    tl_cell_raw     = tl_cell_raw_vec.t();        // 1 x D

    // -- A. AF selection via covariance-propagated composite score --
    tl_init_f = P * tl_cell_raw_vec;              // F

    // Baseline errors for normalisation
    tl_base_resid_vec.set_size(tl_cell_raw.n_elem);
    tl_base_resid_vec = tl_cell_raw - (tl_init_f.t() * spectra);  // 1 x D

    double base_e_fluor_af = dot(w_af, abs(tl_init_f)) + 1e-6;
    double base_e_resid_af = std::sqrt(dot(tl_base_resid_vec, tl_base_resid_vec)) + 1e-6;

    double best_score_af = datum::inf;
    uword  best_idx_af   = 0;

    tl_adj_f.set_size(n_fluors);
    tl_adj_resid.set_size(tl_cell_raw_vec.n_elem);

    for (uword j = 0; j < n_af; ++j) {
      double k_j = dot(tl_cell_raw_vec, r_lib_af.col(j)) / r_dots_af[j];
      if (k_j < 0.0) k_j = 0.0;

      tl_adj_f    = tl_init_f - k_j * v_lib_af.col(j);
      tl_adj_resid = tl_base_resid_vec.t() - k_j * r_lib_af.col(j);

      double e_fluor_j = dot(w_af, abs(tl_adj_f));
      double e_resid_j = std::sqrt(dot(tl_adj_resid, tl_adj_resid));

      double score_j = (e_resid_j / base_e_resid_af) * (e_fluor_j / base_e_fluor_af);

      if (score_j < best_score_af) {
        best_score_af = score_j;
        best_idx_af   = j;
      }
    }

    // -- B. Fluorophore variant optimisation --

    // Build AF-augmented spectra matrix [F+1 x D] in thread-local scratch.
    // AF row is last; its coefficient is the AF intensity.
    tl_af_row = af_spectra.row(best_idx_af);                      // 1 x D
    tl_cell_spectra_final.set_size(n_fluors + 1, spectra.n_cols);
    tl_cell_spectra_final.rows(0, n_fluors - 1) = spectra;
    tl_cell_spectra_final.row(n_fluors)         = tl_af_row;

    // Initial solve against cell_raw, AF included
    tl_unmixed_full = solve(tl_cell_spectra_final.t(), tl_cell_raw_vec,
                            solve_opts::fast).t();                 // 1 x (F+1)
    tl_unmixed_curr = tl_unmixed_full.head(n_fluors);
    double af_intensity = tl_unmixed_full(n_fluors);
    tl_resid = tl_cell_raw - (tl_unmixed_full * tl_cell_spectra_final); // 1 x D

    if (!active_f.empty()) {

      // Build sorted fluorophore order (brightest first), skipping sub-threshold
      for (size_t f : active_f) {
        int midx = precomp[f].master_idx;
        if (pos_thresholds.n_elem > 0 &&
            tl_unmixed_curr[(uword)midx] < pos_thresholds[(uword)midx]) continue;
        tl_fluor_order.push_back({ tl_unmixed_curr[(uword)midx], f });
      }

      std::sort(tl_fluor_order.begin(), tl_fluor_order.end(),
                [](const std::pair<double,size_t>& a,
                   const std::pair<double,size_t>& b) {
                  return a.first > b.first;
                });

      double resid_norm = std::sqrt(dot(tl_resid, tl_resid));

      for (auto const& [abundance, f] : tl_fluor_order) {
        const FluorPrecomp& pc = precomp[f];
        const uword nv = pc.n_variants;
        const uword F1 = pc.w.n_elem;  // F-1 (excludes fluorophore f, not AF)

        if (resid_norm < 1e-12) continue;

        // Other-fluorophore unmixed values (excludes f, not the AF row)
        tl_unmixed_other.set_size(F1);
        {
          uword oi = 0;
          for (uword p = 0; p < n_fluors; ++p)
            if ((int)p != pc.master_idx) tl_unmixed_other[oi++] = tl_unmixed_curr[p];
        }

        double base_e_fluor     = dot(pc.w, abs(tl_unmixed_other.t())) + 1e-6;
        double base_e_resid     = resid_norm + 1e-6;
        double align_scale      = abundance / base_e_resid;
        double base_e_fluor_inv = 1.0 / base_e_fluor;

        // Pre-screen scoring
        // curr_spectrum is updated after each accepted fluorophore
        tl_curr_spectrum = tl_cell_spectra_final.row((uword)pc.master_idx);  // 1 x D

        tl_cell_corrected_vec = (tl_cell_raw - af_intensity * tl_af_row).t(); // D

        tl_scores.set_size(nv);
        tl_d_j_raw.set_size(spectra.n_cols);
        tl_adj_other.set_size(F1);

        for (uword j = 0; j < nv; ++j) {
          tl_d_j_raw  = pc.v_mats.row(j) - tl_curr_spectrum;
          double d_j_norm = std::sqrt(dot(tl_d_j_raw, tl_d_j_raw));

          double p_resid_j = (d_j_norm > 1e-12)
            ? 1.0 - dot(tl_d_j_raw / d_j_norm, tl_resid) * align_scale
          : 1.0;

          double k_j = dot(tl_cell_corrected_vec, pc.r_lib.col(j)) / pc.r_dots[j];
          if (k_j < 0.0) k_j = 0.0;

          tl_adj_other = tl_unmixed_other.t() - k_j * pc.v_lib.col(j);
          double p_fluor_j = dot(pc.w, abs(tl_adj_other)) * base_e_fluor_inv;

          tl_scores[j] = p_resid_j * p_fluor_j;
        }

        // Top-K (ascending = better score)
        int k_eff = std::min(k_opt, (int)nv);
        tl_top_k.resize(nv);
        std::iota(tl_top_k.begin(), tl_top_k.end(), 0);
        std::partial_sort(tl_top_k.begin(), tl_top_k.begin() + k_eff, tl_top_k.end(),
                          [&](uword a, uword b) { return tl_scores[a] < tl_scores[b]; });
        tl_top_k.resize(k_eff);

        double best_score = 1.0;  // baseline: both ratios = 1

        tl_t_other.set_size(F1);

        for (int vi = 0; vi < k_eff; ++vi) {
          uword var_idx = tl_top_k[vi];

          // Trial spectra: substitute variant row, AF row unchanged
          tl_trial_spectra = tl_cell_spectra_final;
          tl_trial_spectra.row((uword)pc.master_idx) = pc.v_mats.row(var_idx);

          // Solve against cell_raw — AF re-estimated jointly
          tl_t_full  = solve(tl_trial_spectra.t(), tl_cell_raw_vec,
                             solve_opts::fast).t();                 // 1 x (F+1)
          tl_t_resid = tl_cell_raw - (tl_t_full * tl_trial_spectra); // 1 x D

          // Other-fluorophore slice (excludes f, not AF)
          {
            uword oi = 0;
            for (uword p = 0; p < n_fluors; ++p)
              if ((int)p != pc.master_idx) tl_t_other[oi++] = tl_t_full[p];
          }

          double e_fluor_t = dot(pc.w, abs(tl_t_other.t())) * base_e_fluor_inv;
          double e_resid_t = std::sqrt(dot(tl_t_resid, tl_t_resid));

          double composite = (e_resid_t / base_e_resid) * e_fluor_t;

          if (composite < best_score) {
            best_score              = composite;
            tl_unmixed_curr         = tl_t_full.head(n_fluors);
            af_intensity            = tl_t_full(n_fluors);
            tl_resid                = tl_t_resid;
            resid_norm              = e_resid_t;
            tl_cell_spectra_final   = tl_trial_spectra;
          }
        }
      }

      // Final solve with all accepted variants, AF included
      tl_final_full   = solve(tl_cell_spectra_final.t(), tl_cell_raw_vec,
                              solve_opts::fast).t();
      tl_unmixed_curr = tl_final_full.head(n_fluors);
      af_intensity    = tl_final_full(n_fluors);
    }

    res(i, span(0, n_fluors - 1)) = tl_unmixed_curr;
    res(i, n_fluors)     = af_intensity;
    res(i, n_fluors + 1) = (double)best_idx_af + 1.0;
  }

  return res;
}
