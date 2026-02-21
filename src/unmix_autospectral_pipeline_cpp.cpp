#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// Helper: Fast Top-K selection (replaces sort_index for small k)
std::vector<uword> find_top_k(const vec& scores, int k) {
  int n = scores.n_elem;
  int k_eff = std::min(k, n);
  if (k_eff <= 0) return {};
  
  std::vector<uword> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  
  // Partial sort only the top K elements
  std::partial_sort(indices.begin(), indices.begin() + k_eff, indices.end(),
                    [&scores](uword i, uword j) { return scores[i] > scores[j]; });
  
  indices.resize(k_eff);
  return indices;
}

// [[Rcpp::export]]
arma::mat unmix_autospectral_pipeline_cpp(
    arma::mat raw_data_in,
    const arma::mat& spectra,
    const arma::mat& af_spectra,
    const CharacterVector& fluor_names,
    const CharacterVector& optimize_fluors,
    const arma::vec& pos_thresholds,
    const List& variants,
    const List& delta_list,
    const List& delta_norms,
    int k_opt = 1,
    int n_threads = 1,
    bool optimize = true,
    bool use_dist0 = true
) {
  mat raw_data = raw_data_in.t();
  uword n_cells = raw_data.n_cols;
  uword n_fluors = spectra.n_rows;
  uword n_af = af_spectra.n_rows;
  
  // --- 1. PRE-CALCULATIONS ---
  mat P = solve(spectra * spectra.t(), spectra);
  mat S_t = spectra.t();
  mat AF_t = af_spectra.t();
  
  mat v_library_af = P * AF_t;
  mat r_library_af = AF_t - (S_t * v_library_af);
  vec r_dots_af(n_af);
  for(uword j = 0; j < n_af; ++j) {
    double d = dot(r_library_af.col(j), r_library_af.col(j));
    r_dots_af[j] = (d == 0) ? 1e-10 : d;
  }
  
  // --- 2. OPTIMIZATION SETUP ---
  std::vector<std::string> cpp_names = as<std::vector<std::string>>(fluor_names);
  std::map<std::string, int> name_to_idx;
  for(size_t i = 0; i < cpp_names.size(); ++i) name_to_idx[cpp_names[i]] = (int)i;
  
  size_t n_var_total = (optimize) ? (size_t)variants.size() : 0;
  std::vector<mat> D_scaled(n_var_total);
  std::vector<mat> v_mats(n_var_total);
  std::vector<int> var_to_master(n_var_total);
  std::vector<size_t> active_opt_indices;
  
  if (optimize) {
    std::vector<std::string> v_names = as<std::vector<std::string>>(variants.names());
    std::vector<std::string> cpp_opt_names = as<std::vector<std::string>>(optimize_fluors);
    
    for (size_t f = 0; f < n_var_total; ++f) {
      mat d_mat = as<mat>(delta_list[f]);
      vec d_norm = as<vec>(delta_norms[f]);
      
      // Pre-scale delta matrices by their norms
      for(uword r = 0; r < d_mat.n_rows; ++r) {
        double n = d_norm[r];
        d_mat.row(r) /= (n > 1e-12 ? n : 1.0);
      }
      D_scaled[f] = d_mat;
      v_mats[f] = as<mat>(variants[f]);
      var_to_master[f] = name_to_idx.count(v_names[f]) ? name_to_idx[v_names[f]] : -1;
      
      bool is_req = false;
      for (const auto& name : cpp_opt_names) if (name == v_names[f]) { is_req = true; break; }
      if (is_req && !d_norm.is_empty() && (max(d_norm) > 1e-12)) {
        active_opt_indices.push_back(f);
      }
    }
  }
  
  mat res(n_cells, n_fluors + 2);
  
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  
  // --- 3. PARALLEL PIPELINE ---
#pragma omp parallel for schedule(dynamic, 64)
  for(uword i = 0; i < n_cells; ++i) {
    // alignas(64) ensures data fits into CPU cache lines
    alignas(64) static thread_local rowvec cell_unmixed, cell_resid_raw, unmixed_curr, resid, t_unmix, t_resid;
    alignas(64) static thread_local mat spectra_curr;
    static thread_local std::vector<std::pair<double, size_t>> fluor_order;
    static thread_local std::vector<uword> m_to_c_row;
    
    fluor_order.clear();
    if(m_to_c_row.size() != n_var_total) m_to_c_row.assign(n_var_total, 0);
    
    vec cell_raw_vec = raw_data.col(i); 
    rowvec cell_raw = cell_raw_vec.t();
    
    // A. AF EXTRACTION
    vec init_f = (P * cell_raw_vec);
    double best_k_af = 0; uword best_idx_af = 0;
    
    if (use_dist0) {
      double min_err = datum::inf;
      for(uword j = 0; j < n_af; ++j) {
        double k = dot(cell_raw_vec, r_library_af.col(j)) / r_dots_af[j];
        double err = sum(abs(init_f - (k * v_library_af.col(j))));
        if(err < min_err) { min_err = err; best_k_af = k; best_idx_af = j; }
      }
    } else {
      rowvec init_res = cell_raw - (init_f.t() * spectra);
      double max_score = -datum::inf;
      for(uword j = 0; j < n_af; ++j) {
        double score = dot(init_res, AF_t.col(j));
        if(score > max_score) { max_score = score; best_idx_af = j; }
      }
      best_k_af = dot(init_res, r_library_af.col(best_idx_af).t()) / r_dots_af[best_idx_af];
    }
    
    cell_resid_raw = cell_raw - (best_k_af * af_spectra.row(best_idx_af));
    cell_unmixed = (P * cell_resid_raw.t()).t();
    
    // B. SPECTRAL OPTIMIZATION
    if (optimize && !active_opt_indices.empty()) {
      mat cell_spectra_final = spectra;
      uvec pos_idx = find(cell_unmixed.t() >= pos_thresholds);
      
      // only proceed if some fluorophores are positive
      if (pos_idx.n_elem > 0) {
        spectra_curr = cell_spectra_final.rows(pos_idx);
        unmixed_curr = solve(spectra_curr.t(), cell_resid_raw.t(), solve_opts::fast).t();
        resid = cell_resid_raw - (unmixed_curr * spectra_curr);
        double err_final = sum(abs(resid));
        double r_n = std::sqrt(dot(resid, resid));
        
        for (size_t f_idx : active_opt_indices) {
          int m_idx = var_to_master[f_idx];
          for(uword p = 0; p < pos_idx.n_elem; ++p) {
            if((int)pos_idx[p] == m_idx) { 
              fluor_order.push_back({unmixed_curr[p], f_idx}); 
              m_to_c_row[f_idx] = p; break; 
            }
          }
        }
        
        // brightest to dimmest
        if(!fluor_order.empty()){
          std::sort(fluor_order.begin(), fluor_order.end(), std::greater<>());
          double inv_rn = 1.0 / (r_n + 1e-12);
          
          for (auto const& pair : fluor_order) {
            size_t f_idx = pair.second; int m_idx = var_to_master[f_idx]; uword r_curr = m_to_c_row[f_idx];
            
            if (r_n > 1e-12) {
              // scoring: Vector-matrix multiplication with pre-scaled delta
              vec scores = (D_scaled[f_idx] * resid.t()) * (unmixed_curr[r_curr] * inv_rn);
              
              // Select top K variants: Only look at best variants, avoid full sort
              std::vector<uword> topK = find_top_k(scores, k_opt);
              
              for (uword v_idx : topK) {
                rowvec backup = spectra_curr.row(r_curr);
                spectra_curr.row(r_curr) = v_mats[f_idx].row(v_idx);
                
                t_unmix = solve(spectra_curr.t(), cell_resid_raw.t(), solve_opts::fast).t();
                t_resid = cell_resid_raw - (t_unmix * spectra_curr);
                
                // Evaluate based on L1 (absolute residuals)
                double current_err = sum(abs(t_resid));
                if (current_err < err_final) {
                  err_final = current_err; unmixed_curr = t_unmix; resid = t_resid;
                  r_n = std::sqrt(dot(resid, resid)); 
                  inv_rn = 1.0 / (r_n + 1e-12);
                  cell_spectra_final.row(m_idx) = spectra_curr.row(r_curr);
                } else { 
                  spectra_curr.row(r_curr) = backup; 
                }
              }
            }
          }
          cell_unmixed = solve(cell_spectra_final.t(), cell_resid_raw.t(), solve_opts::fast).t();
        }
      }
    }
    
    res(i, span(0, n_fluors-1)) = cell_unmixed;
    res(i, n_fluors) = best_k_af;
    res(i, n_fluors + 1) = (double)best_idx_af + 1;
  }
  return res;
}