#include <Rcpp.h>
#include "constants.h"
#include "utils.h"
#include "inbreeding_utils.h"
using namespace Rcpp;

//' Reproductive value
//'
//' Computes the reproductive value
//' 
//' @template ped-arg
//' @template reference-arg
//' @template nameto-arg
//' @template target-arg
//' @param enable_correction Correct reproductive values.
//' @return The input dataframe, plus an additional column with reproductive values for the reference and target populations assumed.
//' @references
//' Hunter DC et al. 2019. Pedigree-based estimation of reproductive value. Journal of Heredity 110 (4): 433-444
// [[Rcpp::export]]
Rcpp::DataFrame reproductive_value(Rcpp::DataFrame ped,
                                   Rcpp::LogicalVector reference,
                                   std::string name_to,
                                   Nullable<Rcpp::LogicalVector> target = R_NilValue,
                                   bool enable_correction = true) {
  
  // Input values
  Rcpp::IntegerVector id = ped[id_col];
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  Rcpp::LogicalVector ped_ref = reference;
  Rcpp::LogicalVector ped_tgt;
  if (target.isNotNull()) ped_tgt = target;
  int N = ped.rows();

  // Map ancestors
  std::vector<std::vector<bool>> vis_ancestor (N, std::vector<bool>(N, false));
  for (int i(0); i < N; ++i) {
    int anc = id[i];
    for (int j(i+1); j < N; ++j) {
      if (dam[j] && ((dam[j] == anc) || vis_ancestor[i][dam[j]-1])) {
        vis_ancestor[i][j] = true;
      }
      if (sire[j] && ((sire[j] == anc) || vis_ancestor[i][sire[j]-1])) {
        vis_ancestor[i][j] = true;
      }
    }
  }

  // Truncate pedigree
  Rcpp::IntegerVector dam_t = ifelse(ped_ref, 0, dam);
  Rcpp::IntegerVector sire_t = ifelse(ped_ref, 0, sire);
  std::vector<int> id_ref; // ids in the reference population
  for (int i (0); i < N; ++i) if (ped_ref[i]) id_ref.push_back(id[i]);
  auto it_younger = std::max_element(std::begin(id_ref), std::end(id_ref));
  int id_younger = *it_younger;
  for (int i (0); i < id_younger; ++i) {
    if (dam_t[i] || sire_t[i]) {
      for (const auto& j: id_ref) {
        if (vis_ancestor[i][j-1]) {
          dam_t[i] = 0;
          sire_t[i] = 0;
          break;
        }
      }
    }
  }
  if (target.isNotNull()) {
    std::vector<int> id_tgt; // ids in the target population
    for (int i (0); i < N; ++i) if (ped_tgt[i]) id_tgt.push_back(id[i]);
    // check that target individual is an ancestor of a reference one
    for (const auto& i: id_tgt) {
      for (const auto& j: id_ref) {
        if (vis_ancestor[i][j-1]) throw std::logic_error("Found reference individual who is a descendant of the target population.");
      }
    }
    // continue truncating the pedigree
    std::vector<int> tmp_tgt;
    for (int i (id_younger); i < N; ++i) if (ped_tgt[i]) tmp_tgt.push_back(id[i]);
    for (int i (id_younger); i < N; ++i) {
      if (!is_in(id[i], id_tgt)) {
        for (const auto& j: id_tgt) {
          if (!vis_ancestor[i][j-1]) {
            dam_t[i] = 0;
            sire_t[i] = 0;
            break;
          }
        }
      }
    }
  }
  Rcpp::DataFrame ped_t = Rcpp::DataFrame::create(Named("id") = id,
                                                  Named("dam") = dam_t,
                                                  Named("sire") = sire_t);
  
  // Compute coancestry and additive relationship matrix
  std::vector<std::vector<double>> f = complete_coancestry(ped_t);
  std::vector<std::vector<double>> A = f;
  for (int i(0); i < N; ++i) {
    for (int j(0); j < i+1; ++j) {
      A[i][j] *= 2.0;
    }
  }
  
  // Uncorrected reproductive value for focal individuals
  Rcpp::NumericVector rp (N);
  for (const auto& i: id_ref) {
    for (int j(i-1); j < N; ++j) {
      rp[i-1] += A[j][i-1];
    }
  }
  
  if (!enable_correction) {
    ped[name_to] = rp;
    return ped;
  } else {
    // Corrected reproductive value
    // Standardization is made using contributions from all individuals to
    // the set of descendants of Ft
    
    double tgc = sum(rp); // total genetic contributions
    std::vector<int> id_descedants; // of id_focal (at Ft+Dt)
    for (int i (id_younger); i < N; ++i) id_descedants.push_back(id[i]);
    for (int i (0); i < id_younger; ++i) {
      if (!is_in(id[i], id_ref)) {
        for (const auto& j: id_descedants) {
          if (vis_ancestor[i][j-1]) {
            for (int k(i); k < N; ++k) {
              tgc += A[k][i];
            }
            break;
          }
        }
      }
    }
    rp = rp / tgc;
    ped[name_to] = rp;
    return ped;
  }
}
