#include <Rcpp.h>
#include "constants.h"
#include "inbreeding_utils.h"
#include "partial.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
//' Partial inbreeding coefficient (core function)
//'
//' Computes partial inbreeding coefficients, Fi(j).
//' A coefficient Fi(j) can be read as the probability of individual i being
//' homozygous for alleles derived from ancestor j
//' 
//' @name Fij_core_i_cpp
//' @template ped-arg
//' @param anc_idx Index of ancestors.
//' @param Fi Vector of inbreeding coefficients.
//' @param mapa Map of ancestors
//' @return A matrix of partial inbreeding coefficients. Fi(j) values can thus be read from row i and column j.
// [[Rcpp::export]]
Rcpp::NumericVector Fij_core_i_cpp(Rcpp::DataFrame ped,
                                  const int& anc_idx,
                                  Rcpp::LogicalVector mapa,
                                  Rcpp::NumericVector Fi) {
  
  // Input
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  int N = ped.rows();

  // Read partial inbreeding coefficients
  Rcpp::NumericVector pi(N);

  // founders partial kinship
  pcoancestry_matrix pkm;
  pkm[std::make_pair(anc_idx, anc_idx)] = 0.5;

  // other ancestors
  for (int j(anc_idx); j < N; ++j) {
    if (mapa[j]) {
      for (int k(0); k < j; ++k) {
        double f_mw = pkm[std::make_pair(dam[j] - 1, k)];
        double f_pw = pkm[std::make_pair(sire[j] - 1, k)];
        pkm[std::make_pair(j, k)] = 0.5 * (f_mw + f_pw);
      }
    }
    // self-coancestry (j==k)
    double f_xj = pkm[std::make_pair(j, anc_idx)];
    double f_mp = pkm[std::make_pair(dam[j] - 1, sire[j] - 1)];
    pkm[std::make_pair(j, j)] = f_xj + 0.5 * f_mp;
    pi[j] = pkm[std::make_pair(dam[j] - 1, sire[j] - 1)];
  }
  return pi;
}
