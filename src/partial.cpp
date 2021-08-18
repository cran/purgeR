#include <Rcpp.h>
#include "constants.h"
#include "inbreeding_utils.h"
#include "genedrop.h"
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
//' @param dam Vector of dam ids.
//' @param sire Vector of sire ids.
//' @param anc_idx Index of ancestors.
//' @param Fi Vector of inbreeding coefficients.
//' @param mapa Map of ancestors
//' @param genedrop Enable genedrop simulation
//' @template seed-arg
//' @return A matrix of partial inbreeding coefficients. Fi(j) values can thus be read from row i and column j.
// [[Rcpp::export]]
Rcpp::NumericVector Fij_core_i_cpp(Rcpp::IntegerVector dam,
                                   Rcpp::IntegerVector sire,
                                   const int& anc_idx,
                                   Rcpp::LogicalVector mapa,
                                   Rcpp::NumericVector Fi,
                                   int genedrop = 0,
                                   Rcpp::Nullable<int> seed = R_NilValue) {
  
  // Input
  int N = dam.length();

  // Read partial inbreeding coefficients
  Rcpp::NumericVector pi(N);

  if (!genedrop) {
    
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
  } else {

    // Set seed (optional(
    if (seed.isNotNull()) {
      int s = Rcpp::as<int>(seed);
      Rcpp::Environment base_env("package:base");
      Rcpp::Function set_seed_r = base_env["set.seed"];
      set_seed_r(std::floor(std::fabs(s)));
    }
   
    // Set vector of dam and sires
    Rcpp::IntegerVector dam_clone = Rcpp::clone(dam);
    Rcpp::IntegerVector sire_clone = Rcpp::clone(sire);
    dam_clone[anc_idx] = 0;
    sire_clone[anc_idx] = 0;
    for (int iter (0); iter < genedrop; ++iter) {

      // each iteration vectors save alleles inherited from dam and sire
      Rcpp::IntegerVector dam_alleles;
      Rcpp::IntegerVector sire_alleles;
      int n_alleles (0);
      for (int i(0); i<N; ++i) {
        if (i == anc_idx) {
          ++n_alleles; // allele ID will start from 1 (not 0) in the reference ancestor
        }
        if (i < anc_idx) {
          dam_alleles.push_back(0);
          sire_alleles.push_back(0);
        } else {
          dam_alleles.push_back(genedrop_sim_allele(dam_clone[i], dam_alleles, sire_alleles, n_alleles));
          sire_alleles.push_back(genedrop_sim_allele(sire_clone[i], dam_alleles, sire_alleles, n_alleles));
          // with the design above, the reference ancestor is always heterozygote for alleles '1' and '2'
          if ((dam_alleles[i] == 1 && sire_alleles[i] == 1) || (dam_alleles[i] == 2 && sire_alleles[i] == 2)) pi[i] += 1.0;
        }
      }
      if (iter % 1000 == 0) Rcpp::checkUserInterrupt(); // check cancelation from user
    }
    pi = pi / (double)genedrop;
  }

  return pi;
}
