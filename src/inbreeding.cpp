#include <Rcpp.h>
#include "constants.h"
#include "utils.h"
#include "inbreeding_utils.h"
#include "genedrop.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]

//' Inbreeding coefficient
//'
//' Computes the standard inbreeding coefficient (\emph{F}).
//' This is the probability that two alleles on a locus are identical by descent (Falconer and Mackay 1996, Wright 1922), calculated from the genealogical coancestry matrix (Malécot 1948).
//' 
//' @template ped-arg
//' @template nameto-arg
//' @return The input dataframe, plus an additional column named "F" with individual inbreeding coefficient values.
//' @encoding UTF-8
//' @references
//' \itemize{
//'   \item{Falconer DS, Mackay TFC. 1996. Introduction to Quantitative Genetics. 4th edition. Longman, Essex, U.K.}
//'   \item{Malécot G, 1948. Les Mathématiques de l’hérédité. Masson & Cie., Paris.}
//'   \item{Wright S. 1922. Coefficients of inbreeding and relationship. The American Naturalist 56: 330-338.}
//'}
// [[Rcpp::export]]
DataFrame F(DataFrame ped, std::string name_to) {
  
  // Check errors
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  
  // Compute coancestry
  coancestry_matrix f_map = minimum_coancestry(ped);
  
  // Calculate inbreeding (F)
  int N = ped.rows();
  Rcpp::NumericVector F;
  for (int i(0); i<N; ++i) F.push_back(f_map[std::make_pair(dam[i], sire[i])]);
  ped[name_to] = F;
  return ped;
}

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
//' Ancestral inbreeding coefficient
//'
//' Computes the ancestral inbreeding coefficient (\emph{Fa}).
//' This is the probability that an allele has been in homozygosity in at least one ancestor (Ballou 1997).
//' A genedrop approach is included to compute unbiased estimates of \emph{Fa} (Baumung et al. 2015).
//' 
//' @template ped-arg
//' @template nameto-arg
//' @template Fi-arg
//' @param genedrop Number of genedrop iterations to run. If set to zero (as default), Ballou's Fa is computed.
//' @template seed-arg
//' @return The input dataframe, plus an additional column named "Fa" with individual ancestral inbreeding coefficient values.
//' @references
//' \itemize{
//'   \item{Ballou JD. 1997. Ancestral inbreeding only minimally affects inbreeding depression in mammalian populations. J Hered. 88:169–178.}
//'   \item{Baumung et al. 2015. GRAIN: A computer program to calculate ancestral and partial inbreeding coefficients using a gene dropping approach. Journal of Animal Breeding and Genetics 132: 100-108.}
//'}
// [[Rcpp::export]]
DataFrame Fa(DataFrame ped,
             Rcpp::NumericVector Fi,
             std::string name_to,
             int genedrop = 0,
             Nullable<int> seed = R_NilValue) {
  
  // Check errors
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  
  if (genedrop && seed.isNotNull()) {
    int s = Rcpp::as<int>(seed);
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(s)));
  }
  
  // Calculate ancestral inbreeding (Fa)
  int N = ped.rows();
  Rcpp::NumericVector Fa;
  
  if (!genedrop) {
    for (int i(0); i<N; ++i) {
      int dam_id = dam[i];
      int sire_id = sire[i];
      double f_dame(0.0);
      double fa_dame(0.0);
      double f_sir (0.0);
      double fa_sir (0.0);
      
      if (dam_id)  {
        f_dame = Fi[dam_id-1];
        fa_dame = Fa[dam_id-1];
      }
      if (sire_id) {
        f_sir  = Fi[sire_id -1];
        fa_sir  = Fa[sire_id -1];
      }
      Fa.push_back(0.5 * (fa_dame + (1.0 - fa_dame) * f_dame + fa_sir + (1.0 - fa_sir) * f_sir));
    }
  } else {
    bool display_progress = true;
    Progress p(genedrop, display_progress);
    for (int iter (0); iter < genedrop; ++iter) {
      
      // each iteration vectors save alleles inherited from dam and sire
      Rcpp::IntegerVector dam_alleles;
      Rcpp::IntegerVector sire_alleles;
      int n_alleles (0);
      for (int i(0); i<N; ++i) {
        dam_alleles.push_back(genedrop_sim_allele(dam[i], dam_alleles, sire_alleles, n_alleles));
        sire_alleles.push_back(genedrop_sim_allele(sire[i], dam_alleles, sire_alleles, n_alleles));
      }
      Rcpp::LogicalVector ibd_dam_alleles;
      Rcpp::LogicalVector ibd_sire_alleles;
      genedrop_ibd(dam, sire, dam_alleles, sire_alleles, ibd_dam_alleles, ibd_sire_alleles);
      if (iter) Fa += genedrop_read_Fa (dam, sire, dam_alleles, sire_alleles, ibd_dam_alleles, ibd_sire_alleles);
      else Fa = genedrop_read_Fa (dam, sire, dam_alleles, sire_alleles, ibd_dam_alleles, ibd_sire_alleles);
      p.increment(); // update progress
      if (iter % 10000 == 0) Rcpp::checkUserInterrupt(); // check cancelation from user
    }
    Fa = Fa / (double)genedrop;
  }
  ped[name_to] = Fa;
  return ped;
}

//' Purged inbreeding coefficient
//'
//' Computes the purged inbreeding coefficient (\emph{g}).
//' This is the probability that two alleles on a locus are identical by descent,
//' but relative to deleterious recessive alleles (García-Dorado 2012). The reduction
//' in \emph{g} relative to standard inbreeding (\emph{F}) is given by an effective purging
//' coefficient (\emph{d}), that measures the strength of the deleterious recessive
//' component in the genome. The coefficient \emph{g} is computed following the methods
//' for pedigrees in García-Dorado (2012) and García-Dorado et al. (2016).
//' 
//' @template ped-arg
//' @param d Purging coefficient (taking values between 0.0 and 0.5).
//' @template nameto-arg
//' @template Fi-arg
//' @return The input dataframe, plus an additional column named "g" followed by the purging coefficient, containing purged inbreeding coefficient values.
//' @encoding UTF-8
//' @references
//' \itemize{
//'   \item{García-Dorado. 2012. Understanding and predicting the fitness decline of shrunk populations: Inbreeding, purging, mutation, and standard selection. Genetics 190: 1-16.}
//'   \item{García-Dorado et al. 2016. Predictive model and software for inbreeding-purging analysis of pedigreed populations. G3 6: 3593-3601.}
//' }
// [[Rcpp::export]]
DataFrame g(Rcpp::DataFrame ped,
            double d,
            Rcpp::NumericVector Fi,
            std::string name_to) {

  Rcpp::IntegerVector id = ped[id_col];
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  
  // Compute pseudo-numbers of generation (to search ancestors)
  std::vector<int> t = coancestry_t(dam, sire);
  
  // Compute the minimum required pairs of ancestors
  int N = ped.rows();
  std::set<std::pair<int, int>> log_pairs;
  std::vector<int> id_std = Rcpp::as<std::vector<int>>(id);
  std::vector<int> dam_std = Rcpp::as<std::vector<int>>(dam);
  std::vector<int> sire_std = Rcpp::as<std::vector<int>>(sire);
  for (int i(N); i>=1; --i) {
    search_ancestors_min_coancestry (id_std, dam_std, sire_std, t, log_pairs, dam[i-1], sire[i-1]);
    if (i % 100 == 0) Rcpp::checkUserInterrupt(); // check cancelation from user
  }
  int Npairs = log_pairs.size();
  
  // Compute purged coancestry
  coancestry_matrix g_map;
  g_map.coancestry.reserve (Npairs);
  for (std::set<std::pair<int, int>>::iterator it=log_pairs.begin(); it!=log_pairs.end(); ++it) {
    int i (std::get<0>(*it));
    int j (std::get<1>(*it));
    if (i == j) {
      g_map[std::make_pair(i, j)] = 0.5 * (1.0 + g_map[std::make_pair(dam[i-1], sire[i-1])]) * (1.0 - 2.0*d*Fi[i-1]);
    } else if (t[i-1] == t[j-1]) {
      double fAC (g_map[std::make_pair(dam[i-1], dam[j-1])]);
      double fAD (g_map[std::make_pair(dam[i-1], sire[j-1])]);
      double fBC (g_map[std::make_pair(sire[i-1], dam[j-1])]);
      double fBD (g_map[std::make_pair(sire[i-1], sire[j-1])]);
      g_map[std::make_pair(i, j)] = 0.25 * (fAC + fAD + fBC + fBD) * (1.0 - d*(Fi[i-1]+Fi[j-1]));
    } else {
      double fAE(0.0), fAH(0.0);
      fAE = g_map[std::make_pair(id[i-1], dam[j-1])];
      fAH = g_map[std::make_pair(id[i-1], sire[j-1])];
      g_map[std::make_pair(i, j)] = 0.5 * (fAE + fAH) * (1.0 - d*Fi[j-1]);
    }
  }
  
  // Calculate purged inbreeding (g)
  Rcpp::NumericVector g;
  for (int i(0); i<N; ++i) g.push_back(g_map[std::make_pair(dam[i], sire[i])]);
  ped[name_to] = g;
  return ped;
}
