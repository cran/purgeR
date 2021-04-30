#include <Rcpp.h>
#include "inbreeding_utils.h"

//' Generation numbers for minimum coancestry
//'
//' Minimum coancestry matrix requires to split data in pseudo-genrations. This is to help to skip the computation of unnecesary
//' kinship coefficients. The general rule is that each individual belongs to only one generation number, and ids cannot belong to
//' the same generation number than its ancestors.
//' 
//' @name coancestry_t
//' @param dam Vector of dam identities.
//' @param sire Vector of sire identities.
//' @return Vector of generation number approximations
std::vector<int> coancestry_t(Rcpp::IntegerVector dam,
                              Rcpp::IntegerVector sire) {
  
  int N = dam.length();
  int gen(0);
  std::vector<int> t;
  for (int i(0); i<N; ++i) {
    if (dam[i] != 0 && t[dam[i]-1]==gen) ++gen;
    else if (sire[i] != 0 && t[sire[i]-1]==gen) ++gen;
    t.push_back(gen);
  }
  return t;
}

//' Search ancestors for minimum coancestry
//'
//' This function searches recursively the ancestors required to compute kinship coefficients under the
//' minimum coancestry algorithm
//' 
//' @name search_ancestors_min_coancestry
//' @param id Vector of id.
//' @param dam Vector of dam.
//' @param sire Vector of sire.
//' @param t Vector of generation numbers (see coancestry_t function).
//' @param log_pairs Set of coancestry coordinates to gather required individual pairs (returned as reference).
//' @param i Individual i coordinate.
//' @param j Individual j coordinate.
void search_ancestors_min_coancestry (const std::vector<int>& id,
                                      const std::vector<int>& dam,
                                      const std::vector<int>& sire,
                                      const std::vector<int>& t,
                                      std::set<std::pair<int, int>>& log_pairs,
                                      const int& i,
                                      const int& j) {
  
  std::pair<int, int> f = std::make_pair(i,j);
  std::pair<int, int> r = std::make_pair(j,i);
  if (( log_pairs.count(f) != 0) || ( log_pairs.count(r) != 0) || !i || !j) return;
  else if (i<j) log_pairs.insert(f);
  else log_pairs.insert(r);
  
  if (t[i-1] == t[j-1]) {
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, dam[i-1], dam[j-1]);
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, dam[i-1], sire[j-1]);
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, sire[i-1], dam[j-1]);
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, sire[i-1], sire[j-1]);
  } else if (t[i-1] < t[j-1]) {
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, id[i-1], dam[j-1]);
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, id[i-1], sire[j-1]);
  } else {
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, id[j-1], dam[i-1]);
    search_ancestors_min_coancestry (id, dam, sire, t, log_pairs, id[j-1], sire[i-1]);
  }
}

//' Minimum coancestry matrix
//'
//' Computes a minimum map with coancestry coefficients for all individuals required to compute inbreeding coefficients.
//' Unnecessary coefficients are skipped. 
//' 
//' @name minimum_coancestry
//' @template ped-arg
//' @return Coancestry matrix (as a dedicated, queryable object)
coancestry_matrix minimum_coancestry(Rcpp::DataFrame ped) {
  
  // Get input
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
    if (i % 100 == 0) Rcpp::checkUserInterrupt(); // check cancellation from user
  }
  int Npairs = log_pairs.size();
  
  coancestry_matrix f_map;
  f_map.coancestry.reserve (Npairs);
  int count (0);
  for (std::set<std::pair<int, int>>::iterator it = log_pairs.begin(); it != log_pairs.end(); ++it) {
    int i (std::get<0>(*it));
    int j (std::get<1>(*it));
    if (i == j) {
      f_map[std::make_pair(i,j)] = 0.5 * (1.0 + f_map[std::make_pair(dam[i-1], sire[i-1])]);
    } else if (t[i-1]==t[j-1]) {
      double fAC (f_map[std::make_pair(dam[i-1], dam[j-1])]);
      double fAD (f_map[std::make_pair(dam[i-1], sire[j-1])]);
      double fBC (f_map[std::make_pair(sire[i-1], dam[j-1])]);
      double fBD (f_map[std::make_pair(sire[i-1], sire[j-1])]);
      f_map[std::make_pair(i,j)] = 0.25 * (fAC + fAD + fBC + fBD);
    } else {
      double fAE(0.0), fAH(0.0);
      fAE = f_map[std::make_pair(id[i-1], dam[j-1])];
      fAH = f_map[std::make_pair(id[i-1], sire[j-1])];
      f_map[std::make_pair(i,j)] = 0.5 * (fAE + fAH );
    }
    if (++count % 100 == 0) Rcpp::checkUserInterrupt(); // check cancelation from user
  }
  return f_map;
}

//' Complete coancestry matrix
//'
//' Computes the complete, triangular matrix, with the pairwise coancestry coefficients of all
//' individuals in the pedigree. 
//' 
//' @name complete_coancestry
//' @template ped-arg
//' @return Coancestry matrix
std::vector<std::vector<double>> complete_coancestry(Rcpp::DataFrame ped) {
  
  // Get input
  Rcpp::IntegerVector id = ped[id_col];
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  int N = ped.rows();
  
  std::vector<std::vector<double>> f;
  for (int i(0); i<N; ++i) {
    std::vector<double> id_i_f;
    for (int j(0); j<(i+1); ++j) {
      if (i==j) id_i_f.push_back(0.5);
      else id_i_f.push_back(0.0);
    }
    f.push_back(id_i_f);
  }
  for (int i(0); i<N; ++i) {
    int id_i = id[i];
    for (int j(0); j<(i+1); ++j) {
      int id_j = id[j];
      if (id_i == id_j) {
        if (dam[id_i-1] == 0 || sire[id_i-1] == 0) continue;
        double F_i = 0.0;
        if (dam[id_i-1] >= sire[id_i-1]) F_i = f[dam[id_i-1]-1][sire[id_i-1]-1];
        else F_i = f[sire[id_i-1]-1][dam[id_i-1]-1];
        f[id_i-1][id_i-1] = 0.5 * (1.0 + F_i);
      } else {
        double f_ijdam = 0.0;
        double f_ijsire = 0.0;
        if (dam[id_i-1] > 0) {
          if (id_j >= dam[id_i-1]) f_ijdam = f[id_j-1][dam[id_i-1]-1];
          else f_ijdam = f[dam[id_i-1]-1][id_j-1];
        }
        if (sire[id_i-1] > 0) {
          if (id_j >= sire[id_i-1]) f_ijsire = f[id_j-1][sire[id_i-1]-1];
          else f_ijsire = f[sire[id_i-1]-1][id_j-1];
        }
        f[id_i-1][id_j-1] = 0.5 * (f_ijdam + f_ijsire);
      }
    }
    if (i % 100 == 0) Rcpp::checkUserInterrupt(); // check cancellation from user
  }
  return f;
}
