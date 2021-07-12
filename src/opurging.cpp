#include <Rcpp.h>
#include "constants.h"
#include "inbreeding_utils.h"
#include "partial.h"

void map_ij_distance (const int&, const int&, Rcpp::IntegerVector, Rcpp::IntegerVector, std::vector<int>&, const int&);
std::vector<std::set<int>> map_ancestor_ibd (Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::NumericVector);

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
//' Opportunity of purging
//'
//' The potential reduction in individual inbreeding load can be
//' estimated by means of the opportunity of purging (\emph{O}) and expressed
//' opportunity of purging (\emph{Oe}) parameters described by Gulisija
//' and Crow (2007). Whereas \emph{O} relates to the total potential reduction
//' of the inbreeding load in an individual, as a consequence of it having
//' inbred ancestors, \emph{Oe} relates to the expressed potential reduction of the
//' inbreeding load. In both cases, these measures are referred to fully recessive,
//' high effect size alleles (e.g. lethals). For complex pedigrees, involving more than one
//' autozygous individual per path from a reference individual to an
//' ancestor, these estimates are estimated following an heuristic approach
//' (see details below).
//' 
//' In simple pedigrees, the opportunity of purging (\emph{O}) and the expressed
//' opportunity of purging (\emph{Oe}) are estimated as in Gulisija and Crow (2007).
//' For complex pedigrees involving more than one autozygous individual per
//' path from an individual to an ancestor, \emph{O} and \emph{Oe} in the closer ancestors
//' need to be discounted for what was already accounted for in their predecessors.
//' To solve this problem, Gulisija and Crow (2007) provide expression to 
//' correct \emph{O} and \emph{Oe} (see equations 21 and 22 in the manuscript).
//' 
//' Here, an heuristic approach is used to prevent the inflation of \emph{O} and \emph{Oe},
//' and avoid the use of additional looped expressions that may result in an
//' excessive computational cost. To do so, when using \code{ip_op(complex = TRUE)}
//' only the contribution of the most recent ancestors in a path will be considered.
//' This may not provide exact values of \emph{O} and \emph{Oe}, but we expect little bias, since
//' more distant ancestors also contribute lesser to \emph{O} and \emph{Oe}.
//' 
//' @template ped-arg
//' @param pi Partial inbreeding matrix
//' @template Fi-arg
//' @param name_O A string naming the new output column for total opportunity of purging (defaults to "O") 
//' @param name_Oe A string naming the new output column for the expressed opportunity of purging (defaults to "Oe")
//' @param compute_O Enable computation of total opportunity of purging (false by default)
//' @param complex Enable correction for complex pedigrees.
//' @return The input dataframe, plus two additional column named "O" and "Oe", containing total and expressed opportunity of purging measures.
//' @encoding UTF-8
//' @references
//' \itemize{
//'   \item{Gulisija D, Crow JF. 2007. Inferring purging from pedigree data. Evolution 61(5): 1043-1051.}
//' }
// [[Rcpp::export]]
Rcpp::DataFrame op(Rcpp::DataFrame ped,
                   Rcpp::NumericMatrix pi,
                   Rcpp::NumericVector Fi,
                   std::string name_O,
                   std::string name_Oe,
                   bool compute_O = false,
                   bool complex = true) {
  
  // Check errors
  Rcpp::IntegerVector id = ped[id_col];
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  int N = ped.rows();
  
  // Search and map ancestors
  std::vector<std::set<int>> ancestors_ibd = map_ancestor_ibd(dam, sire, Fi); // inbred ancestors

  // Compute opportunity of purging measures
  Rcpp::NumericVector O;
  Rcpp::NumericVector E;
  bool display_progress = true;
  Rcpp::Rcerr << "Computing opportunity of purging values... " << std::endl;
  Progress p(N, display_progress);

  Rcpp::Rcout << std::endl;
  // Compute Oe (and O)
  for (int i(0); i<N; ++i) {

    double Ei = 0.0;

    // Compute O (loop over all i's inbred ancestors)
    std::set<int> ancestors_inbred;
    if (compute_O) {
      double Oi = 0.0;
      ancestors_inbred = ancestors_ibd[i];
      for (const auto& ancestor: ancestors_inbred) {
        // Skip if ancestor's inbreeding O and Oe is already accounted by a previous ancestor
        // Following Gulisija & Crow's algorithm, " the opportunity for purging in a closer
        // ancestor is discounted for what was already accounted for in its predecessor"
        // Here, we do not compute the predecessor's O and Oe if a closer ancestor already
        // contributes O and Oe, and is related to that predecessor.
        // This also helps to avoid intensive  computational burden
        bool skip = false;
        std::set<int> ancestors_other = ancestors_inbred;
        ancestors_other.erase(ancestor);
        for (std::set<int>::iterator it = ancestors_other.begin(); it != ancestors_other.end(); ++it) { // loop over i ancestors
          if (pi(*it-1, ancestor-1) > 0.0) {
            skip = true;
            break;
          }
        }
        if (skip & complex) continue;
        std::vector<int> path_n;
        map_ij_distance (id[i], ancestor, dam, sire, path_n, 1);
        double Fancestor_test (Fi[ancestor-1]);
        for (const auto& n: path_n) Oi += pow(0.5, n-1)*Fancestor_test;
      }
      O.push_back(Oi);
    }

    // Compute Oe (loop over all ancestors common to i's parents)
    std::set<int> ancestors_common; // for every individual, it list the common maternal and paternal inbred ancestors
    if (dam[i] && sire[i]) {
      std::set_intersection(ancestors_ibd[dam[i]-1].begin(),
                            ancestors_ibd[dam[i]-1].end(),
                            ancestors_ibd[sire[i]-1].begin(),
                            ancestors_ibd[sire[i]-1].end(),
                            std::inserter(ancestors_common, ancestors_common.begin()));
    }
    for (const auto& ancestor: ancestors_common) {
      // Heuristic to skip complex computations (see above)
      bool skip = false;
      std::set<int> ancestors_other = ancestors_common;
      ancestors_other.erase(ancestor);
      for (std::set<int>::iterator it = ancestors_other.begin(); it != ancestors_other.end(); ++it) { // loop over i ancestors
        if (pi(*it, ancestor-1) > 0.0) {
          skip = true;
          break;
        }
      }
      if (skip & complex) continue;
      // Compute Oe only for recent ancestors
      double Fancestor (Fi[ancestor-1]);
      Ei += pi(i, ancestor-1)*Fancestor;
    }
    E.push_back(2.0*Ei);
    p.increment(); // update progress
    Rcpp::checkUserInterrupt(); // check cancellation from user
  }

  if (compute_O) ped[name_O] = O;
  ped[name_Oe] = E;
  return ped;
}

void map_ij_distance (const int& i_id,
                      const int& j_id,
                      Rcpp::IntegerVector dam,
                      Rcpp::IntegerVector sire,
                      std::vector<int>& distance_map,
                      const int& depth) {
  
  // compute distance between individual i and its ancestor j
  // see ancestors direct descendants
  std::vector<int> descendants;
  for (int i(j_id); i < i_id; ++i) {
    if (dam[i] == j_id || sire[i] == j_id) descendants.push_back(i+1);
  }
  
  // if any direct descendant is i_idx, add distance of +1 to the map,
  // if not, loop over all remaining descendants.
  int Nd = descendants.size();
  for (int i(0); i < Nd; ++i) {
    if (descendants[i] == i_id) distance_map.push_back(depth+1);
    else map_ij_distance(i_id, descendants[i], dam, sire, distance_map, depth+1);
  }
  return;
}

std::vector<std::set<int>> map_ancestor_ibd (Rcpp::IntegerVector dam,
                                             Rcpp::IntegerVector sire,
                                             Rcpp::NumericVector Fi) {
  
  std::vector<std::set<int>> inbred_ancestors;
  
  // Map ancestors
  int N = dam.length();
  for (int i(0); i<N; ++i) {
    std::set<int> ancestors_tmp;
    if (dam[i] && Fi[dam[i]-1] > 0.0) {
      ancestors_tmp.insert(dam[i]);
      std::merge (ancestors_tmp.begin(),
                  ancestors_tmp.end(),
                  inbred_ancestors[dam[i]-1].begin(),
                  inbred_ancestors[dam[i]-1].end(),
                  std::inserter(ancestors_tmp, ancestors_tmp.begin()));
    }
    if (sire[i] && Fi[sire[i]-1] > 0.0) {
      ancestors_tmp.insert(sire[i]);
      std::merge (ancestors_tmp.begin(),
                  ancestors_tmp.end(),
                  inbred_ancestors[sire[i]-1].begin(),
                  inbred_ancestors[sire[i]-1].end(),
                  std::inserter(ancestors_tmp, ancestors_tmp.begin()));
    }
    inbred_ancestors.push_back(ancestors_tmp);
  }
  return inbred_ancestors;
}
