#ifndef INBREEDING_UTILS_H
#define INBREEDING_UTILS_H

#include <memory>
#include "constants.h"

struct pair_hash {
  template <class T1, class T2>
  int operator() (const std::pair<T1, T2> &pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

typedef std::pair<int, int> coord;
typedef std::unordered_map<coord, std::shared_ptr<double>, pair_hash> coancestry_map;

class coancestry_matrix {
  
public:
  coancestry_map coancestry;
  double& operator[] (coord c) {
    if (c.first>c.second) c = std::make_pair (c.second, c.first);
    coancestry_map::iterator it = coancestry.find(c);
    if (it!=coancestry.end()) coancestry[c] = std::make_shared<double> (*it->second);
    else if (c.first==c.second && c.first) coancestry[c] = std::make_shared<double> (0.5);
    else coancestry[c] = std::make_shared<double> (0.0);
    return (*coancestry[c]);
  }
};

class pcoancestry_matrix { // for partial kinship coefficients
  
public:
  coancestry_map coancestry;
  double& operator[] (coord c) {
    if (c.first>c.second) c = std::make_pair (c.second, c.first);
    coancestry_map::iterator it = coancestry.find(c);
    if (it!=coancestry.end()) coancestry[c] = std::make_shared<double> (*it->second);
    else coancestry[c] = std::make_shared<double> (0.0);
    return (*coancestry[c]);
  }
};

coancestry_matrix minimum_coancestry(Rcpp::DataFrame);
std::vector<std::vector<double>> complete_coancestry(Rcpp::DataFrame);
std::vector<int> coancestry_t(Rcpp::IntegerVector, Rcpp::IntegerVector);
void search_ancestors_min_coancestry (const std::vector<int>&, const std::vector<int>&, const std::vector<int>&, const std::vector<int>&, std::set<std::pair<int, int>>&, const int& i, const int& j);

#endif
