#include <Rcpp.h>
#include "constants.h"
#include "inbreeding_utils.h"
using namespace Rcpp;

//' Deviation from Hardy-Weinberg equilibrium
//'
//' Computes the deviation from Hardy-Weinberg equilibrium following Caballero and Toro (2000).
//' 
//' @template ped-arg
//' @template reference-arg
//' @return A numeric value indicating the deviation from Hardy-Weinberg equilibrium.
//' @seealso \code{\link{pop_Ne}}
//' @references
//' \itemize{
//'   \item{Caballero A, Toro M. 2000. Interrelations between effective population size and other pedigree tools for the management of conserved populations. Genet. Res. 75: 331-343.}
//' }
// [[Rcpp::export]]
double hwd(DataFrame ped,
           Nullable<LogicalVector> reference = R_NilValue) {

  // Input values
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  int N = ped.rows();
  
  // Compute coancestry and inbreeding
  std::vector<std::vector<double>> f = complete_coancestry(ped);
  Rcpp::NumericVector F;
  for (int i(0); i<N; ++i) {
    if (dam[i] == 0 || sire[i] == 0) F.push_back(0.0);
    else if (dam[i] >= sire[i]) F.push_back(f[dam[i]-1][sire[i]-1]);
    else F.push_back(f[sire[i]-1][dam[i]-1]);
  }

  // (Optional) Filter results by reference population
  Rcpp::LogicalVector ped_ref;
  if (reference.isNotNull()) ped_ref = reference;
  int Neval = sum(ped_ref);

  // Compute Fis
  double f_ii = 0.0;
  if (reference.isNotNull()) {
    int Neval_f (0);
    for (int i(0); i<N; ++i) {
      if (!ped_ref[i]) continue;
      for (int j(0); j<(i+1); ++j) {
        if (!ped_ref[j]) continue;
        else if (i == j) {
          f_ii += f[i][i];
          ++Neval_f;
        } else {
          f_ii += 2.0*f[i][j];
          Neval_f += 2;
        }
      }
    }
    f_ii /= (double)(Neval_f);
  } else {
    for (int i(0); i<N; ++i) {
      for (int j(0); j<(i+1); ++j) {
        if (i == j) f_ii += f[i][i];
        else f_ii += 2.0*f[i][j];
      }
    }
    f_ii /= (double)(N*N);
  }

  double F_i = 0.0;
  if (reference.isNotNull()) {
    for (int i(0); i<N; ++i) {
      if (ped_ref[i]) F_i += F[i];
    }
    F_i /= (double) Neval;
  } else F_i = mean(F);

  return (F_i-f_ii) / (1.0-f_ii);
}

