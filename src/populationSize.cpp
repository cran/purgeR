#include <Rcpp.h>
#include "constants.h"
#include "utils.h"
#include "error.h"
using namespace Rcpp;

void search_ancestors(Rcpp::DataFrame, int, Rcpp::LogicalVector, Rcpp::LogicalVector);
int sample_allele(int, int);

//' Individuals to be evaluated in purging analyses
//'
//' Returns a boolean vector indicating what individuals are suitable for purging analyses, given a measure of fitness.
//' Individuals with NA values of fitness, and that do not have descendants with non-NA fitness values, are excluded.
//' 
//' @template ped-arg
//' @template reference-arg
//' @param nboot Number of bootstrap iterations (for computing Ng).
//' @template seed-arg
//' @param skip_Ng Skip Ng computation or not (FALSE by default).
//' @return Boolean vector indicating what individuals will be evaluated.
// [[Rcpp::export]]
DataFrame ancestors(Rcpp::DataFrame ped,
                    Rcpp::LogicalVector reference,
                    int nboot = 10000,
                    Nullable<NumericVector> seed = R_NilValue,
                    bool skip_Ng = false) {
  
  // Check input errors
  if (seed.isNotNull()) {
    Rcpp::NumericVector s(seed);
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(s[0])));
  }

  // Setting
  Rcpp::LogicalVector ped_ref = reference;
  Rcpp::LogicalVector eval = !is_na(ped_ref);
  Rcpp::IntegerVector id = ped[id_col];
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  R_xlen_t N = ped.rows();

  // Individuals younger than the youngest in the reference population (RP) are not evaluated.
  int Neval = N;
  for (int i(N); i>=1; --i) {
    if (eval[i-1]) break;
    else --N;
  }
  
  // Compute total number of founders, ancestors, and individuals in the RP
  int Nr = 0; // RP size
  int Nf = 0; // number of founders
  int Na = 0; // number of ancestors
  Rcpp::LogicalVector founders (N, false);
  Rcpp::LogicalVector ancestors (N, false);
  for (int i(Neval); i>=1; --i) {
    if (eval[i-1] && ped_ref[i-1]) {
      ++Nr;
      search_ancestors(ped, i-1, founders, ancestors);
      if (Nr % 50 == 0) Rcpp::checkUserInterrupt(); // check cancelation from user
    }
  }
  Nf = sum(founders);
  Na = sum(ancestors);
  
  // Compute probabilities of gene origin
  Rcpp::NumericVector q (N, 0.0);
  q = ifelse (ped_ref & eval, 1.0, q);
  Rcpp::NumericVector q_init = Rcpp::clone(q);
  for (int i(N); i>=1; --i) {
    if (sire[i-1]) q[sire[i-1]-1] += 0.5*q[id[i-1]-1];
    if (dam[i-1]) q[dam[i-1]-1] += 0.5*q[id[i-1]-1];
  }
  for (int i(0); i<N; ++i) {
    if (dam[i] && !sire[i]) q[i] *= 0.5;
    else if (!dam[i] && sire[i]) q[i] *= 0.5;
  }
  q = q/Nr;
  
  // Effective number of founders
  double Nfe (0.0);
  double sum_qf = sum(ifelse(founders, q, 0.0));
  Rcpp::NumericVector qf = q/sum_qf;
  Nfe = sum(ifelse(founders, qf*qf, 0.0));
  Nfe = 1.0/Nfe;
  
  // Effective number of ancestors
  double Nae (0.0);
  double maxq (0.0);
  int major (N-1);
  for (int i(0); i<N; ++i) {
    if (ancestors[i] && q[i] > maxq) {
      maxq = q[i];
      major = i;
    }
  }
  
  Rcpp::DataFrame copy_ped = Rcpp::clone(ped);
  Rcpp::IntegerVector copy_dam = copy_ped[dam_col];
  Rcpp::IntegerVector copy_sire = copy_ped[sire_col];
  std::vector<bool> ancestors_std = Rcpp::as<std::vector<bool>>(ancestors);
  copy_dam[major] = 0;
  copy_sire[major] = 0;
  ancestors_std[major] = 0;
  std::vector<double> pf (N, 0.0);
  std::vector<bool> used_a (N, false);
  pf[major] = q[major];
  used_a[major] = true;
  int count (0);

  while (is_in(true, ancestors_std)) {

    // control errors
    if (++count > N) {
      Rcpp::Rcerr << ERROR_UNEXPECTED << std::endl;
      return Rcpp::DataFrame::create(
        Named("Nr") = Nr,
        Named("Nf") = Nf,
        Named("Nfe") = Nfe,
        Named("Na") = Na,
        Named("Nae") = NA_REAL,
        Named("Ng") = NA_REAL,
        Named("se_Ng") = NA_REAL
      );
    }

    // compute q
    std::vector<double> tmp_q = Rcpp::as<std::vector<double>>(q_init);
    for (int i(N); i>=1; --i) {
      if (used_a[i-1]) continue;
      if (sire[i-1]) tmp_q[sire[i-1]-1] += 0.5*tmp_q[id[i-1]-1];
      if (dam[i-1]) tmp_q[dam[i-1]-1] += 0.5*tmp_q[id[i-1]-1];
    }
    // compute a
    std::vector<double> a (N, 1.0);
    for (int i(0); i<N; ++i) if (ancestors_std[i]) a[i] = 0.0;
    for (int i(0); i<N; ++i) {
      if (used_a[i]) continue;
      if ((sire[i]-1) && a[sire[i]-1] > 0.0) a[id[i]-1] += 0.5*a[sire[i]-1];
      if ((dam[i]-1) && a[dam[i]-1] > 0.0) a[id[i]-1] += 0.5*a[dam[i]-1];
    }
    // compute p
    std::vector<double> p (tmp_q);
    for (int i(0); i<N; ++i) { p[i] *= (1.0-a[i]); }
    // choose new major
    double maxp (-1.0);
    for (int i(N); i>=1; --i) {
      if (!ancestors_std[i-1]) continue;
      if (p[i-1] > maxp) { maxp = p[i-1]; major = i-1; }
    }
    // final f
    pf[major] = p[major];
    pf[major] /= (double)Nr;
    // clean ped
    copy_dam[major] = 0;
    copy_sire[major] = 0;
    ancestors_std[major] = false;
    used_a[major] = true;

    // check cancelation from user
    Rcpp::checkUserInterrupt();
  }
  for (int i(0); i<N; ++i) {
    if (ancestors[i]) Nae += (pf[i]*pf[i]);
  }
  Nae = 1.0 / Nae;

  // Premature output
  if (skip_Ng) {
    return Rcpp::DataFrame::create(
      Named("Nr") = Nr,
      Named("Nf") = Nf,
      Named("Nfe") = Nfe,
      Named("Na") = Na,
      Named("Nae") = Nae,
      Named("Ng") = NA_REAL,
      Named("se_Ng") = NA_REAL
    );
  }

  // Effective number of founder genomes
  copy_ped = Rcpp::clone(ped);
  Rcpp::NumericVector v_Ng;
  for (int b(0); b < nboot; ++b) {
    // simulate segregation
    Rcpp::IntegerVector abs_freq;
    Rcpp::DataFrame tmp_ped = Rcpp::clone(ped);
    Rcpp::IntegerVector tmp_dam = copy_ped[dam_col];
    Rcpp::IntegerVector tmp_sire = copy_ped[sire_col];
    std::vector<int> al_dam;
    std::vector<int> al_sire;
    int al_count(0);
    for (int i(0); i<N; ++i) {
      if (!tmp_dam[i] && !tmp_sire[i]) {
        al_dam.push_back(al_count++);
        al_sire.push_back(al_count++);
        abs_freq.push_back(0);
        abs_freq.push_back(0);
      } else if (!tmp_dam[i]) {
        abs_freq.push_back(0);
        if (ped_ref[i] & eval[i]) ++abs_freq[abs_freq.size()-1];
        al_dam.push_back(al_count++);
        int al (sample_allele(al_dam[sire[i]-1], al_sire[sire[i]-1]));
        if (ped_ref[i] & eval[i]) ++abs_freq[al];
        al_sire.push_back(al);
      } else if (!tmp_sire[i]) {
        int al (sample_allele(al_dam[dam[i]-1], al_sire[dam[i]-1]));
        al_dam.push_back(al);
        if (ped_ref[i] & eval[i]) ++abs_freq[al];
        abs_freq.push_back(0);
        if (ped_ref[i] & eval[i]) ++abs_freq[abs_freq.size()-1];
        al_sire.push_back(al_count++);
      } else {
        int al1 (sample_allele(al_dam[dam[i]-1], al_sire[dam[i]-1]));
        int al2 (sample_allele(al_dam[sire[i]-1], al_sire[sire[i]-1]));
        al_dam.push_back(al1);
        al_sire.push_back(al2);
        if (ped_ref[i] & eval[i]) {
          ++abs_freq[al1];
          ++abs_freq[al2];
        }
      }
    }
    Rcpp::NumericVector rel_freq;
    int n2f = abs_freq.size();
    for (int i(0); i<n2f; ++i) rel_freq.push_back((double)abs_freq[i]/(2.0*Nr));
    double Ngi (0.0);
    for (int i(0); i<n2f; ++i) Ngi += (rel_freq[i]*rel_freq[i]);
    Ngi = 1.0 / (2.0*Ngi);
    v_Ng.push_back(Ngi);
  }
  double Ng (0.0);
  double Ng_se (0.0);
  for (const auto& i: v_Ng) Ng += i;
  Ng /= double(nboot);
  for (const auto& i: v_Ng) Ng_se += (i-Ng)*(i-Ng);
  Ng_se = sqrt(Ng_se);
  Ng_se /= sqrt(double(nboot));

  // Output
  return Rcpp::DataFrame::create(
    Named("Nr") = Nr,
    Named("Nf") = Nf,
    Named("Nfe") = Nfe,
    Named("Na") = Na,
    Named("Nae") = Nae,
    Named("Ng") = Ng,
    Named("se_Ng") = Ng_se
  );
}

//' Search and individuals' ancestors
//'
//' Recursive function that gathers all founders and ancestors for a given individual
//' 
//' @name search_ancestors
//' @template ped-arg
//' @param i Reference individual (its index, not id).
//' @param fnd Vector of founders (to be returned as reference).
//' @param anc Vector of ancestors (to be returned as reference).
//' @return The sampled allele.
void search_ancestors(Rcpp::DataFrame ped,
                      int i,
                      Rcpp::LogicalVector fnd,
                      Rcpp::LogicalVector anc) {
  
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  int idam = dam[i];
  int isire = sire[i];
  if (!idam && !isire) {
    fnd[i] = true;
    anc[i] = true;
    return;
  } else {
    if (idam && !anc[idam-1]) {
      anc[idam-1] = true;
      search_ancestors(ped, idam-1, fnd, anc);
    }
    if (isire && !anc[isire-1]) {
      anc[isire-1] = true;
      search_ancestors(ped, isire-1, fnd, anc);
    }
  }
}

//' Sample dam or sire inherited allele
//'
//' Given two alleles (one from dam, the other from sire), it samples one at random.
//' 
//' @name sample_allele
//' @param dam_al Dam allele.
//' @param sire_al Sire allele.
//' @return The sampled allele.
int sample_allele (int dam_al, int sire_al) {
  double random = R::runif(0.0, 1.0);
  if (random < 0.5) return dam_al;
  else return sire_al;
}

