#include <Rcpp.h>

//' Simulate alleles for the genedrop method
//'
//' Depending on the parental id and genotype, allele segregation is simulated.
//' 
//' @name genedrop_sim_allele
//' @param pid Parental identity.
//' @param dam_alleles Vector of saved maternal alleles.
//' @param sire_alleles Vector of saved paternal alleles.
//' @param allele_count number identity of the last new allele created.
//' @return A new or sampled allele.
int genedrop_sim_allele (int pid,
                         Rcpp::IntegerVector dam_alleles,
                         Rcpp::IntegerVector sire_alleles,
                         int& allele_count) {
  
  /* LOCAL VARIABLES */
    int allele (0); // simulated allele
  
  // inheritance if the parental is unknown
  if (!pid) {
    allele = allele_count;
    ++allele_count;
    // inheritance if the parental is homozygous
  } else if (dam_alleles[pid-1] == sire_alleles[pid-1]) {
    return dam_alleles[pid-1];
    // inheritance if the parental is heterozygous 
  } else {
    double random = R::runif(0.0, 1.0);
    if (random < 0.5) return dam_alleles[pid-1];
    else return sire_alleles[pid-1];
  }
  return allele;
}

//' Keep track of IBD alleles
//'
//' Use vectors on genedrop-inherited alleles to annotate IBD events.
//' 
//' @name genedrop_ibd
//' @param pid Parental identity.
//' @param dam_alleles Vector of saved maternal alleles.
//' @param sire_alleles Vector of saved paternal alleles.
//' @param allele_count number identity of the last new allele created.
//' @return Vector of maternal and paternal alleles that were in ancestral IBD (passed by reference).
void genedrop_ibd (Rcpp::IntegerVector dam,
                   Rcpp::IntegerVector sire,
                   Rcpp::IntegerVector dam_alleles,
                   Rcpp::IntegerVector sire_alleles,
                   Rcpp::LogicalVector& ibd_dam_alleles,
                   Rcpp::LogicalVector& ibd_sire_alleles) {
  
  int N = Rf_length(dam);
  for (int i(0); i<N; ++i) {
    
    int i_d = dam[i];
    int i_s = sire[i];
    
    // maternal allele
    // if i's mother is unknown, allele is not IBD
    if (i_d == 0) {
      ibd_dam_alleles.push_back(false);
      // if i's mother is homozygous, allele is IBD
  } else if (dam_alleles[i_d-1] == sire_alleles[i_d-1]) {
    ibd_dam_alleles.push_back(true);
    // if i's maternal allele is i's mother maternal allele, and it was IBD
  } else if ((dam_alleles[i] == dam_alleles[i_d-1]) && ibd_dam_alleles[i_d-1]) {
    ibd_dam_alleles.push_back(true);
    // if i's maternal allele is i's mother paternal allele, and it was IBD
  } else if ((dam_alleles[i] == sire_alleles[i_d-1]) && ibd_sire_alleles[i_d-1]) {
    ibd_dam_alleles.push_back(true);
  } else {
    ibd_dam_alleles.push_back(false);
  }
  
  // idem for paternal allele
  if (i_s == 0) {
    ibd_sire_alleles.push_back(false);
  } else if (dam_alleles[i_s-1] == sire_alleles[i_s-1]) {
    ibd_sire_alleles.push_back(true);
  } else if ((sire_alleles[i] == dam_alleles[i_s-1]) && ibd_dam_alleles[i_s-1]) {
    ibd_sire_alleles.push_back(true);
  } else if ((sire_alleles[i] == sire_alleles[i_s-1]) && ibd_sire_alleles[i_s-1]) {
    ibd_sire_alleles.push_back(true);
  } else {
    ibd_sire_alleles.push_back(false);
  }
}
}

//' Read ancestral inbreeding from IBD alleles
//'
//' Return ancestral inbreeding measures from vectors indicating IBD of sampled alleles
//' 
//' @name genedrop_read_Fa
//' @param pid Parental identity.
//' @param dam_alleles Vector of saved maternal alleles.
//' @param sire_alleles Vector of saved paternal alleles.
//' @param allele_count number identity of the last new allele created.
//' @return Vector of maternal and paternal alleles that were in ancestral IBD (passed by reference).
Rcpp::NumericVector genedrop_read_Fa (Rcpp::IntegerVector dam,
                                      Rcpp::IntegerVector sire,
                                      Rcpp::IntegerVector dam_alleles,
                                      Rcpp::IntegerVector sire_alleles,
                                      Rcpp::LogicalVector ibd_dam_alleles,
                                      Rcpp::LogicalVector ibd_sire_alleles) {
  
  int N = Rf_length(ibd_dam_alleles);
  Rcpp::NumericVector Fa (N, 0.0);
  
  /* INDIVIDUALS CARRYING ALLELES THAT WHERE IBD IN AN ANCESTOR HAVE Fa */
    for (int i(0); i<N; ++i) {
      
      int i_d = dam[i];
      int i_s = sire[i];
      
      // if i's mother is homozygous, or has Fa=1, the allele inherited from her is IBD
    if (i_d && ((dam_alleles[i_d-1] == sire_alleles[i_d-1]) || Fa[i_d-1] == 1.0 )) {
      Fa[i] += 0.5;
      // if i's mother is heterozygous, and only one of her alleles was IBD
    } else if (i_d && Fa[i_d-1] == 0.5) {
      // if the allele is the maternal allele of the mother, and is IBD
      if (ibd_dam_alleles[i_d-1] && dam_alleles[i_d-1] == dam_alleles[i]) {
        Fa[i] += 0.5;
        // if the allele is the paternal allele of the mother, and is IBD
      } else if (ibd_sire_alleles[i_d-1] && sire_alleles[i_d-1] == dam_alleles[i]) {
        Fa[i] += 0.5;
      }
    }
  
  // idem for the father
  if (i_s && ((sire_alleles[i_s-1] == dam_alleles[i_s-1]) || Fa[i_s-1] == 1.0 )) {
    Fa[i] += 0.5;
  } else if (i_s && Fa[i_s-1] == 0.5) {
    if (ibd_sire_alleles[i_s-1] && sire_alleles[i_s-1] == sire_alleles[i]) {
      Fa[i] += 0.5;
    } else if (ibd_dam_alleles[i_s-1] && dam_alleles[i_s-1] == dam_alleles[i]) {
      Fa[i] += 0.5;
    }
  }
}
return Fa;
}
