#ifndef GENEDROP_H
#define GENEDROP_H

int genedrop_sim_allele (int, Rcpp::IntegerVector, Rcpp::IntegerVector, int&);
void genedrop_ibd (Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::LogicalVector&, Rcpp::LogicalVector&);
Rcpp::NumericVector genedrop_read_Fa (Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::LogicalVector, Rcpp::LogicalVector);

#endif
