// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// reproductive_value
Rcpp::DataFrame reproductive_value(Rcpp::DataFrame ped, Rcpp::LogicalVector reference, std::string name_to, Nullable<Rcpp::LogicalVector> target, bool enable_correction);
RcppExport SEXP _purgeR_reproductive_value(SEXP pedSEXP, SEXP referenceSEXP, SEXP name_toSEXP, SEXP targetSEXP, SEXP enable_correctionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type reference(referenceSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_to(name_toSEXP);
    Rcpp::traits::input_parameter< Nullable<Rcpp::LogicalVector> >::type target(targetSEXP);
    Rcpp::traits::input_parameter< bool >::type enable_correction(enable_correctionSEXP);
    rcpp_result_gen = Rcpp::wrap(reproductive_value(ped, reference, name_to, target, enable_correction));
    return rcpp_result_gen;
END_RCPP
}
// F
DataFrame F(DataFrame ped, std::string name_to);
RcppExport SEXP _purgeR_F(SEXP pedSEXP, SEXP name_toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_to(name_toSEXP);
    rcpp_result_gen = Rcpp::wrap(F(ped, name_to));
    return rcpp_result_gen;
END_RCPP
}
// Fa
DataFrame Fa(DataFrame ped, Rcpp::NumericVector Fi, std::string name_to, int genedrop, Nullable<int> seed);
RcppExport SEXP _purgeR_Fa(SEXP pedSEXP, SEXP FiSEXP, SEXP name_toSEXP, SEXP genedropSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Fi(FiSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_to(name_toSEXP);
    Rcpp::traits::input_parameter< int >::type genedrop(genedropSEXP);
    Rcpp::traits::input_parameter< Nullable<int> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(Fa(ped, Fi, name_to, genedrop, seed));
    return rcpp_result_gen;
END_RCPP
}
// g
DataFrame g(Rcpp::DataFrame ped, double d, Rcpp::NumericVector Fi, std::string name_to);
RcppExport SEXP _purgeR_g(SEXP pedSEXP, SEXP dSEXP, SEXP FiSEXP, SEXP name_toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Fi(FiSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_to(name_toSEXP);
    rcpp_result_gen = Rcpp::wrap(g(ped, d, Fi, name_to));
    return rcpp_result_gen;
END_RCPP
}
// hwd
double hwd(DataFrame ped, Nullable<LogicalVector> reference);
RcppExport SEXP _purgeR_hwd(SEXP pedSEXP, SEXP referenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< Nullable<LogicalVector> >::type reference(referenceSEXP);
    rcpp_result_gen = Rcpp::wrap(hwd(ped, reference));
    return rcpp_result_gen;
END_RCPP
}
// op
Rcpp::DataFrame op(Rcpp::DataFrame ped, Rcpp::NumericMatrix pi, Rcpp::NumericVector Fi, std::string name_O, std::string name_Oe, std::string sufix, bool compute_O);
RcppExport SEXP _purgeR_op(SEXP pedSEXP, SEXP piSEXP, SEXP FiSEXP, SEXP name_OSEXP, SEXP name_OeSEXP, SEXP sufixSEXP, SEXP compute_OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pi(piSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Fi(FiSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_O(name_OSEXP);
    Rcpp::traits::input_parameter< std::string >::type name_Oe(name_OeSEXP);
    Rcpp::traits::input_parameter< std::string >::type sufix(sufixSEXP);
    Rcpp::traits::input_parameter< bool >::type compute_O(compute_OSEXP);
    rcpp_result_gen = Rcpp::wrap(op(ped, pi, Fi, name_O, name_Oe, sufix, compute_O));
    return rcpp_result_gen;
END_RCPP
}
// Fij_core_i_cpp
Rcpp::NumericVector Fij_core_i_cpp(Rcpp::IntegerVector dam, Rcpp::IntegerVector sire, const int& anc_idx, Rcpp::LogicalVector mapa, Rcpp::NumericVector Fi, int genedrop, Rcpp::Nullable<int> seed);
RcppExport SEXP _purgeR_Fij_core_i_cpp(SEXP damSEXP, SEXP sireSEXP, SEXP anc_idxSEXP, SEXP mapaSEXP, SEXP FiSEXP, SEXP genedropSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dam(damSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type sire(sireSEXP);
    Rcpp::traits::input_parameter< const int& >::type anc_idx(anc_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type mapa(mapaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Fi(FiSEXP);
    Rcpp::traits::input_parameter< int >::type genedrop(genedropSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(Fij_core_i_cpp(dam, sire, anc_idx, mapa, Fi, genedrop, seed));
    return rcpp_result_gen;
END_RCPP
}
// ancestors
DataFrame ancestors(Rcpp::DataFrame ped, Rcpp::LogicalVector reference, Rcpp::IntegerVector rp_idx, int nboot, Nullable<NumericVector> seed, bool skip_Ng);
RcppExport SEXP _purgeR_ancestors(SEXP pedSEXP, SEXP referenceSEXP, SEXP rp_idxSEXP, SEXP nbootSEXP, SEXP seedSEXP, SEXP skip_NgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type reference(referenceSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type rp_idx(rp_idxSEXP);
    Rcpp::traits::input_parameter< int >::type nboot(nbootSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type skip_Ng(skip_NgSEXP);
    rcpp_result_gen = Rcpp::wrap(ancestors(ped, reference, rp_idx, nboot, seed, skip_Ng));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_purgeR_reproductive_value", (DL_FUNC) &_purgeR_reproductive_value, 5},
    {"_purgeR_F", (DL_FUNC) &_purgeR_F, 2},
    {"_purgeR_Fa", (DL_FUNC) &_purgeR_Fa, 5},
    {"_purgeR_g", (DL_FUNC) &_purgeR_g, 4},
    {"_purgeR_hwd", (DL_FUNC) &_purgeR_hwd, 2},
    {"_purgeR_op", (DL_FUNC) &_purgeR_op, 7},
    {"_purgeR_Fij_core_i_cpp", (DL_FUNC) &_purgeR_Fij_core_i_cpp, 7},
    {"_purgeR_ancestors", (DL_FUNC) &_purgeR_ancestors, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_purgeR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
