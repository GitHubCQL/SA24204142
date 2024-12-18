// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// GibbsC
NumericMatrix GibbsC(int N, double a, double b, int n);
RcppExport SEXP _SA24204142_GibbsC(SEXP NSEXP, SEXP aSEXP, SEXP bSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(GibbsC(N, a, b, n));
    return rcpp_result_gen;
END_RCPP
}
// MALA
List MALA(int N, NumericVector q_init, double sigma, int L, Function logTarget, Function glogTarget);
RcppExport SEXP _SA24204142_MALA(SEXP NSEXP, SEXP q_initSEXP, SEXP sigmaSEXP, SEXP LSEXP, SEXP logTargetSEXP, SEXP glogTargetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q_init(q_initSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< Function >::type logTarget(logTargetSEXP);
    Rcpp::traits::input_parameter< Function >::type glogTarget(glogTargetSEXP);
    rcpp_result_gen = Rcpp::wrap(MALA(N, q_init, sigma, L, logTarget, glogTarget));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204142_GibbsC", (DL_FUNC) &_SA24204142_GibbsC, 4},
    {"_SA24204142_MALA", (DL_FUNC) &_SA24204142_MALA, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204142(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
