// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// functionalBART_res
Rcpp::List functionalBART_res(Rcpp::List& x, Rcpp::NumericVector y, size_t p, size_t n, size_t funcP, size_t intvls, Rcpp::List& x_test, size_t n_test, double nu, double lambda, const Rcpp::List& iXinfo, Rcpp::IntegerVector& nc, size_t nkeeptrain, double binaryOffset, int iter, int burn, double alpha, double beta, double tau, bool dart, bool dart_fp, double theta_fp, double omega_fp, double a_fp, double b_fp, double rho_fp, bool aug_fp, bool dart_int, double theta_int, double omega_int, double a_int, double b_int, double rho_int, bool aug_int);
RcppExport SEXP _FuncBART_functionalBART_res(SEXP xSEXP, SEXP ySEXP, SEXP pSEXP, SEXP nSEXP, SEXP funcPSEXP, SEXP intvlsSEXP, SEXP x_testSEXP, SEXP n_testSEXP, SEXP nuSEXP, SEXP lambdaSEXP, SEXP iXinfoSEXP, SEXP ncSEXP, SEXP nkeeptrainSEXP, SEXP binaryOffsetSEXP, SEXP iterSEXP, SEXP burnSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP, SEXP dartSEXP, SEXP dart_fpSEXP, SEXP theta_fpSEXP, SEXP omega_fpSEXP, SEXP a_fpSEXP, SEXP b_fpSEXP, SEXP rho_fpSEXP, SEXP aug_fpSEXP, SEXP dart_intSEXP, SEXP theta_intSEXP, SEXP omega_intSEXP, SEXP a_intSEXP, SEXP b_intSEXP, SEXP rho_intSEXP, SEXP aug_intSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< size_t >::type p(pSEXP);
    Rcpp::traits::input_parameter< size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< size_t >::type funcP(funcPSEXP);
    Rcpp::traits::input_parameter< size_t >::type intvls(intvlsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type x_test(x_testSEXP);
    Rcpp::traits::input_parameter< size_t >::type n_test(n_testSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type iXinfo(iXinfoSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< size_t >::type nkeeptrain(nkeeptrainSEXP);
    Rcpp::traits::input_parameter< double >::type binaryOffset(binaryOffsetSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< bool >::type dart(dartSEXP);
    Rcpp::traits::input_parameter< bool >::type dart_fp(dart_fpSEXP);
    Rcpp::traits::input_parameter< double >::type theta_fp(theta_fpSEXP);
    Rcpp::traits::input_parameter< double >::type omega_fp(omega_fpSEXP);
    Rcpp::traits::input_parameter< double >::type a_fp(a_fpSEXP);
    Rcpp::traits::input_parameter< double >::type b_fp(b_fpSEXP);
    Rcpp::traits::input_parameter< double >::type rho_fp(rho_fpSEXP);
    Rcpp::traits::input_parameter< bool >::type aug_fp(aug_fpSEXP);
    Rcpp::traits::input_parameter< bool >::type dart_int(dart_intSEXP);
    Rcpp::traits::input_parameter< double >::type theta_int(theta_intSEXP);
    Rcpp::traits::input_parameter< double >::type omega_int(omega_intSEXP);
    Rcpp::traits::input_parameter< double >::type a_int(a_intSEXP);
    Rcpp::traits::input_parameter< double >::type b_int(b_intSEXP);
    Rcpp::traits::input_parameter< double >::type rho_int(rho_intSEXP);
    Rcpp::traits::input_parameter< bool >::type aug_int(aug_intSEXP);
    rcpp_result_gen = Rcpp::wrap(functionalBART_res(x, y, p, n, funcP, intvls, x_test, n_test, nu, lambda, iXinfo, nc, nkeeptrain, binaryOffset, iter, burn, alpha, beta, tau, dart, dart_fp, theta_fp, omega_fp, a_fp, b_fp, rho_fp, aug_fp, dart_int, theta_int, omega_int, a_int, b_int, rho_int, aug_int));
    return rcpp_result_gen;
END_RCPP
}
// functionalBART
Rcpp::List functionalBART(Rcpp::List& x, Rcpp::NumericVector y, size_t p, size_t n, size_t funcP, size_t intvls, Rcpp::List& x_test, size_t n_test, const Rcpp::List& iXinfo, Rcpp::IntegerVector& nc, size_t nkeeptrain, double binaryOffset, int iter, int burn, double alpha, double beta, double tau, bool dart, double theta, double omega, double a, double b, double rho, bool aug);
RcppExport SEXP _FuncBART_functionalBART(SEXP xSEXP, SEXP ySEXP, SEXP pSEXP, SEXP nSEXP, SEXP funcPSEXP, SEXP intvlsSEXP, SEXP x_testSEXP, SEXP n_testSEXP, SEXP iXinfoSEXP, SEXP ncSEXP, SEXP nkeeptrainSEXP, SEXP binaryOffsetSEXP, SEXP iterSEXP, SEXP burnSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP, SEXP dartSEXP, SEXP thetaSEXP, SEXP omegaSEXP, SEXP aSEXP, SEXP bSEXP, SEXP rhoSEXP, SEXP augSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< size_t >::type p(pSEXP);
    Rcpp::traits::input_parameter< size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< size_t >::type funcP(funcPSEXP);
    Rcpp::traits::input_parameter< size_t >::type intvls(intvlsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type x_test(x_testSEXP);
    Rcpp::traits::input_parameter< size_t >::type n_test(n_testSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type iXinfo(iXinfoSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< size_t >::type nkeeptrain(nkeeptrainSEXP);
    Rcpp::traits::input_parameter< double >::type binaryOffset(binaryOffsetSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< bool >::type dart(dartSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type aug(augSEXP);
    rcpp_result_gen = Rcpp::wrap(functionalBART(x, y, p, n, funcP, intvls, x_test, n_test, iXinfo, nc, nkeeptrain, binaryOffset, iter, burn, alpha, beta, tau, dart, theta, omega, a, b, rho, aug));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FuncBART_functionalBART_res", (DL_FUNC) &_FuncBART_functionalBART_res, 34},
    {"_FuncBART_functionalBART", (DL_FUNC) &_FuncBART_functionalBART, 24},
    {NULL, NULL, 0}
};

RcppExport void R_init_FuncBART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
