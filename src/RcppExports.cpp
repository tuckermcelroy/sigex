// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// polymul_matt
arma::cube polymul_matt(arma::cube amat, arma::cube bmat);
RcppExport SEXP _sigex_polymul_matt(SEXP amatSEXP, SEXP bmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type amat(amatSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type bmat(bmatSEXP);
    rcpp_result_gen = Rcpp::wrap(polymul_matt(amat, bmat));
    return rcpp_result_gen;
END_RCPP
}
// VARMA_auto
arma::cube VARMA_auto(arma::mat param, int p, int q, int maxlag);
RcppExport SEXP _sigex_VARMA_auto(SEXP paramSEXP, SEXP pSEXP, SEXP qSEXP, SEXP maxlagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type maxlag(maxlagSEXP);
    rcpp_result_gen = Rcpp::wrap(VARMA_auto(param, p, q, maxlag));
    return rcpp_result_gen;
END_RCPP
}
// getEigenValues
arma::cx_vec getEigenValues(arma::mat M);
RcppExport SEXP _sigex_getEigenValues(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(getEigenValues(M));
    return rcpp_result_gen;
END_RCPP
}
// complexExp
arma::cx_vec complexExp(arma::vec x);
RcppExport SEXP _sigex_complexExp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(complexExp(x));
    return rcpp_result_gen;
END_RCPP
}
// polymult
arma::cx_vec polymult(arma::cx_vec a, arma::cx_vec b);
RcppExport SEXP _sigex_polymult(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::cx_vec >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(polymult(a, b));
    return rcpp_result_gen;
END_RCPP
}
// polymul_mat
arma::cube polymul_mat(arma::cube amat, arma::cube bmat);
RcppExport SEXP _sigex_polymul_mat(SEXP amatSEXP, SEXP bmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type amat(amatSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type bmat(bmatSEXP);
    rcpp_result_gen = Rcpp::wrap(polymul_mat(amat, bmat));
    return rcpp_result_gen;
END_RCPP
}
// ar_adjoint
List ar_adjoint(arma::cube poly_array);
RcppExport SEXP _sigex_ar_adjoint(SEXP poly_arraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type poly_array(poly_arraySEXP);
    rcpp_result_gen = Rcpp::wrap(ar_adjoint(poly_array));
    return rcpp_result_gen;
END_RCPP
}
// auto_VARMA
arma::cube auto_VARMA(arma::mat param, int p, int q, int ps, int qs, int season, int grid, int maxlag);
RcppExport SEXP _sigex_auto_VARMA(SEXP paramSEXP, SEXP pSEXP, SEXP qSEXP, SEXP psSEXP, SEXP qsSEXP, SEXP seasonSEXP, SEXP gridSEXP, SEXP maxlagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type ps(psSEXP);
    Rcpp::traits::input_parameter< int >::type qs(qsSEXP);
    Rcpp::traits::input_parameter< int >::type season(seasonSEXP);
    Rcpp::traits::input_parameter< int >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< int >::type maxlag(maxlagSEXP);
    rcpp_result_gen = Rcpp::wrap(auto_VARMA(param, p, q, ps, qs, season, grid, maxlag));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sigex_polymul_matt", (DL_FUNC) &_sigex_polymul_matt, 2},
    {"_sigex_VARMA_auto", (DL_FUNC) &_sigex_VARMA_auto, 4},
    {"_sigex_getEigenValues", (DL_FUNC) &_sigex_getEigenValues, 1},
    {"_sigex_complexExp", (DL_FUNC) &_sigex_complexExp, 1},
    {"_sigex_polymult", (DL_FUNC) &_sigex_polymult, 2},
    {"_sigex_polymul_mat", (DL_FUNC) &_sigex_polymul_mat, 2},
    {"_sigex_ar_adjoint", (DL_FUNC) &_sigex_ar_adjoint, 1},
    {"_sigex_auto_VARMA", (DL_FUNC) &_sigex_auto_VARMA, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_sigex(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
