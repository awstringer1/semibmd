// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// knotindex
int knotindex(double x, Eigen::VectorXd t);
RcppExport SEXP _semibmd_knotindex(SEXP xSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(knotindex(x, t));
    return rcpp_result_gen;
END_RCPP
}
// deBoor
double deBoor(double x, int k, Eigen::VectorXd t, Eigen::VectorXd beta, int p);
RcppExport SEXP _semibmd_deBoor(SEXP xSEXP, SEXP kSEXP, SEXP tSEXP, SEXP betaSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(deBoor(x, k, t, beta, p));
    return rcpp_result_gen;
END_RCPP
}
// reflect
double reflect(double xt, double lb, double ub);
RcppExport SEXP _semibmd_reflect(SEXP xtSEXP, SEXP lbSEXP, SEXP ubSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type xt(xtSEXP);
    Rcpp::traits::input_parameter< double >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< double >::type ub(ubSEXP);
    rcpp_result_gen = Rcpp::wrap(reflect(xt, lb, ub));
    return rcpp_result_gen;
END_RCPP
}
// bounded_newton
double bounded_newton(Function g, Function gprime, Eigen::VectorXd bounds, double xt, double eps, int maxitr);
RcppExport SEXP _semibmd_bounded_newton(SEXP gSEXP, SEXP gprimeSEXP, SEXP boundsSEXP, SEXP xtSEXP, SEXP epsSEXP, SEXP maxitrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type g(gSEXP);
    Rcpp::traits::input_parameter< Function >::type gprime(gprimeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< double >::type xt(xtSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxitr(maxitrSEXP);
    rcpp_result_gen = Rcpp::wrap(bounded_newton(g, gprime, bounds, xt, eps, maxitr));
    return rcpp_result_gen;
END_RCPP
}
// get_gamma
Eigen::VectorXd get_gamma(Eigen::VectorXd beta);
RcppExport SEXP _semibmd_get_gamma(SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_gamma(beta));
    return rcpp_result_gen;
END_RCPP
}
// Ux_cpp
double Ux_cpp(double x, Eigen::VectorXd beta, Eigen::VectorXd knots, int k, double fx0, double sigmaest, double A);
RcppExport SEXP _semibmd_Ux_cpp(SEXP xSEXP, SEXP betaSEXP, SEXP knotsSEXP, SEXP kSEXP, SEXP fx0SEXP, SEXP sigmaestSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type fx0(fx0SEXP);
    Rcpp::traits::input_parameter< double >::type sigmaest(sigmaestSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(Ux_cpp(x, beta, knots, k, fx0, sigmaest, A));
    return rcpp_result_gen;
END_RCPP
}
// Uxd_cpp
double Uxd_cpp(double x, Eigen::VectorXd beta, Eigen::VectorXd knots, int k, double sigmaest);
RcppExport SEXP _semibmd_Uxd_cpp(SEXP xSEXP, SEXP betaSEXP, SEXP knotsSEXP, SEXP kSEXP, SEXP sigmaestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaest(sigmaestSEXP);
    rcpp_result_gen = Rcpp::wrap(Uxd_cpp(x, beta, knots, k, sigmaest));
    return rcpp_result_gen;
END_RCPP
}
// Vx_cpp
double Vx_cpp(double x, Eigen::MatrixXd Vbeta, Eigen::VectorXd knots, int k, double x0, double sigmaest);
RcppExport SEXP _semibmd_Vx_cpp(SEXP xSEXP, SEXP VbetaSEXP, SEXP knotsSEXP, SEXP kSEXP, SEXP x0SEXP, SEXP sigmaestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Vbeta(VbetaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type sigmaest(sigmaestSEXP);
    rcpp_result_gen = Rcpp::wrap(Vx_cpp(x, Vbeta, knots, k, x0, sigmaest));
    return rcpp_result_gen;
END_RCPP
}
// get_bmd_cpp
double get_bmd_cpp(Eigen::VectorXd beta, Eigen::VectorXd knots, Eigen::VectorXd bounds, double x0, double sigmaest, double A, double eps, int maxitr);
RcppExport SEXP _semibmd_get_bmd_cpp(SEXP betaSEXP, SEXP knotsSEXP, SEXP boundsSEXP, SEXP x0SEXP, SEXP sigmaestSEXP, SEXP ASEXP, SEXP epsSEXP, SEXP maxitrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type sigmaest(sigmaestSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxitr(maxitrSEXP);
    rcpp_result_gen = Rcpp::wrap(get_bmd_cpp(beta, knots, bounds, x0, sigmaest, A, eps, maxitr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_semibmd_knotindex", (DL_FUNC) &_semibmd_knotindex, 2},
    {"_semibmd_deBoor", (DL_FUNC) &_semibmd_deBoor, 5},
    {"_semibmd_reflect", (DL_FUNC) &_semibmd_reflect, 3},
    {"_semibmd_bounded_newton", (DL_FUNC) &_semibmd_bounded_newton, 6},
    {"_semibmd_get_gamma", (DL_FUNC) &_semibmd_get_gamma, 1},
    {"_semibmd_Ux_cpp", (DL_FUNC) &_semibmd_Ux_cpp, 7},
    {"_semibmd_Uxd_cpp", (DL_FUNC) &_semibmd_Uxd_cpp, 5},
    {"_semibmd_Vx_cpp", (DL_FUNC) &_semibmd_Vx_cpp, 6},
    {"_semibmd_get_bmd_cpp", (DL_FUNC) &_semibmd_get_bmd_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_semibmd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
