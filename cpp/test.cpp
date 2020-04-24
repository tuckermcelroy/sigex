// Use the RcppArmadillo package 
// Requires different header file from Rcpp.h 
#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec getRe(arma::cx_vec x) {
  return arma::real(x);
}
