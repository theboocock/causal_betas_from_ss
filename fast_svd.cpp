#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec baseSVD(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "standard");
  return S;
}

// [[Rcpp::export]]
arma::vec dcSVD(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "dc");
  return S;
}

