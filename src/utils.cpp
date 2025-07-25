#include "utils.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Implementation of my_rinvgamma
NumericVector my_rinvgamma(int n, double shape, double rate) {
  Function ff("rinvgamma", Environment::namespace_env("invgamma"));
  NumericVector res = ff(n, Named("shape") = shape, _["rate"] = rate);
  return res;
}

// Implementation of my_rinvwishart
// [[Rcpp::export]]
arma::mat my_rinvwishart(double nu, arma::mat S) {
  return arma::iwishrnd(S, nu);
}

// Helper function for multivariate normal sampling with robust Cholesky
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  sigma = 0.5 * (sigma + sigma.t());
  arma::mat L;
  bool success = arma::chol(L, sigma, "lower");
  if (!success) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, sigma);
    eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
    L = eigvec * arma::diagmat(arma::sqrt(eigval));
  }
  return arma::repmat(mu, 1, n).t() + Y * L.t();
}

// [[Rcpp::export]]
Rcpp::NumericVector my_rdirichlet(int n, Rcpp::NumericVector alpha) {
  if (n != 1) Rcpp::stop("my_rdirichlet currently supports n = 1 only.");
  int K = alpha.size();
  Rcpp::NumericVector draws(K);
  double sum = 0.0;
  for (int k = 0; k < K; ++k) {
    draws[k] = R::rgamma(alpha[k], 1.0);
    sum += draws[k];
  }
  // Normalize
  for (int k = 0; k < K; ++k) draws[k] /= sum;
  return draws;
}

