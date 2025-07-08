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
arma::mat my_rinvwishart(double nu, arma::mat S) {
  Function ff("rinvwishart");
  SEXP result = ff(Named("nu") = nu, _["S"] = S);
  arma::mat res = as<arma::mat>(result);
  return res;
}

// Helper function for multivariate normal sampling
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
