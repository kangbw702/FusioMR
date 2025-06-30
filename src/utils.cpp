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
