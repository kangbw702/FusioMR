#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

// Function declarations
Rcpp::NumericVector my_rinvgamma(int n, double shape, double rate);
Rcpp::NumericVector my_rdirichlet(int n, Rcpp::NumericVector alpha);

#endif
