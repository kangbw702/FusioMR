#' @keywords internal
"_PACKAGE"

#' FusioMR: Fusion Methods for Mendelian Randomization
#'
#' @description
#' The FusioMR package implements FusioMR framework for Mendelian Randomization
#' analysis using Bayesian approaches with Gibbs sampling. It provides
#' functions for instrumental variable selection and causal effect estimation
#' from GWAS summary statistics.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{fusiomr}}}{Main function for causal effect estimation with given summary statistics}
#' }
#'
#' @name FusioMR-package
#' @aliases FusioMR
#' @useDynLib FusioMR, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm pnorm quantile sd
NULL
