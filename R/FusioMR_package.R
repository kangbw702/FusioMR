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
#'   \item{\code{\link{summary_stats_selected}}}{Implement instrumental variable selection with given p-value}
#' }
#'
#' @section Types of framework:
#' \describe{
#'   \item{seso_uhp_only}{Single exposure, single outcome with uncorrelated horizontal pleiotropy only}
#'   \item{seso_nohp}{Single exposure, single outcome with no horizontal pleiotropy}
#'   \item{multi_uhp_only}{Multiple exposures with uncorrelated horizontal pleiotropy (not yet implemented)}
#'   \item{multi_nohp}{Multiple exposures with no horizontal pleiotropy (not yet implemented)}
#' }
#'
#' @name FusioMR-package
#' @aliases FusioMR
#' @useDynLib FusioMR, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm pnorm quantile sd
NULL
