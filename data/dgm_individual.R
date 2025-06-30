#' Generate Simulated GWAS Data
#'
#' @description
#' This function generates simulated individual-level data for two-sample
#' Mendelian Randomization analysis with potential horizontal pleiotropy.
#'
#' @param m Integer. Number of genetic variants (SNPs) to simulate.
#' @param nx Integer. Sample size for the exposure study.
#' @param ny Integer. Sample size for the outcome study.
#' @param a_gamma Numeric. Lower bound for IV-to-exposure effect sizes.
#' @param b_gamma Numeric. Upper bound for IV-to-exposure effect sizes.
#' @param a_f Numeric. Lower bound for MAF.
#' @param b_f Numeric. Upper bound for MAF.
#' @param a_alpha Numeric. Lower bound for horizontal pleiotropy effects.
#' @param b_alpha Numeric. Upper bound for horizontal pleiotropy effects.
#' @param theta Numeric. True causal effect of exposure on outcome.
#'
#' @details
#' The data generation model follows the structure:
#' \deqn{X = \gamma^T G_X + \epsilon_X}
#' \deqn{Y = \alpha^T G_Y + \theta X_Y + \epsilon_Y}
#'
#' where:
#' \itemize{
#'   \item \eqn{G_X, G_Y} are standardized genotype matrices for exposure and outcome studies
#'   \item \eqn{\gamma} represents valid IV-to-exposure effects (instrument strength)
#'   \item \eqn{\alpha} represents horizontal pleiotropy effects
#'   \item \eqn{\theta} is the causal effect
#'   \item \eqn{\epsilon_X, \epsilon_Y} are standard normal error terms
#' }
#'
#' Genotypes are simulated from binomial distributions with varying MAFs,
#' then standardized to have mean 0 and variance 1 for each SNP.
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{Numeric vector of length \code{nx}. Simulated exposure values for exposure study.}
#'   \item{Y}{Numeric vector of length \code{ny}. Simulated outcome values for outcome study.}
#'   \item{gx}{Numeric matrix of dimension \code{m x nx}. Standardized genotype matrix for exposure study.}
#'   \item{gy}{Numeric matrix of dimension \code{m x ny}. Standardized genotype matrix for outcome study.}
#'   \item{gamma}{Numeric vector of length \code{m}. True IV-to-exposure effect sizes.}
#'   \item{alpha}{Numeric vector of length \code{m}. True horizontal pleiotropy effects.}
#'   \item{f}{Numeric vector of length \code{m}. Minor allele frequencies used for simulation.}
#'   \item{theta}{Numeric scalar. True causal effect used in simulation.}
#' }
#'
#' @examples
#' # Basic usage with default parameters
#' set.seed(123)
#' sim_data <- dgm_individual(
#'   m = 50,           # 50 genetic variants
#'   nx = 10000,       # 10,000 individuals in exposure study
#'   ny = 10000,       # 10,000 individuals in outcome study
#'   a_gamma = 0.1,    # IV effects between 0.1 and 0.3
#'   b_gamma = 0.3,
#'   a_f = 0.05,       # MAF between 5% and 45%
#'   b_f = 0.45,
#'   a_alpha = -0.1,   # Pleiotropy effects between -0.1 and 0.1
#'   b_alpha = 0.1,
#'   theta = 0.5       # True causal effect of 0.5
#' )
#'
#' # Check dimensions
#' length(sim_data$X)  # Should be 10000
#' length(sim_data$Y)  # Should be 10000
#' dim(sim_data$gx)    # Should be 50 x 10000
#' dim(sim_data$gy)    # Should be 50 x 10000
#'
#' # Scenario with no horizontal pleiotropy
#' sim_valid <- dgm_individual(
#'   m = 20, nx = 5000, ny = 5000,
#'   a_gamma = 0.2, b_gamma = 0.4,
#'   a_f = 0.1, b_f = 0.4,
#'   a_alpha = 0, b_alpha = 0,  # No pleiotropy
#'   theta = 0.3
#' )
#'
#' @export
#' @author [Your Name]
#' @seealso \code{\link{dgm_summary}} for summary-level GWAS data generation
dgm_individual <- function(m, nx, ny, a_gamma, b_gamma, a_f, b_f,
                           a_alpha, b_alpha, theta) {

  # Input validation
  if (!is.numeric(c(m, nx, ny)) || any(c(m, nx, ny) <= 0) || any(c(m, nx, ny) != floor(c(m, nx, ny)))) {
    stop("m, nx, and ny must be positive integers")
  }

  if (!is.numeric(c(a_gamma, b_gamma, a_f, b_f, a_alpha, b_alpha, theta))) {
    stop("All effect size and frequency parameters must be numeric")
  }

  if (a_gamma >= b_gamma) {
    stop("a_gamma must be less than b_gamma")
  }

  if (a_f >= b_f || a_f < 0 || b_f > 1) {
    stop("MAF bounds must satisfy 0 <= a_f < b_f <= 1")
  }

  if (a_alpha >= b_alpha) {
    stop("a_alpha must be less than b_alpha")
  }

  # Generate IV-to-exposure effects (instrument strength)
  gamma <- stats::runif(m, min = a_gamma, max = b_gamma)

  # Generate minor allele frequencies
  f <- stats::runif(m, min = a_f, max = b_f)

  # Generate genotype matrices [0, 1, 2] for m SNPs x n individuals
  # Each row represents a SNP, each column an individual
  gx <- replicate(nx, stats::rbinom(n = m, size = 2, prob = f))
  gy <- replicate(ny, stats::rbinom(n = m, size = 2, prob = f))

  # Standardize genotypes for each SNP (row-wise)
  # This ensures each SNP has mean 0 and variance 1
  gx <- t(apply(gx, 1, function(x) as.numeric(scale(x))))
  gy <- t(apply(gy, 1, function(x) as.numeric(scale(x))))

  # Handle constant genotypes (e.g., all individuals have same genotype)
  # Replace NaN from standardization with zeros
  gx[is.nan(gx)] <- 0
  gy[is.nan(gy)] <- 0

  # Generate horizontal pleiotropy effects
  alpha <- stats::runif(m, min = a_alpha, max = b_alpha)

  # Generate exposure phenotypes
  # X = gamma^T * G + error
  Xx <- as.numeric(gamma %*% gx) + stats::rnorm(nx, mean = 0, sd = 1)
  Xy <- as.numeric(gamma %*% gy) + stats::rnorm(ny, mean = 0, sd = 1)

  # Generate outcome phenotypes
  # Y = alpha^T * G + theta * X + error
  Y <- as.numeric(alpha %*% gy) + theta * Xy + stats::rnorm(ny, mean = 0, sd = 1)

  # Return comprehensive results
  return(list(
    X = Xx,
    Y = Y,
    gx = gx,
    gy = gy,
    gamma = gamma,
    alpha = alpha,
    f = f,
    theta = theta
  ))
}

