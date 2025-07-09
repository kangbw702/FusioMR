# R/helper_functions.R

#' Get indices for selected IVs
#'
#' @param b_exp Effect estimates for IV-exposure associations
#' @param se_exp Standard errors for IV-exposure associations
#' @param p_threshold P-value threshold for selection
#' @return Logical vector indicating selected IVs
#' @importFrom stats pnorm
#' @export
get_sel_idx <- function(b_exp, se_exp, p_threshold) {
  z_scores <- abs(b_exp / se_exp)
  p_values <- 2 * (1 - pnorm(z_scores))
  sel_idx <- p_values < p_threshold
  return(sel_idx)
}

#' Calculate empirical credible interval
#'
#' @param estimates Vector of post estimates
#' @param alpha Significance level (default: 0.05)
#' @return Vector of length 2 with lower and upper bounds
#' @importFrom stats quantile
#' @export
get_empirical_ci <- function(estimates, alpha = 0.05) {
  return(stats::quantile(estimates, c(alpha/2, 1-alpha/2)))
}

#' Calculate normal approximation confidence interval
#'
#' @param mean Estimated mean
#' @param se Standard error
#' @param alpha Significance level (default: 0.05)
#' @return Vector of length 2 with lower and upper bounds
#' @importFrom stats qnorm
#' @export
get_normal_ci <- function(mean, se, alpha = 0.05) {
  z_val <- stats::qnorm(1 - alpha/2)
  return(c(mean - z_val * se, mean + z_val * se))
}

#' Get summary statistics from posterior estimates
#'
#' @param estimates Vector of posterior MCMC estimations
#' @return List with beta_est, beta_se, beta_pval, ci_emp, and ci_normal
#' @importFrom stats pnorm sd
#' @export
get_summary <- function(estimates) {
  est_mean <- mean(estimates)
  est_se <- stats::sd(estimates)
  est_pval <- 2 * (1 - stats::pnorm(abs(est_mean / est_se)))
  ci_emp <- get_empirical_ci(estimates)
  ci_normal <- get_normal_ci(est_mean, est_se)

  res_summary <- list(
    beta_est = est_mean,
    beta_se = est_se,
    beta_pval = est_pval,
    ci_emp = ci_emp,
    ci_normal = ci_normal
  )

  return(res_summary)
}

#' Print summary results
#'
#' @param res_summary Result summary object from get_summary()
#' @return NULL
#' @export
print_summary <- function(res_summary) {
  cat(sprintf("Estimated Causal Effect (Beta): %.4f\n", res_summary$beta_est))
  cat(sprintf("Standard Error: %.4f\n", res_summary$beta_se))
  cat(sprintf("P-value: %.4f\n", res_summary$beta_pval))
  cat(sprintf("95%% Empirical CI: [%.4f, %.4f]\n", res_summary$ci_emp[1], res_summary$ci_emp[2]))
  cat(sprintf("95%% Normal CI: [%.4f, %.4f]\n", res_summary$ci_normal[1], res_summary$ci_normal[2]))
}
