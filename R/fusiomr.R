#' FusioMR: Main Function for MR analysis
#'
#' @description
#' A Flexible, Unified and verSatile Mendelian Randomization framework
#'
#' @param summary_stats_raw A data frame with 4 columns: b_exp, se_exp, b_out, se_out
#' @param type A string specifying model type
#' @param p_value_threshold P-value threshold for IV selection (default: 1e-3)
#' @param niter Number of MCMC iterations (default: 20000)
#' @param burnin_prop Proportion of iterations to discard as burn-in (default: 0.5)
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_estimate}{Posterior mean of causal effect}
#'   \item{beta_se}{Posterior standard error}
#'   \item{pval}{P-value for causal effect}
#'   \item{beta_ci_emp}{95 percent empirical credible interval}
#'   \item{beta_ci_normal}{95 percent normal-approximation confidence interval}
#'   \item{beta_samples}{MCMC samples (post burn-in)}
#'   \item{type}{Model type used}
#'   \item{n_ivs}{Number of instrumental variables selected}
#'   \item{selected_stats}{Selected summary statistics}
#'   \item{niter}{Number of MCMC iterations}
#'   \item{burnin_prop}{Burn-in proportion used}
#' }
#'
#' @export
#'
fusiomr <- function(summary_stats_raw,
                    type,
                    p_value_threshold = 1e-3,
                    niter = 20000,
                    burnin_prop = 0.5) {

  # Input validation
  if (!is.data.frame(summary_stats_raw) || ncol(summary_stats_raw) != 4) {
    stop("Wrong input!! \n summary_stats_raw must be a data frame with 4 columns: b_exp, se_exp, b_out, se_out")
  }
  required_cols <- c("b_exp", "se_exp", "b_out", "se_out")
  if (!all(required_cols %in% names(summary_stats_raw))) {
    stop("Wrong input!! \n summary_stats_raw must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # types
  if (!type %in% c("seso_uhp_only", "seso_nohp", "multi_uhp_only", "multi_nohp")) {
    stop("Wrong input!! \n Only 4 option for type input: seso_uhp_only, seso_nohp, multi_uhp_only, multi_nohp")
  }

  # niter
  if (niter <= 0) {
    stop("Wrong input!! \n niter must be positive integer.")
  }

  # burnin_prop
  if (burnin_prop < 0 || burnin_prop >= 1) {
    stop("Wrong input!! \n burnin_prop must be between 0 and 1")
  }

  cat("=== FusioMR Analysis Start ===\n \n")
  cat("Model type:", type, "\n")
  cat("P-value threshold:", p_value_threshold, "\n")
  cat("MCMC iterations:", niter, "\n\n")

  # Step 1: IV Selection
  cat("Step 1: Instrumental Variable Selection...\n")
  summary_stats_selected_result <- summary_stats_selected(summary_stats_raw, p_threshold = p_value_threshold)

  Gamma_hat <- summary_stats_selected_result$b_out
  gamma_hat <- summary_stats_selected_result$b_exp
  s2_hat_Gamma <- summary_stats_selected_result$se_out^2
  s2_hat_gamma <- summary_stats_selected_result$se_exp^2

  # Step 2: Gibbs Sampling
  cat("\nStep 2: Gibbs Sampling...\n")
  if (type == "seso_uhp_only") {
    gibbs_beta_est <- gibbs_seso_uhp_only(niter, Gamma_hat, gamma_hat, s2_hat_Gamma, s2_hat_gamma)
  } else if (type == "seso_nohp") {
    gibbs_beta_est <- gibbs_seso_nohp(niter, Gamma_hat, gamma_hat, s2_hat_Gamma)
  } else if (type == "multi_uhp_only") {
    stop("multi_uhp_only gibbs has not been implemented yet")
  } else if (type == "multi_nohp") {
    stop("multi_nohp gibbs has not been implemented yet")
  }

  # Process MCMC results
  burnin_n <- floor(niter * burnin_prop)
  beta_post_burnin <- gibbs_beta_est[(burnin_n + 1):niter]

  # Causal effect estimate
  beta_estimate <- mean(beta_post_burnin)

  # Standard Error
  beta_se <- stats::sd(beta_post_burnin)

  # P-value calculation
  beta_pvalue <- 2 * exp(stats::pnorm(abs(beta_estimate) / beta_se, lower.tail = FALSE, log.p = TRUE))

  # Empirical credible interval (quantile-based)
  beta_ci_lower_emp <- stats::quantile(beta_post_burnin, 0.025)
  beta_ci_upper_emp <- stats::quantile(beta_post_burnin, 0.975)

  # Normal approximation confidence interval (critical value approach)
  beta_ci_lower_normal <- beta_estimate - 1.96 * beta_se
  beta_ci_upper_normal <- beta_estimate + 1.96 * beta_se

  # Store all the results
  result <- list(
    # Main results
    beta_estimate = beta_estimate,
    beta_se = beta_se,
    beta_pvalue = beta_pvalue,
    beta_ci_emp = c(beta_ci_lower_emp, beta_ci_upper_emp),
    beta_ci_normal = c(beta_ci_lower_normal, beta_ci_upper_normal),

    # MCMC samples
    beta_samples = beta_post_burnin,

    # Method information
    type = type,
    n_ivs = nrow(summary_stats_selected_result),

    # Selected data
    selected_stats = summary_stats_selected_result,

    # MCMC settings
    niter = niter,
    burnin_prop = burnin_prop
  )

  # Print summary
  cat("\nStep 3: Results Summary...\n")
  cat("Causal effect estimate (beta):", round(beta_estimate, 4), "\n")
  cat("Standard error:", round(beta_se, 4), "\n")
  cat("P-value:", formatC(beta_pvalue, format = "e", digits = 3), "\n")
  cat("95% Empirical credible interval: [", round(beta_ci_lower_emp, 4), ",", round(beta_ci_upper_emp, 4), "]\n")
  cat("95% Normal approximation CI: [", round(beta_ci_lower_normal, 4), ",", round(beta_ci_upper_normal, 4), "]\n")
  cat("Number of IVs used:", nrow(summary_stats_selected_result), "\n\n")
  cat("=== FusioMR Analysis Completed ===\n\n")

  class(result) <- "fusiomr"
  return(result)
}
