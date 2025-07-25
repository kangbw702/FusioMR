#' FusioMR: Main Function for Mendelian Randomization Analysis
#'
#' @description
#' A Flexible, Unified and verSatile Mendelian Randomization framework.
#' Supports both single and multiple outcome analyses with optional horizontal pleiotropy.
#'
#' @param b_exp Numeric matrix or vector of effect estimates for instrumental variables
#'   on exposure. Must have exactly 1 column and same number of rows as other inputs.
#' @param se_exp Numeric matrix or vector of standard errors for IV-exposure effect
#'   estimates. Must have exactly 1 column, same dimensions as b_exp, and all positive values.
#' @param b_out Numeric matrix or vector of effect estimates for instrumental variables
#'   on outcome(s). Must have 1 column (single outcome) or 2 columns (multiple outcomes),
#'   with same number of rows as b_exp.
#' @param se_out Numeric matrix or vector of standard errors for IV-outcome effect
#'   estimates. Must have same dimensions as b_out and all positive values.
#' @param CHP Logical indicating whether to include horizontal pleiotropy correction
#'   (default: FALSE).
#' @param p_value_threshold Numeric value between 0 and 1 for IV selection threshold
#'   (default: 1e-3). Used to filter instruments based on IV-exposure association strength.
#' @param niter Positive integer specifying number of Gibbs sampling iterations
#'   (default: 20000).
#' @param burnin_prop Numeric value between 0 and 1 specifying proportion of iterations
#'   to discard as burn-in (default: 0.5).
#'
#' \itemize{
#'   \item Model 1: Single Exposure Single Outcome (No Horizontal Pleiotropy)
#'   \item Model 2: Single Exposure Single Outcome (With Horizontal Pleiotropy)
#'   \item Model 3: Single Exposure Multiple Outcomes (No Horizontal Pleiotropy)
#'   \item Model 4: Single Exposure Multiple Outcomes (With Horizontal Pleiotropy)
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{est}{Numeric vector of posterior mean(s) of causal effect. Single value for
#'     single outcome, vector of length 2 for multiple outcomes.}
#'   \item{se}{Numeric vector of posterior standard error(s) of causal effect. Single value
#'     for single outcome, vector of length 2 for multiple outcomes.}
#'   \item{pval}{Numeric vector of p-value(s) for testing null hypothesis of no causal effect.
#'     Single value for single outcome, vector of length 2 for multiple outcomes.}
#' }
#'
#' @export
#'
fusiomr <- function(b_exp,
                    se_exp,
                    b_out,
                    se_out,
                    CHP = FALSE,
                    p_value_threshold = 1e-3,
                    niter = 20000,
                    burnin_prop = 0.5) {

  # Input validation
  cat("\n=== FusioMR: Mendelian Randomization Analysis ===\n")
  cat("Validating inputs...\n")

  # Check if inputs are numeric
  if (!is.numeric(b_exp) || !is.numeric(se_exp) || !is.numeric(b_out) || !is.numeric(se_out)) {
    stop("Wrong input!! \n All effect estimates and standard errors must be numeric.")
  }

  if (is.vector(b_exp)) b_exp <- as.matrix(b_exp)
  if (is.vector(se_exp)) se_exp <- as.matrix(se_exp)
  if (is.vector(b_out)) b_out <- as.matrix(b_out)
  if (is.vector(se_out)) se_out <- as.matrix(se_out)

  # Check dimensions
  if (nrow(b_exp) != nrow(se_exp) || nrow(b_exp) != nrow(b_out) || nrow(b_exp) != nrow(se_out)) {
    stop("Wrong input!! \n All inputs must have the same number of rows.")
  }

  # Check that exposure data has only 1 column
  if (ncol(b_exp) != ncol(se_exp)) {
    stop("Wrong input!! \n IV-Exposure summary statistics (b_exp, se_exp) must have the same number of columns.")
  }

  # Check that outcome data has 1 or 2 columns
  if (ncol(b_out) != ncol(se_out)) {
    stop("Wrong input!! \n IV-Outcome summary statistics (b_out, se_out) must have the same number of columns.")
  }

  if (!(ncol(b_exp) %in% c(1, 2))) {
    stop("Wrong input!! \n IV-Exposure summary statistics (b_out, se_out) must have 1 or 2 columns.")
  }

  if (!(ncol(b_out) %in% c(1, 2))) {
    stop("Wrong input!! \n IV-Outcome summary statistics (b_out, se_out) must have 1 or 2 columns.")
  }

  # Check parameter ranges
  if (p_value_threshold < 0 || p_value_threshold >= 1) {
    stop("Wrong input!! \n p_value_threshold must be between 0 and 1")
  }

  if (niter <= 0) {
    stop("Wrong input!! \n niter must be positive integer.")
  }

  if (burnin_prop < 0 || burnin_prop >= 1) {
    stop("Wrong input!! \n burnin_prop must be between 0 and 1")
  }

  # Check for missing values
  if (any(is.na(b_exp)) || any(is.na(se_exp)) || any(is.na(b_out)) || any(is.na(se_out))) {
    stop("Wrong input!! \n Missing values are not allowed.")
  }

  # Check for standard errors
  if (any(se_exp <= 0) || any(se_out <= 0)) {
    stop("Wrong input!! \n Standard errors must be positive.")
  }

  cat("Input validation completed.\n")


  n_outcomes <- ncol(b_out)
  n_exposure <- ncol(b_exp)
  n_ivs <- nrow(b_exp)

  # Model Implementation
  results <- list()
  if (n_exposure == 1) { # Single Exposure
    if (n_outcomes == 1) { # Single outcome models
      # IV selection
      cat("Performing instrumental variable selection...\n")
      sel_ivs_idx <- get_sel_idx(b_exp, se_exp, p_value_threshold)
      n_sel <- sum(sel_ivs_idx)
      cat(sprintf("Selected %d out of %d instruments (p < %g)\n", n_sel, n_ivs, p_value_threshold))

      if (n_sel < 3) {
        warning("Less than 3 instruments selected. Results may be unreliable.")
      }

      # Select summary statistics
      b_exp_sel <- b_exp[sel_ivs_idx, , drop = FALSE]
      se_exp_sel <- se_exp[sel_ivs_idx, , drop = FALSE]
      b_out_sel <- b_out[sel_ivs_idx, , drop = FALSE]
      se_out_sel <- se_out[sel_ivs_idx, , drop = FALSE]
      if (!CHP) {
        # Model 1: seso without CHP
        cat("\n--- Running Model: Single Exposure Single Outcome (No Horizontal Pleiotropy) ---\n")

        # Gibbs
        cat("Running Gibbs sampling...\n")
        cat(sprintf("Iterations: %d, Burn-in: %d\n", niter, floor(niter * burnin_prop)))
        beta_est <- gibbs_seso_nohp(niter, b_out_sel, b_exp_sel, se_out_sel^2, se_exp_sel^2)
        burnin <- floor(niter * burnin_prop)
        beta_est <- beta_est[(burnin + 1):niter]

        # Get results statistics
        res_summary <- get_summary(beta_est)
        # Print results summary
        cat("\n--- RESULTS---\n")
        print_summary(res_summary)

        # Return
        results <- list(
          est = res_summary$beta_est,
          se = res_summary$beta_se,
          pval = res_summary$beta_pval
        )
      } else {
        # Model 2: seso with CHP
        cat("\n--- Running Model: Single Exposure Single Outcome (With Horizontal Pleiotropy) ---\n")

        # Gibbs
        cat("Running Gibbs sampling...\n")
        cat(sprintf("Iterations: %d, Burn-in: %d\n", niter, floor(niter * burnin_prop)))
        beta_est <- gibbs_seso_uhp_only(niter, b_out_sel, b_exp_sel, se_out_sel^2, se_exp_sel^2)
        burnin <- floor(niter * burnin_prop)
        beta_est <- beta_est[(burnin + 1):niter]

        # Get results statistics
        res_summary <- get_summary(beta_est)
        # Print results summary
        cat("\n--- RESULTS---\n")
        print_summary(res_summary)

        # Return
        results <- list(
          est = res_summary$beta_est,
          se = res_summary$beta_se,
          pval = res_summary$beta_pval
        )
      }
    } else { # Multiple outcome models (n_outcomes == 2)
      if (!CHP) {
        # Model 3: semo without CHP
        cat("\n--- Running Model: Single Exposure Multiple Outcomes (No Horizontal Pleiotropy) ---\n")

        # Gibbs
        cat("Running Gibbs sampling...\n")
        cat(sprintf("Iterations: %d, Burn-in: %d\n", niter, floor(niter * burnin_prop)))

        beta_est <- gibbs_semo_nohp(niter, b_out[, 1], b_out[, 2], b_exp,
                                    se_out[, 1], se_out[, 2], se_exp)
        burnin <- floor(niter * burnin_prop)
        beta_est1 <- beta_est$beta_1[(burnin + 1):niter]
        beta_est2 <- beta_est$beta_2[(burnin + 1):niter]

        # Get results statistics
        res_summary1 <- get_summary(beta_est1)
        res_summary2 <- get_summary(beta_est2)
        # Print results summary
        cat("\n--- RESULTS for Outcome1---\n")
        print_summary(res_summary1)
        cat("\n--- RESULTS for Outcome2---\n")
        print_summary(res_summary2)

        # Return
        results <- list(
          est = c(res_summary1$beta_est, res_summary2$beta_est),
          se = c(res_summary1$beta_se, res_summary2$beta_se),
          pval = c(res_summary1$beta_pval, res_summary2$beta_pval)
        )
      } else {
        # Model 4: semo with uhp
        cat("\n--- Running Model: Single Exposure Multiple Outcomes (With Horizontal Pleiotropy)---\n")
        cat("has not implemented yet")
      }
    }
  } else { # Multiple Exposure
    if (n_outcomes == 2) { # Multiple Outcomes
      # Model 5: memo uhp+chp
      cat("\n--- Running Model: Multiple Exposure Multiple Outcomes---\n")
      # Gibbs
      cat("Running Gibbs sampling...\n")
      cat(sprintf("Iterations: %d, Burn-in: %d\n", niter, floor(niter * burnin_prop)))
      burnin <- floor(niter * burnin_prop)
      res <- gibbs_memo_joint(niter, b_out[, 1], b_out[, 2], b_exp[, 1], b_exp[, 2],
                              se_out[, 1], se_out[, 2], se_exp[, 1], se_exp[, 2])

      # Debug the results
      # debug_gibbs_results(res, niter, burnin_prop)

      eta_true_1 <- rep(0, n_ivs)
      eta_true_2 <- rep(0, n_ivs)

      # Use debug version of label_flip_joint
      post_res <- label_flip_joint(res, eta_true_1, eta_true_2, niter, burnin_prop)
      print_summary_memo(post_res)
      results <- list(
        est = c(post_res$beta_est1, post_res$beta_est2),
        se = c(post_res$beta_se1, post_res$beta_se2),
        pval = c(post_res$beta_pval1, post_res$beta_pval2)
      )
    }
  }
  cat("=== Analysis Complete ===\n \n")
  return(results)
}

