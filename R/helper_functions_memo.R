#' A helper function to drop burnin proportion
#'
#' Removes burn-in samples from MCMC results to obtain post-burnin posterior samples
#' for subsequent analysis.
#'
#' @param res A list containing MCMC results with components:
#'   \itemize{
#'     \item \code{pst_tk}: Matrix of posterior samples for joint probabilities.
#'     \item \code{beta_1_tk}: Vector of posterior samples for causal effect of outcome1
#'     \item \code{beta_2_tk}: Vector of posterior samples for causal effect of outcome2
#'     \item \code{eta_1_tk}: Vector of posterior samples for eta1
#'     \item \code{eta_2_tk}: Vector of posterior samples for eta2
#'     \item \code{alpha_1_tk}: Vector of posterior samples for alpha1
#'     \item \code{alpha_2_tk}: Vector of posterior samples for alpha2
#'   }
#' @param niter An integer represents the total number of MCMC iterations.
#' @param burnin_prop A number between 0 and 1 representing burn-in proportion.
#' @return A list with the same structure as \code{res}, containing only post-burnin samples.
#'
#' @export
get_post_burnin_res <- function(res, niter, burnin_prop) {
  burnin <- floor(niter * burnin_prop)
  post_res <- list(pst_tk = res$pst_tk[(burnin+1):niter, , drop = FALSE],
                   beta_1_tk = res$beta_1_tk[(burnin+1):niter],
                   beta_2_tk = res$beta_2_tk[(burnin+1):niter],
                   eta_1_tk = res$eta_1_tk[(burnin+1):niter],
                   eta_2_tk = res$eta_2_tk[(burnin+1):niter],
                   alpha_1_tk = res$alpha_1_tk[(burnin+1):niter],
                   alpha_2_tk = res$alpha_2_tk[(burnin+1):niter])
  return(post_res)
}

#' A helper function to post-processing for multiple exposure, multiple outcome joint model.
#'
#' This function performs label-switching correction for memo model by
#' adjusting parameter estimates based on posterior probabilities.
#'
#' @param res Object containing posterior samples from MCMC estimation. Must include
#'   fields: \code{pst_tk}, \code{beta_1_tk}, \code{alpha_1_tk}, \code{beta_2_tk},
#'   and \code{alpha_2_tk}.
#' @param eta_true_1 Numeric. True value of first eta parameter indicating if IVs are valid
#' @param eta_true_2 Numeric. True value of second eta parameter indicating if IVs are valid
#' @param niter Integer. Total number of MCMC iterations.
#' @param burnin_prop Numeric. Proportion of iterations to discard as burn-in
#'   (between 0 and 1).
#'
#' @return A list containing statistical summaries for both corrected parameters:
#' \describe{
#'   \item{beta_est1}{Posterior mean for causal effect of outcome1}
#'   \item{beta_se1}{Posterior standard deviation for outcome1}
#'   \item{beta_est2}{Posterior mean for causal effect of outcome2}
#'   \item{beta_se2}{Posterior standard deviation for outcome2}
#'   \item{beta_pval1}{Two-sided p-value for beta1 using normal approximation}
#'   \item{beta_pval2}{Two-sided p-value for beta2 using normal approximation}
#'   \item{ci_emp1}{95\% empirical credible interval for beta1}
#'   \item{ci_emp2}{95\% empirical credible interval for beta2}
#' }
#'
#' @details
#' The function performs label-switching correction using the following algorithm:
#'
#' 1. Extracts post-burn-in posterior samples using \code{get_post_burnin_res()}
#' 2. Calculates marginal probabilities:
#'    - qq1 = mean(p10_post + p11_post) for first dimension
#'    - qq2 = mean(p01_post + p11_post) for second dimension
#' 3. Applies correction based on 0.5 threshold:
#'    - If qq1 > 0.5: Uses beta_1_tk + alpha_1_tk (label-switched)
#'    - If qq1 ≤ 0.5: Uses beta_1_tk (original labels)
#'    - Same logic applies for second dimension
#' 4. Computes posterior statistics and p-values using normal approximation
#'
#' The label-switching occurs when the marginal probability exceeds 0.5, indicating
#' that the labels may have been flipped during MCMC sampling.
#'
#' @export
label_flip_joint <- function(res, eta_true_1, eta_true_2, niter, burnin_prop) {
  post_res = get_post_burnin_res(res, niter, burnin_prop)
  p00_post = post_res$pst_tk[, 1]
  p01_post = post_res$pst_tk[, 2]
  p10_post = post_res$pst_tk[, 3]
  p11_post = post_res$pst_tk[, 4]
  qq1 = mean(p10_post + p11_post)
  qq2 = mean(p01_post + p11_post)
  if (qq1 > 0.5) {
    b1_mean = mean(post_res$beta_1_tk + post_res$alpha_1_tk)
    b1_sd = sd(post_res$beta_1_tk + post_res$alpha_1_tk)
    bci1 = quantile(post_res$beta_1_tk + post_res$alpha_1_tk, probs=c(0.025,0.975), na.rm=T)
  }
  if (qq1 <=0.5) {
    b1_mean = mean(post_res$beta_1_tk)
    b1_sd = sd(post_res$beta_1_tk)
    bci1 = quantile(post_res$beta_1_tk, probs=c(0.025,0.975), na.rm=T)
  }
  if (qq2 > 0.5) {
    b2_mean = mean(post_res$beta_2_tk + post_res$alpha_2_tk)
    b2_sd = sd(post_res$beta_2_tk + post_res$alpha_2_tk)
    bci2 = quantile(post_res$beta_2_tk + post_res$alpha_2_tk, probs=c(0.025,0.975), na.rm=T)
  }
  if (qq2 <=0.5) {
    b2_mean = mean(post_res$beta_2_tk)
    b2_sd = sd(post_res$beta_2_tk)
    bci2 = quantile(post_res$beta_2_tk, probs=c(0.025,0.975), na.rm=T)
  }
  est_pval1 <- 2 * (1 - stats::pnorm(abs(b1_mean / b1_sd)))
  est_pval2 <- 2 * (1 - stats::pnorm(abs(b2_mean / b2_sd)))
  return(list(beta_est1 = b1_mean,
              beta_se1 = b1_sd,
              beta_est2 = b2_mean,
              beta_se2 = b2_sd,
              beta_pval1 = est_pval1,
              beta_pval2 = est_pval2,
              ci_emp1 = bci1,
              ci_emp2 = bci2))
}

#' Print Summary of Causal Effect Results
#'
#' Prints a formatted summary of causal effect estimation results for outcomes.
#'
#' @param post_res A list containing summary results from \code{\link{label_flip_joint}} with components:
#'   \itemize{
#'     \item \code{beta_est1}, \code{beta_est2}: Point estimates
#'     \item \code{beta_se1}, \code{beta_se2}: Standard errors
#'     \item \code{beta_pval1}, \code{beta_pval2}: P-values
#'     \item \code{ci_emp1}, \code{ci_emp2}: 95% credible intervals
#'   }
#'
#' @return No return value. Prints formatted results to console.
#'
#' @details Provides a clean, formatted output displaying:
#'   \itemize{
#'     \item Estimated causal effects (beta coefficients)
#'     \item Standard errors
#'     \item Two-sided p-values
#'     \item 95% empirical credible intervals
#'   }
#'   Results are displayed separately for each outcome with clear section headers.
#'
#' @export
print_summary_memo <- function(post_res) {
  cat("\n--- RESULTS for Outcome1---\n")
  cat(sprintf("Estimated Causal Effect (Beta): %.4f\n", post_res$beta_est1))
  cat(sprintf("Standard Error: %.4f\n", post_res$beta_se1))
  cat(sprintf("P-value: %.4f\n", post_res$beta_pval1))
  cat(sprintf("95%% Empirical CI: [%.4f, %.4f]\n", post_res$ci_emp1[1], post_res$ci_emp1[2]))
  cat("\n--- RESULTS for Outcome2---\n")
  cat(sprintf("Estimated Causal Effect (Beta): %.4f\n", post_res$beta_est2))
  cat(sprintf("Standard Error: %.4f\n", post_res$beta_se2))
  cat(sprintf("P-value: %.4f\n", post_res$beta_pval2))
  cat(sprintf("95%% Empirical CI: [%.4f, %.4f]\n", post_res$ci_emp2[1], post_res$ci_emp2[2]))
}
