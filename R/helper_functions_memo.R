get_post_burnin_res <- function(res, niter, burnin_prop) {
  burnin <- floor(niter * burnin_prop)
  post_res <- list(pst_tk = res$pst_tk[(burnin+1):niter],
                   beta_1_tk = res$beta_1_tk[(burnin+1):niter],
                   beta_2_tk = res$beta_2_tk[(burnin+1):niter],
                   eta_1_tk = res$eta_1_tk[(burnin+1):niter],
                   eta_2_tk = res$eta_2_tk[(burnin+1):niter],
                   alpha_1_tk = res$alpha_1_tk[(burnin+1):niter],
                   alpha_2_tk = res$alpha_2_tk[(burnin+1):niter])

  return(post_res)
}

label_flip_joint <- function(res, eta_true_1, eta_true_2) {
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
  est_pval2 <- 2 * (1 - stats::pnorm(abs(b2_mean / b2_se)))

  return(list(beta_est1 = b1_mean,
              beta_se1 = b1_sd,
              beta_est2 = b2_mean,
              beta_se2 = b2_sd,
              beta_pval1 = est_pval1,
              beta_pval2 = est_pval2,
              ci_emp1 = bci1,
              ci_emp2 = bci2))
}


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
