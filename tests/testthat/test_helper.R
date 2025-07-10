# tests/testthat/test_helper_functions.R
library(testthat)
library(FusioMR)

test_that("get_sel_idx works correctly", {
  # Create test data
  b_exp <- matrix(c(0.1, 0.05, 0.2, 0.01), ncol = 1)
  se_exp <- matrix(c(0.02, 0.03, 0.01, 0.05), ncol = 1)

  result1 <- get_sel_idx(b_exp, se_exp, 0.05)
  result2 <- get_sel_idx(b_exp, se_exp, 0.001)
  expect_type(result1, "logical")
  expect_type(result2, "logical")
  expect_length(result1, 4)
  expect_length(result2, 4)
  expect_true(sum(result2) <= sum(result1))
})

test_that("get_empirical_ci works correctly", {
  # Create test samples
  samples <- rnorm(1000, mean = 0.5, sd = 0.1)

  ci_default <- get_empirical_ci(samples)
  expect_length(ci_default, 2)
  expect_true(ci_default[1] < ci_default[2])
  ci_90 <- get_empirical_ci(samples, alpha = 0.1)
  expect_length(ci_90, 2)
  expect_true(ci_90[1] < ci_90[2])
  expect_true((ci_90[2] - ci_90[1]) < (ci_default[2] - ci_default[1]))
})

test_that("get_normal_ci works correctly", {
  # Test normal CI calculation
  mean_val <- 0.5
  se_val <- 0.1

  ci_default <- get_normal_ci(mean_val, se_val)
  expect_length(ci_default, 2)
  expect_true(ci_default[1] < ci_default[2])
  ci_90 <- get_normal_ci(mean_val, se_val, alpha = 0.1)
  expect_length(ci_90, 2)
  expect_true((ci_90[2] - ci_90[1]) < (ci_default[2] - ci_default[1]))
})

test_that("get_summary works correctly", {
  # Create test samples
  set.seed(123)
  samples <- rnorm(1000, mean = 0.5, sd = 0.1)
  result <- get_summary(samples)
  expect_type(result, "list")
  expected_names <- c("beta_est", "beta_se", "beta_pval", "ci_emp", "ci_normal")
  expect_named(result, expected_names)
  expect_true(abs(result$beta_est - 0.5) < 0.05)
  expect_true(abs(result$beta_se - 0.1) < 0.02)
  expect_true(result$beta_pval >= 0 && result$beta_pval <= 1)

  expect_length(result$ci_emp, 2)
  expect_length(result$ci_normal, 2)
  expect_true(result$ci_emp[1] < result$ci_emp[2])
  expect_true(result$ci_normal[1] < result$ci_normal[2])
})

test_that("print_summary works correctly", {
  # Create test object
  summary_obj <- list(
    beta_est = 0.5,
    beta_se = 0.1,
    beta_pval = 0.001,
    ci_emp = c(0.1, 0.2),
    ci_normal = c(0.1, 0.2)
  )

  # Test that it prints expected content
  expect_output(print_summary(summary_obj), "Estimated Causal Effect")
  expect_output(print_summary(summary_obj), "0.5000")
  expect_output(print_summary(summary_obj), "Standard Error")
  expect_output(print_summary(summary_obj), "P-value")
  expect_output(print_summary(summary_obj), "Empirical CI")
  expect_output(print_summary(summary_obj), "Normal CI")
})


test_that("small samples work correctly", {
  # Test with very small samples
  small_samples <- rnorm(10, 0.5, 0.1)
  result_small <- get_summary(small_samples)
  expect_type(result_small, "list")
})
