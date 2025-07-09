# tests/testthat/test-fusiomr.R
library(testthat)
library(FusioMR)

test_that("fusiomr basic functionality works", {
  # create test data
  set.seed(123)
  n_ivs <- 100  # Very small for fast testing
  b_exp <- rnorm(n_ivs, 0, 0.1)
  se_exp <- rep(0.01, n_ivs)  # Fixed SE for simplicity
  true_beta <- 0.3
  b_out <- true_beta * b_exp + rnorm(n_ivs, 0, 0.01)
  se_out <- rep(0.01, n_ivs)

  # test seso
  result <- fusiomr(b_exp, se_exp, b_out, se_out,
                    CHP = FALSE,
                    niter = 100,
                    p_value_threshold = 0.01)

  expect_type(result, "list")
  expect_true(all(c("est", "se", "pval") %in% names(result)))
  expect_length(result$est, 1)
  expect_length(result$se, 1)
  expect_length(result$pval, 1)
  expect_true(is.finite(result$est))
  expect_true(is.finite(result$se))
  expect_true(result$se > 0)
  expect_true(result$pval >= 0 && result$pval <= 1)
})

test_that("input validation works", {
  # Quick validation tests
  expect_error(fusiomr("not_numeric", 1, 1, 1))
  expect_error(fusiomr(c(1,2), c(1), c(1,2), c(1,2)))
  expect_error(fusiomr(c(1,2), c(1,2), c(1,2), c(-1,2)))
})

test_that("multiple outcomes work", {
  # create test data
  set.seed(456)
  n_ivs <- 100
  b_exp <- rnorm(n_ivs, 0, 0.1)
  se_exp <- rep(0.01, n_ivs)
  b_out <- cbind(0.2 * b_exp + rnorm(n_ivs, 0, 0.01),
                 0.4 * b_exp + rnorm(n_ivs, 0, 0.01))
  se_out <- matrix(0.01, n_ivs, 2)

  result <- fusiomr(b_exp, se_exp, b_out, se_out,
                    CHP = FALSE, niter = 50,
                    p_value_threshold = 0.01)
  expect_length(result$est, 2)
  expect_length(result$se, 2)
  expect_length(result$pval, 2)
})
