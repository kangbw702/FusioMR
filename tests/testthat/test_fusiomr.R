# Helper function to generate test data
generate_test_data <- function(n_instruments = 10, n_outcomes = 1) {
  causal_effect = 0.05
  set.seed(123)

  b_exp <- matrix(rnorm(n_instruments, 0, 0.1), ncol = 1)
  se_exp <- matrix(runif(n_instruments, 0.01, 0.05), ncol = 1)

  if (n_outcomes == 1) {
    b_out <- matrix(rnorm(n_instruments, causal_effect * b_exp, 0.1), ncol = 1)
    se_out <- matrix(runif(n_instruments, 0.01, 0.05), ncol = 1)
  } else if (n_outcomes == 2) {
    b_out1 <- rnorm(n_instruments, causal_effect * b_exp, 0.1)
    b_out2 <- rnorm(n_instruments, (causal_effect * 0.8) * b_exp, 0.1)
    b_out <- cbind(b_out1, b_out2)
    se_out1 <- runif(n_instruments, 0.01, 0.05)
    se_out2 <- runif(n_instruments, 0.01, 0.05)
    se_out <- cbind(se_out1, se_out2)
  }
  return(list(b_exp = b_exp, se_exp = se_exp, b_out = b_out, se_out = se_out))
}

# Input Validation Tests
test_that("Input validation works correctly", {
  data <- generate_test_data(n_instruments = 10, n_outcomes = 1)

  # Test numeric inputs
  expect_error(
    fusiomr("invalid", data$se_exp, data$b_out, data$se_out),
    "All effect estimates and standard errors must be numeric"
  )
  expect_error(
    fusiomr(data$b_exp, "invalid", data$b_out, data$se_out),
    "All effect estimates and standard errors must be numeric"
  )

  # Test dimension
  expect_error(
    fusiomr(data$b_exp[1:30, , drop = FALSE], data$se_exp, data$b_out, data$se_out),
    "All inputs must have the same number of rows"
  )

  # Test ncols
  b_exp_wrong <- cbind(data$b_exp, data$b_exp)
  expect_error(
    fusiomr(b_exp_wrong, data$se_exp, data$b_out, data$se_out),
    "IV-Exposure summary statistics \\(b_exp, se_exp\\) must have exactly 1 column"
  )

  se_out_wrong <- cbind(data$se_out, data$se_out)
  expect_error(
    fusiomr(data$b_exp, data$se_exp, data$b_out, se_out_wrong),
    "IV-Outcome summary statistics \\(b_out, se_out\\) must have the same number of columns"
  )

  b_out_wrong <- cbind(data$b_out, data$b_out, data$b_out)
  se_out_wrong <- cbind(data$se_out, data$se_out, data$se_out)
  expect_error(
    fusiomr(data$b_exp, data$se_exp, b_out_wrong, se_out_wrong),
    "IV-Outcome summary statistics \\(b_out, se_out\\) must have 1 or 2 columns"
  )

  # Test parameter ranges
  expect_error(
    fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, p_value_threshold = -0.1),
    "p_value_threshold must be between 0 and 1"
  )
  expect_error(
    fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, p_value_threshold = 1.1),
    "p_value_threshold must be between 0 and 1"
  )
  expect_error(
    fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, niter = -100),
    "niter must be positive integer"
  )
  expect_error(
    fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, burnin_prop = -0.1),
    "burnin_prop must be between 0 and 1"
  )
  expect_error(
    fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, burnin_prop = 1.1),
    "burnin_prop must be between 0 and 1"
  )

  # Test missing values
  b_exp_na <- data$b_exp
  b_exp_na[1, 1] <- NA
  expect_error(
    fusiomr(b_exp_na, data$se_exp, data$b_out, data$se_out),
    "Missing values are not allowed"
  )

  # Test standard errors
  se_exp_zero <- data$se_exp
  se_exp_zero[1, 1] <- 0
  expect_error(
    fusiomr(data$b_exp, se_exp_zero, data$b_out, data$se_out),
    "Standard errors must be positive"
  )
  se_out_negative <- data$se_out
  se_out_negative[1, 1] <- -0.01
  expect_error(
    fusiomr(data$b_exp, data$se_exp, data$b_out, se_out_negative),
    "Standard errors must be positive"
  )
})


# Test Model 1
test_that("Model 1 works correctly", {
  skip_if_not(exists("gibbs_seso_nohp"), message = "gibbs_seso_nohp function not available")
  skip_if_not(exists("get_sel_idx"), message = "get_sel_idx function not available")
  skip_if_not(exists("get_summary"), message = "get_summary function not available")
  skip_if_not(exists("print_summary"), message = "print_summary function not available")

  data <- generate_test_data(n_instruments = 10, n_outcomes = 1)

  result <- fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out)
  expect_type(result, "list")
  expect_named(result, c("est", "se", "pval"))
  expect_type(result$est, "double")
  expect_type(result$se, "double")
  expect_type(result$pval, "double")
  expect_length(result$est, 1)
  expect_length(result$se, 1)
  expect_length(result$pval, 1)
  expect_true(is.finite(result$est))
  expect_true(is.finite(result$se))
  expect_true(is.finite(result$pval))
  expect_true(result$se > 0)
  expect_true(result$pval >= 0 && result$pval <= 1)
})

# Test Model 2
test_that("Model 2 works correctly", {
  skip_if_not(exists("gibbs_seso_uhp_only"), message = "gibbs_seso_uhp_only function not available")
  skip_if_not(exists("get_sel_idx"), message = "get_sel_idx function not available")
  skip_if_not(exists("get_summary"), message = "get_summary function not available")
  skip_if_not(exists("print_summary"), message = "print_summary function not available")

  data <- generate_test_data(n_instruments = 10, n_outcomes = 1)

  result <- fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, CHP = TRUE)
  expect_type(result, "list")
  expect_named(result, c("est", "se", "pval"))
  expect_type(result$est, "double")
  expect_type(result$se, "double")
  expect_type(result$pval, "double")
  expect_length(result$est, 1)
  expect_length(result$se, 1)
  expect_length(result$pval, 1)
  expect_true(is.finite(result$est))
  expect_true(is.finite(result$se))
  expect_true(is.finite(result$pval))
  expect_true(result$se > 0)
  expect_true(result$pval >= 0 && result$pval <= 1)
})

# Test Model 3
test_that("Model 3 works correctly", {
  skip_if_not(exists("gibbs_semo_nohp"), message = "gibbs_semo_nohp function not available")
  skip_if_not(exists("get_summary"), message = "get_summary function not available")
  skip_if_not(exists("print_summary"), message = "print_summary function not available")

  data <- generate_test_data(n_instruments = 10, n_outcomes = 2)

  result <- fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out)
  expect_type(result, "list")
  expect_named(result, c("est", "se", "pval"))
  expect_type(result$est, "double")
  expect_type(result$se, "double")
  expect_type(result$pval, "double")
  expect_length(result$est, 2)
  expect_length(result$se, 2)
  expect_length(result$pval, 2)
  expect_true(all(is.finite(result$est)))
  expect_true(all(is.finite(result$se)))
  expect_true(all(is.finite(result$pval)))
  expect_true(all(result$se > 0))
  expect_true(all(result$pval >= 0 & result$pval <= 1))
})

# Test Model 4
test_that("Model 4 shows not implemented message", {
  data <- generate_test_data(n_instruments = 10, n_outcomes = 2)

  expect_output(
    fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out,
            CHP = TRUE, niter = 100, burnin_prop = 0.3),
    "has not implemented yet"
  )
})

# Test 7: Different Parameter Values
test_that("Different parameter values work correctly", {
  skip_if_not(exists("gibbs_seso_nohp"), message = "gibbs_seso_nohp function not available")
  skip_if_not(exists("get_sel_idx"), message = "get_sel_idx function not available")
  skip_if_not(exists("get_summary"), message = "get_summary function not available")
  skip_if_not(exists("print_summary"), message = "print_summary function not available")

  data <- generate_test_data(n_instruments = 10, n_outcomes = 1, causal_effect = 0.2)

  # Test with different p_value_threshold
  result1 <- fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out,
                     p_value_threshold = 0.1, niter = 100, burnin_prop = 0.4)
  result2 <- fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out,
                     p_value_threshold = 0.001, niter = 100, burnin_prop = 0.4)
  expect_type(result1, "list")
  expect_type(result2, "list")
  expect_true(is.finite(result1$est))
  expect_true(is.finite(result2$est))

  # Test with different niter
  result3 <- fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, niter = 100)
  expect_type(result3, "list")
  expect_true(is.finite(result3$est))

  # Test with different burnin_prop
  result4 <- fusiomr(data$b_exp, data$se_exp, data$b_out, data$se_out, niter = 500)
  expect_type(result4, "list")
  expect_true(is.finite(result4$est))
})

# Test Few Instruments
test_that("Warning is issued when few instruments are selected", {
  skip_if_not(exists("get_sel_idx"), message = "get_sel_idx function not available")

  n_instruments <- 10
  set.seed(456)
  b_exp <- matrix(rnorm(n_instruments, 0, 0.01), ncol = 1)
  se_exp <- matrix(runif(n_instruments, 0.05, 0.1), ncol = 1)
  b_out <- matrix(rnorm(n_instruments, 0, 0.1), ncol = 1)
  se_out <- matrix(runif(n_instruments, 0.05, 0.1), ncol = 1)
  expect_warning(
    fusiomr(b_exp, se_exp, b_out, se_out,
            p_value_threshold = 1e-5, niter = 100, burnin_prop = 0.3),
    "Less than 3 instruments selected"
  )
})


# Run all tests
cat("Running comprehensive tests for fusiomr function...\n")
test_dir(".", pattern = "test_fusiomr.R")
cat("All tests completed!\n")
