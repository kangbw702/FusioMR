test_that("fusiomr basic functionality works", {
  # Create sample data
  set.seed(123)
  n_snps <- 100
  test_data <- data.frame(
    b_exp = rnorm(n_snps, 0, 0.1),
    se_exp = runif(n_snps, 0.01, 0.05),
    b_out = rnorm(n_snps, 0, 0.1),
    se_out = runif(n_snps, 0.01, 0.05)
  )

  # Test the function
  result <- fusiomr(test_data, type = "seso_nohp", niter = 1000)

  expect_s3_class(result, "fusiomr")
  expect_true(is.numeric(result$beta_estimate))
  expect_true(length(result$beta_ci) == 2)
})

test_that("fusiomr works with valid input", {
  # Create sample summary statistics
  set.seed(42)
  summary_stats <- data.frame(
    b_exp = rnorm(50, mean = 0, sd = 0.1),
    se_exp = runif(50, min = 0.01, max = 0.05),
    b_out = rnorm(50, mean = 0, sd = 0.1),
    se_out = runif(50, min = 0.01, max = 0.05)
  )

  # Test with seso_uhp_only
  result <- fusiomr(summary_stats_raw = summary_stats,
                    type = "seso_uhp_only",
                    p_value_threshold = 0.05,
                    niter = 1000,
                    burnin_prop = 0.5)

  expect_s3_class(result, "fusiomr")
  expect_true(is.numeric(result$beta_estimate))
  expect_true(is.numeric(result$beta_se))
  expect_equal(length(result$beta_ci), 2)
  expect_equal(result$type, "seso_uhp_only")
})

test_that("summary_stats_selected works correctly", {
  # Create sample data
  set.seed(456)
  summary_stats <- data.frame(
    b_exp = c(0.5, 0.1, 0.01, 0.8),  # Mix of strong and weak effects
    se_exp = c(0.1, 0.05, 0.02, 0.15),
    b_out = rnorm(4),
    se_out = runif(4, 0.01, 0.05)
  )

  # Test IV selection
  result <- summary_stats_selected(summary_stats, p_threshold = 0.05)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 4)
  expect_true(nrow(result) <= nrow(summary_stats))
})
