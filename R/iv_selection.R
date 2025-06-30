#' IV Selection Function
#'
#' @param summary_stats_raw A data frame with 4 columns: b_exp, se_exp, b_out, se_out
#' @param p_threshold P-value threshold for IV selection (default: 1e-3)
#'
#' @return A data frame with selected summary statistics (same 4-column format)
#'
#' @export

summary_stats_selected <- function(summary_stats_raw, p_threshold = 1e-3) {
  # Input validation
  if (!is.data.frame(summary_stats_raw) || ncol(summary_stats_raw) != 4) {
    stop("Wrong input!! \n Provided raw summary statistics must be a data frame with 4 columns: b_exp, se_exp, b_out, se_out")
  }

  required_cols <- c("b_exp", "se_exp", "b_out", "se_out")
  if (!all(required_cols %in% names(summary_stats_raw))) {
    stop("Wrong input!! \n Provided raw summary statistics must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Calculate Z-scores and p-values for exposure
  z_gamma <- abs(summary_stats_raw$b_exp) / summary_stats_raw$se_exp
  p_gamma <- 2 * (1 - pnorm(z_gamma))

  # Select IVs based on p-value threshold
  iv_selected <- p_gamma < p_threshold
  n_selected <- sum(iv_selected)

  if (n_selected == 0) {
    stop("No instrumental variables meet the p-value threshold of ", p_threshold)
  }

  # Create GWAS summary statistics data frame for selected ivs
  selected_stats <- summary_stats_raw[iv_selected, ]
  cat("the number of rows of selected stats", nrow(selected_stats), "\n")
  cat("the number of cols of selected stats", ncol(selected_stats), "\n")

  rownames(selected_stats) <- NULL  # Reset row names
  cat("Selected", n_selected, "IVs out of", nrow(summary_stats_raw), "SNPs\n")

  return (selected_stats)
}
