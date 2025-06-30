---
title: "README"
output: html_document
---

# FusioMR

FusioMR implements fusion framework for Mendelian Randomization analysis.

## Installation

You can install the development version of FusioMR from [GitHub](https://github.com/kangbw702/FusioMR/) with:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install FusioMR from your GitHub repository
devtools::install_github("kangbw702/FusioMR")
```

## Quick Start

```r
library(FusioMR)

# Create example summary statistics
summary_stats <- data.frame(
  b_exp = rnorm(100, 0, 0.1),
  se_exp = runif(100, 0.01, 0.05),
  b_out = rnorm(100, 0, 0.1),  
  se_out = runif(100, 0.01, 0.05)
)

# Run FusioMR analysis
result <- fusiomr(
  summary_stats = summary_stats,
  type = "seso_uhp_only",
  p_value_threshold = 1e-3,
  niter = 10000,
  burnin_prop = 0.5
)

# View results
print(result)
summary(result)
```

## Methods

FusioMR currently supports the following model types:

- **seso_uhp_only**: Single exposure, single outcome with uncorrelated horizontal pleiotropy only
- **seso_nohp**: Single exposure, single outcome with no horizontal pleiotropy
- **multi_uhp_only**: Multiple exposures with uncorrelated horizontal pleiotropy (coming soon)
- **multi_nohp**: Multiple exposures with no horizontal pleiotropy (coming soon)

## Input Data Format

The input `summary_stats` should be a data frame with four columns:

- `b_exp`: Effect sizes for exposure (beta coefficients)
- `se_exp`: Standard errors for exposure effects
- `b_out`: Effect sizes for outcome (beta coefficients)  
- `se_out`: Standard errors for outcome effects

Each row represents a genetic variant (SNP).

## Output

The `fusiomr()` function returns a list containing:

- `beta_estimate`: Posterior mean of the causal effect
- `beta_se`: Posterior standard error
- `beta_ci`: 95% credible interval
- `beta_samples`: MCMC samples (post burn-in)
- `type`: Model type used
- `n_ivs`: Number of instrumental variables selected
- Additional metadata about the analysis

## Citation

If you use FusioMR in your research, please cite:

```
[...citation information]
```

## License

This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions and support, please contact [...] or open an issue on GitHub.

