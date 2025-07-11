# FusioMR

> **F**lexible, **U**nified and ver**S**atile Mendel**I**an Rand**O**mization framework for GWAS analysis.

A R package for Mendelian Randomization analysis that supports both single and multiple outcome analyses with optional correlated horizontal pleiotropy.

## Overview

With provided IV-exposure effect estimates, IV-outcome effect estimates, standard errors for IV-exposure effect, standard errors for IV-outcome effect, 
FusioMR returns a estimated causal effect.

This R package supports four different analytical models:

- **Model 1**: Single Exposure, Single Outcome, No Correlated Horizontal Pleiotropy
- **Model 2**: Single Exposure, Single Outcome, With Correlated Horizontal Pleiotropy 
- **Model 3**: Single Exposure, Multiple Outcomes, No Correlated Horizontal Pleiotropy
- **Model 4**: Single Exposure, Multiple Outcomes, With Correlated Horizontal Pleiotropy


## Function Reference

###`fusiomr()`

**Purpose**: Perform causal effect estimation with provided GWAS summary statistics.

**Arguments**:

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `b_exp` | numeric vector/matrix | IV-exposure effect estimates | required |
| `se_exp` | numeric vector/matrix | IV-exposure standard errors | required |
| `b_out` | numeric vector/matrix | IV-outcome effect estimates  | required |
| `se_out` | numeric vector/matrix | IV-outcome standard errors  | required |
| `CHP` | logical | Include horizontal pleiotropy correction | `FALSE` |
| `p_value_threshold` | numeric | IV selection p-value threshold | `1e-3` |
| `niter` | integer | Number of Gibbs sampling iterations | `20000` |
| `burnin_prop` | numeric | Proportion of iterations to discard as burn-in | `0.5` |

**Returns**: 
The `fusiomr()` function returns a list containing:

- `est`: Posterior mean(s) of causal effect(s)
- `se`: Posterior standard error(s) 
- `pval`: P-value(s)



## Input Data Requirement
- `b_exp`: Must have exactly 1 column and same number of rows as other inputs.  
- `se_exp`: Must have exactly 1 column, same dimensions as b_exp, and all positive values.  
- `b_out`: Must have 1 column (single outcome) or 2 columns (multiple outcomes), with same number of rows as b_exp.  
- `se_out`: Must have same dimensions as b_out and all positive values.   
- `CHP`: Either FALSE or TRUE.  
- `p_value_threshold`: Must be positive.  
- `niter`: Must be positive integer specifying number of Gibbs sampling iterations.   
- `burnin_prop`:Must be numeric value between 0 and 1. 
 

## Installation

You can install the development version of FusioMR from GitHub with:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install FusioMR from your GitHub repository
devtools::install_github("kangbw702/FusioMR")
```

## Quick Start

### Single Outcome Analysis

```r
library(FusioMR)

# Prepare your data
b_exp <- c(0.12, 0.08, 0.15, 0.09)    # IV-exposure effects
se_exp <- c(0.01, 0.01, 0.02, 0.01)   # IV-exposure standard errors
b_out <- c(0.05, 0.03, 0.07, 0.04)    # IV-outcome effects  
se_out <- c(0.02, 0.02, 0.03, 0.02)   # IV-outcome standard errors

# Run single outcome analysis
result <- fusiomr(
  b_exp = b_exp,
  se_exp = se_exp, 
  b_out = b_out,
  se_out = se_out
)

# View results
print(result)
#> $est
#> [1] 0.4333415
#> 
#> $se  
#> [1] 0.1003104
#> 
#> $pval
#> [1] 1.56024e-05
```

### Multiple Outcomes Analysis

```r
library(FusioMR)

# Prepare your input data
b_exp <- c(0.12, 0.08, 0.15, 0.09)    # IV-exposure effects
se_exp <- c(0.01, 0.01, 0.02, 0.01)   # IV-exposure standard errors
b_out_multi <- cbind(                 # IV-outcome effects 
  outcome1 = c(0.05, 0.03, 0.07, 0.04),
  outcome2 = c(0.08, 0.06, 0.11, 0.07)
)
se_out_multi <- cbind(                # IV-outcome standard errors 
  outcome1 = c(0.02, 0.02, 0.03, 0.02),
  outcome2 = c(0.025, 0.02, 0.03, 0.025)
)

# Run multiple outcomes analysis
result_multi <- fusiomr(
  b_exp = b_exp,
  se_exp = se_exp,
  b_out = b_out_multi,
  se_out = se_out_multi
)

# View results for both outcomes
print(result_multi)
#> $est
#> [1] 0.4466227 0.7525042
#> 
#> $se
#> [1] 0.1798687 0.1836945
#> 
#> $pval
#> [1] 1.298748e-02 4.131123e-05
```

## Advanced Usage

### Custom Parameters

```r
result <- fusiomr(
  b_exp = b_exp, 
  se_exp = se_exp,
  b_out = b_out, 
  se_out = se_out,
  CHP = TRUE,              # Include correlated horizontal pleiotropy
  niter = 50000,           # More iterations
  burnin_prop = 0.3,       # Shorter burn-in
  p_value_threshold = 5e-8 # More stringent threshold
)

result <- fusiomr(
  b_exp = b_exp, 
  se_exp = se_exp, 
  b_out = b_out, 
  se_out = se_out,
  niter = 5000,            # Fewer iterations
  burnin_prop = 0.2,       # Shorter burn-in
  p_value_threshold = 1e-4 # More lenient threshold
)
```

## Examples with Real Data

### Application1
### Application2
### Application3


## Citation



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For questions and support, please contact [kbw@uchicago.edu], [sfeng56@uchicago.edu]
