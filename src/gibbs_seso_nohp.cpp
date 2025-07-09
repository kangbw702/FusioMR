#include <RcppArmadillo.h>
#include "utils.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector gibbs_seso_nohp(int niter,
                              NumericVector Gamma_hat,
                              NumericVector gamma_hat,
                              NumericVector s2_hat_Gamma,
                              NumericVector s2_hat_gamma) {
  // define K
  int K = gamma_hat.size();
  NumericVector beta_tk(niter);

  // MCMC trace
  NumericMatrix gamma_tk(niter, K);
  NumericVector sigma2_gamma_tk(niter);

  // setup priors
  double a_gamma = std::max(2.0, K / 4.0);
  double sum_gamma_var = 0.0;
  for (int k = 0; k < K; k++) {
    sum_gamma_var += std::max(0.0, gamma_hat[k] * gamma_hat[k] - s2_hat_gamma[k]);
  }
  double mean_gamma_var = sum_gamma_var / K;
  double b_gamma = std::max(1e-6, mean_gamma_var * (a_gamma - 1.0));

  // initialization
  double beta_cur = 0.0;
  NumericVector gamma_cur(K);
  double sigma2_gamma_cur = 1.0;
  for (int k = 0; k < K; k++) {
    gamma_cur[k] = rnorm(1, 0, 0.1)[0];
  }

  beta_tk[0] = beta_cur;
  gamma_tk(0, _) = gamma_cur;
  sigma2_gamma_tk[0] = sigma2_gamma_cur;

  // Gibbs
  for (int iter = 0; iter < (niter-1); iter++) {
    // update gamma_k
    for (int k = 0; k < K; k++) {
      double Ak_gamma = beta_cur * beta_cur / s2_hat_Gamma[k] + 1.0 / s2_hat_gamma[k] + 1.0 / sigma2_gamma_cur;
      double Bk_gamma = beta_cur * Gamma_hat[k] / s2_hat_Gamma[k] + gamma_hat[k] / s2_hat_gamma[k];
      gamma_cur[k] = rnorm(1, Bk_gamma / Ak_gamma, sqrt(1.0 / Ak_gamma))[0];
    }
    gamma_tk(iter + 1, _) = gamma_cur;

    // update sigma2_gamma
    double rate_gamma = std::max(1e-6, b_gamma + 0.5 * sum(pow(gamma_cur, 2)));
    sigma2_gamma_cur = my_rinvgamma(1, a_gamma + K / 2.0, rate_gamma)[0];
    // add bounds
    if (sigma2_gamma_cur < 1e-10) sigma2_gamma_cur = 1e-10;
    if (sigma2_gamma_cur > 10) sigma2_gamma_cur = 10;
    sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;

    // update beta
    vec U_beta = Gamma_hat;
    vec W_beta = gamma_cur;
    mat Omega_hat_Gamma(K, K, fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      Omega_hat_Gamma(ii, ii) = 1.0 / s2_hat_Gamma[ii];
    }

    mat A_beta_mat = trans(W_beta) * Omega_hat_Gamma * W_beta;
    double A_beta = std::max(A_beta_mat(0,0), 1e-8);
    mat temp = trans(W_beta) * Omega_hat_Gamma * U_beta;
    double mu_beta = temp(0,0) / A_beta;
    beta_cur = rnorm(1, mu_beta, sqrt(1.0 / A_beta))[0];
    beta_tk[iter + 1] = beta_cur;
  }

  return beta_tk;
}
