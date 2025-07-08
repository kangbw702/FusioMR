#include <RcppArmadillo.h>
#include "utils.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Implementation of my_rinvwishart
arma::mat my_rinvwishart(double nu, arma::mat S) {
  Function ff("rinvwishart");
  SEXP result = ff(Named("nu") = nu, _["S"] = S);
  arma::mat res = as<arma::mat>(result);
  return res;
}

// Helper function for multivariate normal sampling
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


List gibbs_semo_uhp_only(int niter,
                         NumericVector Gamma_hat_1,
                         NumericVector Gamma_hat_2,
                         NumericVector gamma_hat,
                         NumericVector s2_hat_Gamma_1,
                         NumericVector s2_hat_Gamma_2,
                         NumericVector s2_hat_gamma) {

  // prior setup
  int K = gamma_hat.size();
  arma::vec gamma_hat_vec = as<arma::vec>(gamma_hat);
  arma::vec Gamma_hat_1_vec = as<arma::vec>(Gamma_hat_1);
  arma::vec Gamma_hat_2_vec = as<arma::vec>(Gamma_hat_2);
  arma::vec s2_gamma_vec = as<arma::vec>(s2_hat_gamma);
  arma::vec s2_Gamma_1_vec = as<arma::vec>(s2_hat_Gamma_1);
  arma::vec s2_Gamma_2_vec = as<arma::vec>(s2_hat_Gamma_2);

  double a_prior = 5.0;
  double a_gamma_prior = std::max(a_prior, K / 2.0);
  double a_Gamma_prior = std::max(a_prior, K / 2.0);

  // Gamma prior setup
  double var_gamma = arma::var(gamma_hat_vec);
  double mean_se_gamma = arma::mean(s2_gamma_vec);
  double b_gamma_prior = std::max(1e-3, var_gamma - mean_se_gamma * mean_se_gamma) * (a_gamma_prior - 1.0);

  // Gamma matrix prior setup
  arma::vec var_Gamma = {arma::var(Gamma_hat_1_vec), arma::var(Gamma_hat_2_vec)};
  arma::vec mean_se_Gamma = {arma::mean(s2_Gamma_1_vec), arma::mean(s2_Gamma_2_vec)};
  arma::vec sigma2_theta_prior_mean(2);
  arma::vec b_theta_prior(2);

  for (int i = 0; i < 2; i++) {
    sigma2_theta_prior_mean[i] = std::max(1e-6, var_Gamma[i] - mean_se_Gamma[i] * mean_se_Gamma[i]);
    b_theta_prior[i] = sigma2_theta_prior_mean[i] * (a_Gamma_prior - 1.0);
  }

  double m_Gamma = std::max(a_prior, K / 2.0);
  arma::mat V_Gamma = arma::diagmat(sigma2_theta_prior_mean) * (m_Gamma - 3.0);

  // Initialize starting values
  arma::vec beta_cur = {0.0, 0.0};
  arma::vec gamma_cur = gamma_hat_vec;
  arma::mat theta_cur(K, 2, arma::fill::zeros);  // theta corresponds to Gamma
  arma::mat Sigma_theta_cur = arma::diagmat(sigma2_theta_prior_mean);
  double sigma2_gamma_cur = var_gamma;

  // MCMC trace
  NumericVector beta_1_tk(niter);
  NumericVector beta_2_tk(niter);
  NumericMatrix gamma_tk(niter, K);
  NumericMatrix Gamma_1_tk(niter, K);
  NumericMatrix Gamma_2_tk(niter, K);
  NumericVector sigma2_gamma_tk(niter);
  NumericVector sigma2_theta1_tk(niter);
  NumericVector sigma2_theta2_tk(niter);

  // Store initial values
  beta_1_tk[0] = beta_cur[0];
  beta_2_tk[0] = beta_cur[1];
  sigma2_gamma_tk[0] = sigma2_gamma_cur;
  sigma2_theta1_tk[0] = Sigma_theta_cur(0,0);
  sigma2_theta2_tk[0] = Sigma_theta_cur(1,1);

  for (int k = 0; k < K; k++) {
    gamma_tk(0, k) = gamma_cur[k];
    Gamma_1_tk(0, k) = theta_cur(k, 0);
    Gamma_2_tk(0, k) = theta_cur(k, 1);
  }

  // Gibbs
  for (int iter = 0; iter < (niter-1); iter++) {
    // update gamma_k
    for (int k = 0; k < K; k++) {
      double Ak_gamma = pow(beta_cur[0],2)/s2_Gamma_1_vec[k] +
        pow(beta_cur[1],2)/s2_Gamma_2_vec[k] +
        1.0/s2_gamma_vec[k] + 1.0/sigma2_gamma_cur;
      double Bk_gamma = beta_cur[0]*(Gamma_hat_1_vec[k]-theta_cur(k,0))/s2_Gamma_1_vec[k] +
        beta_cur[1]*(Gamma_hat_2_vec[k]-theta_cur(k,1))/s2_Gamma_2_vec[k] +
        gamma_hat_vec[k]/s2_gamma_vec[k];
      gamma_cur[k] = R::rnorm(Bk_gamma/Ak_gamma, sqrt(1.0/Ak_gamma));
    }

    for (int k = 0; k < K; k++) {
      gamma_tk(iter + 1, k) = gamma_cur[k];
    }

    // update theta_k
    for (int k = 0; k < K; k++) {
      arma::mat S_hat_Gamma_k(2,2, arma::fill::zeros);
      S_hat_Gamma_k(0,0) = s2_Gamma_1_vec[k];
      S_hat_Gamma_k(1,1) = s2_Gamma_2_vec[k];

      arma::mat precision_post_theta_mat = arma::inv(S_hat_Gamma_k) + arma::inv(Sigma_theta_cur);
      arma::vec temp = {Gamma_hat_1_vec[k] - beta_cur[0]*gamma_cur[k],
                        Gamma_hat_2_vec[k] - beta_cur[1]*gamma_cur[k]};
      arma::vec mean_post_theta = arma::inv(precision_post_theta_mat) * arma::inv(S_hat_Gamma_k) * temp;
      arma::mat theta_sample = mvrnormArma(1, mean_post_theta, arma::inv(precision_post_theta_mat));
      theta_cur(k,0) = theta_sample(0,0);
      theta_cur(k,1) = theta_sample(0,1);
    }

    // update sigma2_gamma
    sigma2_gamma_cur = my_rinvgamma(1, a_gamma_prior + K/2.0,
                                    b_gamma_prior + 0.5*arma::sum(arma::pow(gamma_cur,2)))[0];
    sigma2_gamma_cur = std::max(1e-4, std::min(sigma2_gamma_cur, 10.0));
    sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;

    // update Sigma_theta
    arma::mat S_theta(2,2, arma::fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      S_theta = S_theta + arma::trans(theta_cur.row(ii)) * theta_cur.row(ii);
    }

    S_theta = S_theta / K;
    double m_post_theta = m_Gamma + K;
    arma::mat V_post_theta = S_theta * K + V_Gamma;
    arma::mat Sigma_theta_sample = my_rinvwishart(m_post_theta, V_post_theta);

    // add bounds
    double rho_theta = Sigma_theta_sample(1,0)/sqrt(Sigma_theta_sample(0,0))/sqrt(Sigma_theta_sample(1,1));
    Sigma_theta_sample(0,0) = std::max(1e-8, std::min(Sigma_theta_sample(0,0), 10.0));
    Sigma_theta_sample(1,1) = std::max(1e-8, std::min(Sigma_theta_sample(1,1), 10.0));

    Sigma_theta_cur(0,0) = Sigma_theta_sample(0,0);
    Sigma_theta_cur(1,1) = Sigma_theta_sample(1,1);
    Sigma_theta_cur(0,1) = sqrt(Sigma_theta_sample(0,0) * Sigma_theta_sample(1,1)) * rho_theta;
    Sigma_theta_cur(1,0) = Sigma_theta_cur(0,1);

    sigma2_theta1_tk[iter+1] = Sigma_theta_cur(0,0);
    sigma2_theta2_tk[iter+1] = Sigma_theta_cur(1,1);

    // update betas
    arma::mat Omega_hat_Gamma_1(K, K, arma::fill::zeros);
    arma::mat Omega_hat_Gamma_2(K, K, arma::fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      Omega_hat_Gamma_1(ii, ii) = 1.0/s2_Gamma_1_vec[ii];
      Omega_hat_Gamma_2(ii, ii) = 1.0/s2_Gamma_2_vec[ii];
    }
    arma::vec U_beta_1(K), U_beta_2(K), W_beta_1(K), W_beta_2(K);

    for (int ii = 0; ii < K; ii++) {
      U_beta_1[ii] = Gamma_hat_1_vec[ii] - theta_cur(ii,0);
      U_beta_2[ii] = Gamma_hat_2_vec[ii] - theta_cur(ii,1);
      W_beta_1[ii] = gamma_cur[ii];
      W_beta_2[ii] = gamma_cur[ii];
    }
    double A_beta_1 = std::max(arma::as_scalar(arma::trans(W_beta_1) * Omega_hat_Gamma_1 * W_beta_1), 1e-6);
    double A_beta_2 = std::max(arma::as_scalar(arma::trans(W_beta_2) * Omega_hat_Gamma_2 * W_beta_2), 1e-6);
    double mu_beta_1 = arma::as_scalar(arma::trans(W_beta_1) * Omega_hat_Gamma_1 * U_beta_1) / A_beta_1;
    double mu_beta_2 = arma::as_scalar(arma::trans(W_beta_2) * Omega_hat_Gamma_2 * U_beta_2) / A_beta_2;

    beta_cur[0] = R::rnorm(mu_beta_1, sqrt(1.0/A_beta_1));
    beta_cur[1] = R::rnorm(mu_beta_2, sqrt(1.0/A_beta_2));
    beta_1_tk[iter+1] = beta_cur[0];
    beta_2_tk[iter+1] = beta_cur[1];

    // Store Gamma values
    for (int k = 0; k < K; k++) {
      Gamma_1_tk(iter+1, k) = theta_cur(k, 0);
      Gamma_2_tk(iter+1, k) = theta_cur(k, 1);
    }
  }

  List res = List::create(
    Named("beta_1") = beta_1_tk,
    Named("beta_2") = beta_2_tk,
  );
  return res;
}
