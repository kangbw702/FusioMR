#include <RcppArmadillo.h>
#include "utils.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat my_rinvwishart(double nu, const arma::mat S);
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
Rcpp::NumericVector my_rdirichlet(int n, Rcpp::NumericVector alpha);

// [[Rcpp::export]]
List gibbs_memo_joint(int niter,
                      NumericVector Gamma_hat_1,
                      NumericVector Gamma_hat_2,
                      NumericVector gamma_hat_1,
                      NumericVector gamma_hat_2,
                      NumericVector s_hat_Gamma_1,
                      NumericVector s_hat_Gamma_2,
                      NumericVector s_hat_gamma_1,
                      NumericVector s_hat_gamma_2
) {
  int K = gamma_hat_1.size();
  NumericVector s2_hat_Gamma_1 = pow(s_hat_Gamma_1,2);
  NumericVector s2_hat_Gamma_2 = pow(s_hat_Gamma_2,2);
  NumericVector s2_hat_gamma_1 = pow(s_hat_gamma_1,2);
  NumericVector s2_hat_gamma_2 = pow(s_hat_gamma_2,2);

  if (Gamma_hat_1.size() != K || Gamma_hat_2.size() != K ||
      s2_hat_Gamma_1.size() != K || s2_hat_Gamma_2.size() != K ||
      s2_hat_gamma_1.size() != K || s2_hat_gamma_2.size() != K) {
    Rcpp::stop("Input vectors must all be the same length");
  }

  // Initial values setup
  NumericVector eta_1_init(K, 0.0);
  NumericVector eta_2_init(K, 0.0);
  double alpha_1_init = 1.0;
  double alpha_2_init = 1.0;
  double beta_1_init = 0.0;
  double beta_2_init = 0.0;
  double a_prior = 5.0;
  double m_gamma = std::max(a_prior, static_cast<double>(K)/2.0);
  double m_theta = std::max(a_prior, static_cast<double>(K)/2.0);

  double sigma2_gamma_prior_mean_1 = std::max(1e-3, (var(gamma_hat_1) - pow(mean(s_hat_gamma_1), 2)));
  double sigma2_gamma_prior_mean_2 = std::max(1e-3, (var(gamma_hat_2) - pow(mean(s_hat_gamma_2), 2)));
  arma::mat Sigma_gamma_init = arma::eye<arma::mat>(2, 2);
  Sigma_gamma_init(0,0) = sigma2_gamma_prior_mean_1;
  Sigma_gamma_init(1,1) = sigma2_gamma_prior_mean_2;

  double sigma2_theta_prior_mean_1 = std::max(1e-6, (var(Gamma_hat_1) - pow(mean(s_hat_Gamma_1), 2)));
  double sigma2_theta_prior_mean_2 = std::max(1e-6, (var(Gamma_hat_2) - pow(mean(s_hat_Gamma_2), 2)));
  arma::mat Sigma_theta_init = arma::eye<arma::mat>(2, 2);
  Sigma_theta_init(0,0) = sigma2_theta_prior_mean_1;
  Sigma_theta_init(1,1) = sigma2_theta_prior_mean_2;

  arma::mat V_gamma = Sigma_gamma_init * (m_gamma - 3);
  arma::mat V_theta = Sigma_theta_init * (m_theta - 3);

  double rho_eta = 0;
  double q_chp1 = 0;
  double q_chp2 = 0;
  double p00 = std::max(0.0, rho_eta * sqrt(q_chp1 * q_chp2 * (1 - q_chp1) * (1 - q_chp2)) + (1 - q_chp1) * (1 - q_chp2));
  double p10 = std::max(0.0, 1 - q_chp2 - p00);
  double p01 = std::max(0.0, 1 - q_chp1 - p00);
  double p11 = std::max(0.0, q_chp1 + q_chp2 + p00 - 1);

  // MCMC trace
  NumericMatrix theta_1_tk(niter, K);
  NumericMatrix theta_2_tk(niter, K);
  NumericMatrix gamma_1_tk(niter, K);
  NumericMatrix gamma_2_tk(niter, K);
  NumericMatrix eta_1_tk(niter, K);
  NumericMatrix eta_2_tk(niter, K);
  NumericMatrix pst_tk(niter, 4);
  NumericVector alpha_1_tk(niter);
  NumericVector alpha_2_tk(niter);
  NumericVector beta_1_tk(niter);
  NumericVector beta_2_tk(niter);
  NumericMatrix q_post_1_tk(niter, K);
  NumericMatrix q_post_2_tk(niter, K);
  NumericVector cc = NumericVector::create(1.0, 1.0, 1.0, 1.0);
  NumericVector pst_init = NumericVector::create(p00, p01, p10, p11);

  // Initialize the starting values
  for (int j = 0; j < K; ++j) {
    theta_1_tk(0, j) = 0.0;
    theta_2_tk(0, j) = 0.0;
    gamma_1_tk(0, j) = 0.0;
    gamma_2_tk(0, j) = 0.0;
    eta_1_tk(0, j) = eta_1_init[j];
    eta_2_tk(0, j) = eta_2_init[j];
  }
  for (int j = 0; j < 4; ++j) pst_tk(0, j) = pst_init[j];
  alpha_1_tk[0] = alpha_1_init;
  alpha_2_tk[0] = alpha_2_init;
  beta_1_tk[0] = beta_1_init;
  beta_2_tk[0] = beta_2_init;

  std::vector<double> sigma2_gamma1_tk(niter);
  std::vector<double> sigma2_gamma2_tk(niter);
  std::vector<double> sigma2_theta1_tk(niter);
  std::vector<double> sigma2_theta2_tk(niter);
  sigma2_gamma1_tk[0] = Sigma_gamma_init(0,0);
  sigma2_gamma2_tk[0] = Sigma_gamma_init(1,1);
  sigma2_theta1_tk[0] = Sigma_theta_init(0,0);
  sigma2_theta2_tk[0] = Sigma_theta_init(1,1);

  arma::mat theta_cur(K, 2, fill::zeros);
  arma::mat gamma_cur(K, 2, fill::zeros);
  arma::mat eta_cur(K, 2, fill::zeros);
  for (int k = 0; k < K; ++k) {
    eta_cur(k,0) = eta_1_init[k];
    eta_cur(k,1) = eta_2_init[k];
  }
  arma::vec alpha_cur = {alpha_1_init, alpha_2_init};
  arma::vec beta_cur = {beta_1_init, beta_2_init};
  arma::mat Sigma_gamma_cur = Sigma_gamma_init;
  arma::mat Sigma_theta_cur = Sigma_theta_init;
  NumericVector q_post_1_cur(K);
  NumericVector q_post_2_cur(K);
  NumericVector pst_cur = clone(pst_init);

  // Gibbs
  for (int iter = 0; iter < (niter-1); iter++) {
    for (int k = 0; k < K; k++) {
      // Update gamma_k
      arma::mat A_k(2, 2, fill::zeros);
      A_k(0,0) = beta_cur[0] + alpha_cur[0] * eta_cur(k,0);
      A_k(1,1) = beta_cur[1] + alpha_cur[1] * eta_cur(k,1);
      arma::mat S_hat_Gamma_k_inv(2, 2, fill::zeros);
      S_hat_Gamma_k_inv(0,0) = 1.0 / s2_hat_Gamma_1[k];
      S_hat_Gamma_k_inv(1,1) = 1.0 / s2_hat_Gamma_2[k];
      arma::mat S_hat_gamma_k_inv(2, 2, fill::zeros);
      S_hat_gamma_k_inv(0,0) = 1.0 / s2_hat_gamma_1[k];
      S_hat_gamma_k_inv(1,1) = 1.0 / s2_hat_gamma_2[k];
      arma::mat Sigma_gamma_cur_inv = inv_sympd(Sigma_gamma_cur);
      arma::mat precision_mat_post_gamma = A_k.t() * S_hat_Gamma_k_inv * A_k + S_hat_gamma_k_inv + Sigma_gamma_cur_inv;
      arma::vec temp1 = {Gamma_hat_1[k] - theta_cur(k,0), Gamma_hat_2[k] - theta_cur(k,1)};
      arma::vec temp2 = {gamma_hat_1[k], gamma_hat_2[k]};
      arma::vec mu_post_nom_gamma = A_k.t() * S_hat_Gamma_k_inv * temp1 + S_hat_gamma_k_inv * temp2;
      arma::mat cov_post_gamma = inv_sympd(precision_mat_post_gamma);
      arma::vec mu_post_gamma = cov_post_gamma * mu_post_nom_gamma;
      arma::mat gamma_sample = mvrnormArma(1, mu_post_gamma, cov_post_gamma);
      gamma_cur(k,0) = gamma_sample(0,0);
      gamma_cur(k,1) = gamma_sample(0,1);

      // Update theta_k
      arma::mat Sigma_theta_cur_inv = inv_sympd(Sigma_theta_cur);
      arma::mat precision_mat_post_theta = S_hat_Gamma_k_inv + Sigma_theta_cur_inv;
      arma::vec temp3 = {Gamma_hat_1[k] - beta_cur[0]*gamma_cur(k,0) - alpha_cur[0]*gamma_cur(k,0)*eta_cur(k,0),
                         Gamma_hat_2[k] - beta_cur[1]*gamma_cur(k,1) - alpha_cur[1]*gamma_cur(k,1)*eta_cur(k,1)};
      arma::vec mu_post_nom_theta = S_hat_Gamma_k_inv * temp3;
      arma::mat cov_post_theta = inv_sympd(precision_mat_post_theta);
      arma::vec mu_post_theta = cov_post_theta * mu_post_nom_theta;
      arma::mat theta_sample = mvrnormArma(1, mu_post_theta, cov_post_theta);
      theta_cur(k,0) = theta_sample(0,0);
      theta_cur(k,1) = theta_sample(0,1);
    }

    // Update Sigma_gamma
    arma::mat S_gamma(2, 2, fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      arma::vec gamma_ii = gamma_cur.row(ii).t();
      S_gamma += gamma_ii * gamma_ii.t();
    }
    S_gamma /= K;
    double m_post_gamma = m_gamma + K;
    arma::mat V_post_gamma = S_gamma * K + V_gamma;
    V_post_gamma = 0.5 * (V_post_gamma + V_post_gamma.t());
    arma::mat Sigma_gamma_sample = my_rinvwishart(m_post_gamma, V_post_gamma);
    double sigma11_gamma = std::max(1e-4, std::min(10.0, Sigma_gamma_sample(0,0)));
    double sigma22_gamma = std::max(1e-4, std::min(10.0, Sigma_gamma_sample(1,1)));
    double rho_gamma = Sigma_gamma_sample(0,1) / (sqrt(sigma11_gamma) * sqrt(sigma22_gamma));
    rho_gamma = std::max(-0.99, std::min(0.99, rho_gamma));
    Sigma_gamma_cur(0,0) = sigma11_gamma;
    Sigma_gamma_cur(1,1) = sigma22_gamma;
    Sigma_gamma_cur(0,1) = rho_gamma * sqrt(sigma11_gamma * sigma22_gamma);
    Sigma_gamma_cur(1,0) = Sigma_gamma_cur(0,1);
    sigma2_gamma1_tk[iter+1] = Sigma_gamma_cur(0,0);
    sigma2_gamma2_tk[iter+1] = Sigma_gamma_cur(1,1);

    // Update Sigma_theta
    arma::mat S_theta(2, 2, fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      arma::vec theta_ii = theta_cur.row(ii).t();
      S_theta += theta_ii * theta_ii.t();
    }
    S_theta /= K;
    double m_post_theta = m_theta + K;
    arma::mat V_post_theta = S_theta * K + V_theta;
    V_post_theta = 0.5 * (V_post_theta + V_post_theta.t());
    arma::mat Sigma_theta_sample = my_rinvwishart(m_post_theta, V_post_theta);
    double sigma11_theta = std::max(1e-8, std::min(10.0, Sigma_theta_sample(0,0)));
    double sigma22_theta = std::max(1e-8, std::min(10.0, Sigma_theta_sample(1,1)));
    double rho_theta = Sigma_theta_sample(0,1) / (sqrt(sigma11_theta) * sqrt(sigma22_theta));
    rho_theta = std::max(-0.99, std::min(0.99, rho_theta));
    Sigma_theta_cur(0,0) = sigma11_theta;
    Sigma_theta_cur(1,1) = sigma22_theta;
    Sigma_theta_cur(0,1) = rho_theta * sqrt(sigma11_theta * sigma22_theta);
    Sigma_theta_cur(1,0) = Sigma_theta_cur(0,1);

    sigma2_theta1_tk[iter+1] = Sigma_theta_cur(0,0);
    sigma2_theta2_tk[iter+1] = Sigma_theta_cur(1,1);

    // Update eta_k
    for (int k = 0; k < K; k++) {
      // eta_1
      double fk1_part1 = -0.5 * pow(Gamma_hat_1[k] - theta_cur(k,0) - (beta_cur[0] + alpha_cur[0]) * gamma_cur(k,0), 2) / s2_hat_Gamma_1[k];
      double fk1_part2 = pow(pst_cur[2], (eta_cur(k,1) == 0 ? 1.0 : 0.0)) * pow(pst_cur[3], (eta_cur(k,1) == 1 ? 1.0 : 0.0));
      double fk1 = exp(fk1_part1) * fk1_part2;
      double fk0_part1 = -0.5 * pow(Gamma_hat_1[k] - theta_cur(k,0) - beta_cur[0] * gamma_cur(k,0), 2) / s2_hat_Gamma_1[k];
      double fk0_part2 = pow(pst_cur[0], (eta_cur(k,1) == 0 ? 1.0 : 0.0)) * pow(pst_cur[1], (eta_cur(k,1) == 1 ? 1.0 : 0.0));
      double fk0 = exp(fk0_part1) * fk0_part2;
      double q1_post = fk1 / (fk1 + fk0);
      if (!R_finite(q1_post) || q1_post < 1e-10) q1_post = 0.0;
      if (q1_post > 1.0 - 1e-10) q1_post = 1.0;
      eta_cur(k,0) = R::rbinom(1, q1_post);
      q_post_1_cur[k] = q1_post;
      eta_1_tk(iter+1, k) = eta_cur(k,0);

      // eta_2
      double gk1_part1 = -0.5 * pow(Gamma_hat_2[k] - theta_cur(k,1) - (beta_cur[1] + alpha_cur[1]) * gamma_cur(k,1), 2) / s2_hat_Gamma_2[k];
      double gk1_part2 = pow(pst_cur[1], (eta_cur(k,0) == 0 ? 1.0 : 0.0)) * pow(pst_cur[3], (eta_cur(k,0) == 1 ? 1.0 : 0.0));
      double gk1 = exp(gk1_part1) * gk1_part2;
      double gk0_part1 = -0.5 * pow(Gamma_hat_2[k] - theta_cur(k,1) - beta_cur[1] * gamma_cur(k,1), 2) / s2_hat_Gamma_2[k];
      double gk0_part2 = pow(pst_cur[0], (eta_cur(k,0) == 0 ? 1.0 : 0.0)) * pow(pst_cur[2], (eta_cur(k,0) == 1 ? 1.0 : 0.0));
      double gk0 = exp(gk0_part1) * gk0_part2;
      double q2_post = gk1 / (gk1 + gk0);
      if (!R_finite(q2_post) || q2_post < 1e-10) q2_post = 0.0;
      if (q2_post > 1.0 - 1e-10) q2_post = 1.0;
      eta_cur(k,1) = R::rbinom(1, q2_post);
      q_post_2_cur[k] = q2_post;
      eta_2_tk(iter+1, k) = eta_cur(k,1);
    }

    q_post_1_tk(iter+1, _) = q_post_1_cur;
    q_post_2_tk(iter+1, _) = q_post_2_cur;

    // Update pst
    int n00 = 0, n01 = 0, n10 = 0, n11 = 0;
    for (int k = 0; k < K; k++) {
      if (eta_cur(k,0) == 0 && eta_cur(k,1) == 0) n00++;
      else if (eta_cur(k,0) == 0 && eta_cur(k,1) == 1) n01++;
      else if (eta_cur(k,0) == 1 && eta_cur(k,1) == 0) n10++;
      else if (eta_cur(k,0) == 1 && eta_cur(k,1) == 1) n11++;
    }
    NumericVector dir_par = {n00 + cc[0], n01 + cc[1], n10 + cc[2], n11 + cc[3]};
    pst_cur = my_rdirichlet(1, dir_par);
    for (int j = 0; j < 4; j++) {
      pst_tk(iter+1, j) = pst_cur[j];
    }

    // Update alpha and beta
    arma::vec gamma_1_vec(K), gamma_2_vec(K);
    arma::vec eta_1_vec(K), eta_2_vec(K);
    for (int k = 0; k < K; k++) {
      gamma_1_vec[k] = gamma_cur(k,0);
      gamma_2_vec[k] = gamma_cur(k,1);
      eta_1_vec[k] = eta_cur(k,0);
      eta_2_vec[k] = eta_cur(k,1);
    }

    // Update alpha
    arma::vec U_alpha_1(K), U_alpha_2(K), W_alpha_1(K), W_alpha_2(K);
    for (int k = 0; k < K; k++) {
      U_alpha_1[k] = Gamma_hat_1[k] - theta_cur(k,0) - beta_cur[0] * gamma_cur(k,0);
      U_alpha_2[k] = Gamma_hat_2[k] - theta_cur(k,1) - beta_cur[1] * gamma_cur(k,1);
      W_alpha_1[k] = eta_cur(k,0) * gamma_cur(k,0);
      W_alpha_2[k] = eta_cur(k,1) * gamma_cur(k,1);
    }
    // alpha_1
    double sum_W_alpha_1 = arma::sum(W_alpha_1);
    if (sum_W_alpha_1 < 1e-10) {
      alpha_cur[0] = 0.0;
    } else {
      double A_alpha_1 = 0.0;
      double B_alpha_1 = 0.0;
      for (int k = 0; k < K; k++) {
        A_alpha_1 += W_alpha_1[k] * W_alpha_1[k] / s2_hat_Gamma_1[k];
        B_alpha_1 += W_alpha_1[k] * U_alpha_1[k] / s2_hat_Gamma_1[k];
      }
      double mu_alpha_1 = B_alpha_1 / A_alpha_1;
      double var_alpha_1 = 1.0 / A_alpha_1;
      alpha_cur[0] = R::rnorm(mu_alpha_1, sqrt(var_alpha_1));
    }
    // alpha_2
    double sum_W_alpha_2 = arma::sum(W_alpha_2);
    if (sum_W_alpha_2 < 1e-10) {
      alpha_cur[1] = 0.0;
    } else {
      double A_alpha_2 = 0.0;
      double B_alpha_2 = 0.0;
      for (int k = 0; k < K; k++) {
        A_alpha_2 += W_alpha_2[k] * W_alpha_2[k] / s2_hat_Gamma_2[k];
        B_alpha_2 += W_alpha_2[k] * U_alpha_2[k] / s2_hat_Gamma_2[k];
      }
      double mu_alpha_2 = B_alpha_2 / A_alpha_2;
      double var_alpha_2 = 1.0 / A_alpha_2;
      alpha_cur[1] = R::rnorm(mu_alpha_2, sqrt(var_alpha_2));
    }

    alpha_1_tk[iter+1] = alpha_cur[0];
    alpha_2_tk[iter+1] = alpha_cur[1];

    // Update beta
    arma::vec U_beta_1(K), U_beta_2(K), W_beta_1(K), W_beta_2(K);
    for (int k = 0; k < K; k++) {
      U_beta_1[k] = Gamma_hat_1[k] - theta_cur(k,0) - alpha_cur[0] * eta_cur(k,0) * gamma_cur(k,0);
      U_beta_2[k] = Gamma_hat_2[k] - theta_cur(k,1) - alpha_cur[1] * eta_cur(k,1) * gamma_cur(k,1);
      W_beta_1[k] = gamma_cur(k,0);
      W_beta_2[k] = gamma_cur(k,1);
    }
    double A_beta_1 = 0.0;
    double B_beta_1 = 0.0;
    for (int k = 0; k < K; k++) {
      A_beta_1 += W_beta_1[k] * W_beta_1[k] / s2_hat_Gamma_1[k];
      B_beta_1 += W_beta_1[k] * U_beta_1[k] / s2_hat_Gamma_1[k];
    }
    double mu_beta_1 = B_beta_1 / A_beta_1;
    double var_beta_1 = 1.0 / A_beta_1;
    beta_cur[0] = R::rnorm(mu_beta_1, sqrt(var_beta_1));
    double A_beta_2 = 0.0;
    double B_beta_2 = 0.0;
    for (int k = 0; k < K; k++) {
      A_beta_2 += W_beta_2[k] * W_beta_2[k] / s2_hat_Gamma_2[k];
      B_beta_2 += W_beta_2[k] * U_beta_2[k] / s2_hat_Gamma_2[k];
    }
    double mu_beta_2 = B_beta_2 / A_beta_2;
    double var_beta_2 = 1.0 / A_beta_2;
    beta_cur[1] = R::rnorm(mu_beta_2, sqrt(var_beta_2));

    beta_1_tk[iter+1] = beta_cur[0];
    beta_2_tk[iter+1] = beta_cur[1];

    // Store current values
    for (int k = 0; k < K; k++) {
      theta_1_tk(iter+1, k) = theta_cur(k,0);
      theta_2_tk(iter+1, k) = theta_cur(k,1);
      gamma_1_tk(iter+1, k) = gamma_cur(k,0);
      gamma_2_tk(iter+1, k) = gamma_cur(k,1);
    }
  }

  List res = List::create(
    Named("beta_1_tk") = beta_1_tk,
    Named("beta_2_tk") = beta_2_tk,
    Named("alpha_1_tk") = alpha_1_tk,
    Named("alpha_2_tk") = alpha_2_tk,
    Named("eta_1_tk") = eta_1_tk,
    Named("eta_2_tk") = eta_2_tk,
    Named("pst_tk") = pst_tk
  );
  return res;
}
