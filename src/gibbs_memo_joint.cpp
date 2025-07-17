#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Implementation of my_rinvwishart
// [[Rcpp::export]]
arma::mat my_rinvwishart(double nu, arma::mat S) {
  return arma::iwishrnd(S, nu);
}

// Helper function for multivariate normal sampling with robust Cholesky
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  sigma = 0.5 * (sigma + sigma.t());
  arma::mat L;
  bool success = arma::chol(L, sigma, "lower");
  if (!success) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, sigma);
    eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
    L = eigvec * arma::diagmat(arma::sqrt(eigval));
  }
  return arma::repmat(mu, 1, n).t() + Y * L.t();
}

// [[Rcpp::export]]
NumericVector my_rdirichlet(int n, NumericVector alpha) {
  Function ff("rdirichlet");
  NumericVector res = ff(Named("n") = n, _["alpha"] = alpha);
  return res;
}

// Gibbs sampling

// [[Rcpp::export]]
List gibbs_memo_joint(int niter,
                      vec Gamma_hat_1,
                      vec Gamma_hat_2,
                      vec gamma_hat_1,
                      vec gamma_hat_2,
                      vec s2_hat_Gamma_1,
                      vec s2_hat_Gamma_2,
                      vec s2_hat_gamma_1,
                      vec s2_hat_gamma_2,
                      double rho_eta,
                      double q_chp1,
                      double q_chp2
) {
  int K = gamma_hat_1.size();
  if (Gamma_hat_1.size() != K || Gamma_hat_2.size() != K ||
      s2_hat_Gamma_1.size() != K || s2_hat_Gamma_2.size() != K || s2_hat_gamma_1.size() != K || s2_hat_gamma_2.size() != K) {
    Rcpp::stop("Input vectors to gibbs_semo_nohp() must all be the same length");
  }

  // initial values setup
  NumericVector eta_1_init(K, 0.0);
  NumericVector eta_2_init(K, 0.0);
  double alpha_1_init = 1.0;
  double alpha_2_init = 1.0;
  double beta_1_init = 0.0;
  double beta_2_init = 0.0;
  double a_prior = 5.0;
  double m_gamma = std::max(a_prior, K/2.0);
  double m_theta = std::max(a_prior, K/2.0);
  double sigma2_gamma_prior_mean_1 = std::max(1e-3, (var(gamma_hat_1) - pow(mean(s2_hat_gamma_1), 2)));
  double sigma2_gamma_prior_mean_2 = std::max(1e-3, (var(gamma_hat_2) - pow(mean(s2_hat_gamma_2), 2)));
  arma::mat Sigma_gamma_init = arma::eye<arma::mat>(2, 2);
  Sigma_gamma_init(0,0) = sigma2_gamma_prior_mean_1;
  Sigma_gamma_init(1,1) = sigma2_gamma_prior_mean_2;
  double sigma2_theta_prior_mean_1 = std::max(1e-6, (var(Gamma_hat_1) - pow(mean(s2_hat_Gamma_1), 2)));
  double sigma2_theta_prior_mean_2 = std::max(1e-6, (var(Gamma_hat_2) - pow(mean(s2_hat_Gamma_2), 2)));
  arma::mat Sigma_theta_init = arma::eye<arma::mat>(2, 2);
  Sigma_theta_init(0,0) = sigma2_theta_prior_mean_1;
  Sigma_theta_init(1,1) = sigma2_theta_prior_mean_2;
  arma::mat V_gamma = Sigma_gamma_init * (m_gamma - 3);
  arma::mat V_theta = Sigma_theta_init * (m_theta - 3);
  double p00 = std::max(0.0, rho_eta * sqrt(q_chp1 * q_chp2 * (1 - q_chp1) * (1 - q_chp2)) + (1 - q_chp1) * (1 - q_chp2));
  double p10 = std::max(0.0,1 - q_chp2 - p00);
  double p01 = std::max(0.0,1 - q_chp1 - p00);
  double p11 = std::max(0.0,q_chp1 + q_chp2 + p00 - 1);

  // MCMC trace
  NumericMatrix theta_1_tk(niter, K);
  NumericMatrix theta_2_tk(niter, K);
  NumericMatrix gamma_1_tk(niter, K);
  NumericMatrix gamma_2_tk(niter, K);
  NumericMatrix eta_1_tk(niter, K);
  NumericMatrix eta_2_tk(niter, K);
  NumericMatrix pst_tk(niter, 4);
  NumericVector alpha_1_tk(niter, NA_REAL);
  NumericVector alpha_2_tk(niter, NA_REAL);
  NumericVector beta_1_tk(niter, NA_REAL);
  NumericVector beta_2_tk(niter, NA_REAL);
  NumericMatrix q_post_1_tk(niter, K);
  NumericMatrix q_post_2_tk(niter, K);
  NumericVector cc = NumericVector::create(1.0, 1.0, 1.0, 1.0);
  NumericVector pst_init = NumericVector::create(p00, p01, p10, p11);

  // set initial values
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

  // set current values as initial values
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
    // Update gamma_k and theta_k
    for (int k = 0; k < K; k++) {
      arma::mat A_k(2, 2, fill::zeros);
      A_k(0,0) = beta_cur[0] + alpha_cur[0] * eta_cur(k,0);
      A_k(1,1) = beta_cur[1] + alpha_cur[1] * eta_cur(k,1);
      arma::mat S_hat_Gamma_k(2, 2, fill::zeros);
      S_hat_Gamma_k(0,0) = s2_hat_Gamma_1[k];
      S_hat_Gamma_k(1,1) = s2_hat_Gamma_2[k];
      arma::mat S_hat_gamma_k(2, 2, fill::zeros);
      S_hat_gamma_k(0,0) = s2_hat_gamma_1[k];
      S_hat_gamma_k(1,1) = s2_hat_gamma_2[k];

      // Update gamma_k
      arma::mat precision_mat_post_gamma = A_k.t() * (1/S_hat_Gamma_k) * A_k + (1/S_hat_gamma_k) + (1/Sigma_gamma_cur);
      arma::vec temp1 = {Gamma_hat_1[k] - theta_cur(k,0), Gamma_hat_2[k] - theta_cur(k,1)};
      arma::vec temp2 = {gamma_hat_1[k], gamma_hat_2[k]};
      arma::vec mu_post_nom_gamma = A_k.t() * (1/S_hat_Gamma_k) * temp1 + (1/S_hat_gamma_k) * temp2;
      arma::vec mu_post_gamma = (1/precision_mat_post_gamma) * mu_post_nom_gamma;
      arma::mat gamma_sample = mvrnormArma(1, mu_post_gamma, (1/precision_mat_post_gamma));
      gamma_cur(k,0) = gamma_sample(0,0);
      gamma_cur(k,1) = gamma_sample(0,1);

      // Update theta_k
      arma::mat precision_mat_post_theta = (1/S_hat_Gamma_k) + (1/Sigma_theta_cur);
      arma::vec temp3 = {Gamma_hat_1[k] - beta_cur[0]*gamma_cur(k,0) - alpha_cur[0]*gamma_cur(k,0)*eta_cur(k,0),
                         Gamma_hat_2[k] - beta_cur[1]*gamma_cur(k,1) - alpha_cur[1]*gamma_cur(k,1)*eta_cur(k,1)};
      arma::vec mu_post_nom_theta = (1/S_hat_Gamma_k) * temp3;
      arma::vec mu_post_theta = (1/precision_mat_post_theta) * mu_post_nom_theta;
      arma::mat theta_sample = mvrnormArma(1, mu_post_theta, (1/precision_mat_post_theta));
      theta_cur(k,0) = theta_sample(0,0);
      theta_cur(k,1) = theta_sample(0,1);
    }

    // Update Sigma_gamma
    arma::mat S_gamma(2, 2, fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      arma::rowvec rowk = gamma_cur.row(ii);
      S_gamma += rowk.t() * rowk;
    }
    S_gamma /= K;
    double m_post_gamma = m_gamma + K;
    arma::mat V_post_gamma = S_gamma * K + V_gamma;
    arma::mat Sigma_gamma_sample = my_rinvwishart(m_post_gamma, V_post_gamma);
    // Add bounds
    double rho_gamma = Sigma_gamma_sample(1,0) / (sqrt(Sigma_gamma_sample(0,0)) * sqrt(Sigma_gamma_sample(1,1)));
    double sigma11_gamma = std::max(1e-4, std::min(10.0, Sigma_gamma_sample(0,0)));
    double sigma22_gamma = std::max(1e-4, std::min(10.0, Sigma_gamma_sample(1,1)));
    Sigma_gamma_cur(0,0) = sigma11_gamma;
    Sigma_gamma_cur(1,1) = sigma22_gamma;
    Sigma_gamma_cur(0,1) = sigma11_gamma * sigma22_gamma * rho_gamma;
    Sigma_gamma_cur(1,0) = Sigma_gamma_cur(0,1);
    sigma2_gamma1_tk[iter+1] = Sigma_gamma_cur(0,0);
    sigma2_gamma2_tk[iter+1] = Sigma_gamma_cur(1,1);

    // Update Sigma_theta
    arma::mat S_theta(2, 2, fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      arma::rowvec rowk = theta_cur.row(ii);
      S_theta += rowk.t() * rowk;
    }
    S_theta /= K;
    double m_post_theta = m_theta + K;
    arma::mat V_post_theta = S_theta * K + V_theta;
    arma::mat Sigma_theta_sample = my_rinvwishart(m_post_theta, V_post_theta);
    // Add bounds
    double rho_theta = Sigma_theta_sample(1,0) / (sqrt(Sigma_theta_sample(0,0)) * sqrt(Sigma_theta_sample(1,1)));
    double sigma11_theta = std::max(1e-8, std::min(10.0, Sigma_theta_sample(0,0)));
    double sigma22_theta = std::max(1e-8, std::min(10.0, Sigma_theta_sample(1,1)));
    Sigma_theta_cur(0,0) = sigma11_theta;
    Sigma_theta_cur(1,1) = sigma22_theta;
    Sigma_theta_cur(0,1) = sigma11_theta * sigma22_theta * rho_theta;
    Sigma_theta_cur(1,0) = Sigma_theta_cur(0,1);
    sigma2_theta1_tk[iter+1] = Sigma_theta_cur(0,0);
    sigma2_theta2_tk[iter+1] = Sigma_theta_cur(1,1);

    // Update eta_k (latent indicators)
    for (int k = 0; k < K; k++) {
      double fk1_part1 = -0.5 * pow(Gamma_hat_1[k] - theta_cur(k,0) - (beta_cur[0] + alpha_cur[0]) * gamma_cur(k,0), 2) / s2_hat_Gamma_1[k];
      double fk1_part2 = pow(pst_cur[2], (eta_cur(k,1) == 0 ? 1.0 : 0.0)) * pow(pst_cur[3], (eta_cur(k,1) == 1 ? 1.0 : 0.0));
      double fk1 = exp(fk1_part1) * fk1_part2;
      double fk0_part1 = -0.5 * pow(Gamma_hat_1[k] - theta_cur(k,0) - beta_cur[0] * gamma_cur(k,0), 2) / s2_hat_Gamma_1[k];
      double fk0_part2 = pow(pst_cur[0], (eta_cur(k,1) == 0 ? 1.0 : 0.0)) * pow(pst_cur[1], (eta_cur(k,1) == 1 ? 1.0 : 0.0));
      double fk0 = exp(fk0_part1) * fk0_part2;
      double q1_post = fk1 / (fk1 + fk0);
      // Add bounds
      if (q1_post <= 0.1) q1_post = 0.0;
      if (q1_post >= 0.9) q1_post = 1.0;
      eta_cur(k,0) = R::rbinom(1, q1_post);
      q_post_1_cur[k] = q1_post;
      eta_1_tk(iter+1, k) = eta_cur(k,0);

      double gk1_part1 = -0.5 * pow(Gamma_hat_2[k] - theta_cur(k,1) - (beta_cur[1] + alpha_cur[1]) * gamma_cur(k,1), 2) / s2_hat_Gamma_2[k];
      double gk1_part2 = pow(pst_cur[1], (eta_cur(k,0) == 0 ? 1.0 : 0.0)) * pow(pst_cur[3], (eta_cur(k,0) == 1 ? 1.0 : 0.0));
      double gk1 = exp(gk1_part1) * gk1_part2;
      double gk0_part1 = -0.5 * pow(Gamma_hat_2[k] - theta_cur(k,1) - beta_cur[1] * gamma_cur(k,1), 2) / s2_hat_Gamma_2[k];
      double gk0_part2 = pow(pst_cur[0], (eta_cur(k,0) == 0 ? 1.0 : 0.0)) * pow(pst_cur[2], (eta_cur(k,0) == 1 ? 1.0 : 0.0));
      double gk0 = exp(gk0_part1) * gk0_part2;
      double q2_post = gk1 / (gk1 + gk0);
      // Add bounds
      if (q2_post <= 0.1) q2_post = 0.0;
      if (q2_post >= 0.9) q2_post = 1.0;
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
    arma::mat Omega_hat_Gamma_1(K, K, fill::zeros);
    arma::mat Omega_hat_Gamma_2(K, K, fill::zeros);
    for (int ii = 0; ii < K; ii++) {
      Omega_hat_Gamma_1(ii, ii) = 1.0 / s2_hat_Gamma_1[ii];
      Omega_hat_Gamma_2(ii, ii) = 1.0 / s2_hat_Gamma_2[ii];
    }

    // Update alpha
    arma::vec U_alpha_1(K), U_alpha_2(K), W_alpha_1(K), W_alpha_2(K);
    for (int ii = 0; ii < K; ii++) {
      U_alpha_1[ii] = Gamma_hat_1[ii] - theta_cur(ii,0) - beta_cur[0] * gamma_cur(ii,0);
      U_alpha_2[ii] = Gamma_hat_2[ii] - theta_cur(ii,1) - beta_cur[1] * gamma_cur(ii,1);
      W_alpha_1[ii] = eta_cur(ii,0) * gamma_cur(ii,0);
      W_alpha_2[ii] = eta_cur(ii,1) * gamma_cur(ii,1);
    }

    if (sum(W_alpha_1) == 0) {
      alpha_cur[0] = 0.0;
    } else {
      arma::mat A_alpha_1 = W_alpha_1.t() * Omega_hat_Gamma_1 * W_alpha_1;
      arma::mat mu_alpha_1 = W_alpha_1.t() * Omega_hat_Gamma_1 * U_alpha_1 / A_alpha_1;
      alpha_cur[0] = R::rnorm(mu_alpha_1(0,0), sqrt(1.0 / A_alpha_1(0,0)));
    }

    if (sum(W_alpha_2) == 0) {
      alpha_cur[1] = 0.0;
    } else {
      arma::mat A_alpha_2 = W_alpha_2.t() * Omega_hat_Gamma_2 * W_alpha_2;
      arma::mat mu_alpha_2 = W_alpha_2.t() * Omega_hat_Gamma_2 * U_alpha_2 / A_alpha_2;
      alpha_cur[1] = R::rnorm(mu_alpha_2(0,0), sqrt(1.0 / A_alpha_2(0,0)));
    }
    alpha_1_tk[iter+1] = alpha_cur[0];
    alpha_2_tk[iter+1] = alpha_cur[1];

    // Update beta
    arma::vec U_beta_1(K), U_beta_2(K), W_beta_1(K), W_beta_2(K);
    for (int ii = 0; ii < K; ii++) {
      U_beta_1[ii] = Gamma_hat_1[ii] - theta_cur(ii,0) - alpha_cur[0] * eta_cur(ii,0) * gamma_cur(ii,0);
      U_beta_2[ii] = Gamma_hat_2[ii] - theta_cur(ii,1) - alpha_cur[1] * eta_cur(ii,1) * gamma_cur(ii,1);
      W_beta_1[ii] = gamma_cur(ii,0);
      W_beta_2[ii] = gamma_cur(ii,1);
    }
    arma::mat A_beta_1 = W_beta_1.t() * Omega_hat_Gamma_1 * W_beta_1;
    arma::mat A_beta_2 = W_beta_2.t() * Omega_hat_Gamma_2 * W_beta_2;
    arma::mat mu_beta_1 = W_beta_1.t() * Omega_hat_Gamma_1 * U_beta_1 / A_beta_1;
    arma::mat mu_beta_2 = W_beta_2.t() * Omega_hat_Gamma_2 * U_beta_2 / A_beta_2;
    beta_cur[0] = R::rnorm(mu_beta_1(0,0), sqrt(1.0 / A_beta_1(0,0)));
    beta_cur[1] = R::rnorm(mu_beta_2(0,0), sqrt(1.0 / A_beta_2(0,0)));
    beta_1_tk[iter+1] = beta_cur[0];
    beta_2_tk[iter+1] = beta_cur[1];

    // Store current values in trace
    for (int k = 0; k < K; k++) {
      theta_1_tk(iter+1, k) = theta_cur(k,0);
      theta_2_tk(iter+1, k) = theta_cur(k,1);
      gamma_1_tk(iter+1, k) = gamma_cur(k,0);
      gamma_2_tk(iter+1, k) = gamma_cur(k,1);
    }
  }

  List res = List::create(
    Named("beta_1") = beta_1_tk,
    Named("beta_2") = beta_2_tk,
    Named("pst_tk") = pst_tk,
    Named("eta_1_tk") = eta_1_tk,
    Named("eta_2_tk") = eta_2_tk,
    Named("alpha_1_tk") = alpha_1_tk,
    Named("alpha_2_tk") = alpha_2_tk
  );
  return res;
}
