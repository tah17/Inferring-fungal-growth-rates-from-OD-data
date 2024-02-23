//
// GP model for nuclear count data
//
functions {
  vector gp_pred_rng(real[] x2,  // prediction times
                     vector y1,  // observed nuclear count data
                     real[] x1,
                     real alpha,
                     real rho,
                     real sigma,
                     real delta) {

    //
    // Calculates the analytical GP posterior in the same way recommended by the Stan manual
    // see https://mc-stan.org/docs/stan-users-guide/fit-gp.html section "Analytical form of joint predictive inference"
    //
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] L_K;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;

      K = cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(square(sigma), N1));
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      f2_mu = (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      cov_f2 = cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }

  vector gp_deriv_rng(real[] x2,  // derivative times
                      vector y1,  // observed nuclear count data
                      real[] x1,  // observed times
                      real alpha,
                      real rho,
                      real sigma,
                      real delta) {
    //
    // Calculates the analytical GP posterior of derivative
    // where the predicted mean and covariance are detailed in https://doi.org/10.1038/ncomms13766
    //
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] g;
    {
      matrix[N1, N1] L_K;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N2, N2] k_x2_x2;
      matrix[N1, N2] v_pred;
      matrix[N2, N2] vT_pred_v_pred;
      vector[N2] g_mu;
      matrix[N2, N2] cov_g;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      matrix[N1, N2] dx2_k_x1_x2;

      K = cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(square(sigma), N1));
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      for (i in 1:N1) {
        for (j in 1:N2) {
          dx2_k_x1_x2[i,j] = (x1[i]-x2[j])*inv_square(rho)*k_x1_x2[i,j];
        }
      }
      g_mu = (dx2_k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, dx2_k_x1_x2);
      vT_pred_v_pred = v_pred' * v_pred;
      k_x2_x2 = cov_exp_quad(x2, alpha, rho);
      for (i in 1:N2) {
        for (j in 1:N2) {
          cov_g[i, j] = inv_square(rho)*k_x2_x2[i, j]*(1 - inv_square(rho)*square(x2[i] - x2[j])) - vT_pred_v_pred[i, j];
        }
      }
      diag_delta = diag_matrix(rep_vector(delta, N2));
      g = multi_normal_rng(g_mu, cov_g + diag_delta);
    }
    return g;
  }
}

data {
  int<lower=1> N_obs;  // number of obs
  int<lower=0> N_miss;  // number of missing obs
  int<lower=0> N_pred;  // number of predictions
  int<lower=1> R_t;  // number of hours experiment observed for
  int<lower=1, upper=N_obs+N_miss+N_pred> y_train_idx[N_obs];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_missing_idx[N_miss];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_test_idx[N_pred];
  int<lower=1, upper=R_t> x_idx[N_obs+N_miss+N_pred];  // indices of time covariate in dataframe
  real x[R_t];  // time covariate
  vector[N_obs] y_obs_train; // input data of nuclear count values
  vector[N_pred] y_obs_test; // input data of nuclear count test values
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N = N_obs+N_miss+N_pred;
  int<lower=0> N_log_lik = N*include_likelihood;  // transformed data that sets the data size to 0 if the likelihood is not in the model
  real<lower=0> delta=1e-9;
}

parameters {
  vector[R_t] f_tilde;
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[R_t, R_t] K = cov_exp_quad(x, alpha, rho) + diag_matrix(rep_vector(delta, R_t));
  matrix[R_t, R_t] L_K = cholesky_decompose(K);
  vector[R_t] f = L_K * f_tilde;
  vector[N_obs] log_y = log(y_obs_train);
}

model {
  rho ~ gamma(10, 1);
  alpha ~ normal(0, 0.5);
  sigma ~ std_normal();
  f_tilde ~ std_normal();
  if (include_likelihood) log_y ~ normal(f[x_idx[y_train_idx]], sigma);
}

generated quantities {
  real y_tot[N];
  real log_lik[N_log_lik];
  real y_rep[N_obs];
  real y_pred[N_pred];
  // vector[N_pred] f_pred;
  vector[R_t] g_pred;
  real growth_rate;

  y_tot[y_train_idx] = lognormal_rng(f[x_idx[y_train_idx]], sigma);
  y_tot[y_missing_idx] = lognormal_rng(f[x_idx[y_missing_idx]], sigma);
  if (include_likelihood) {
    for (i in 1:N_obs) log_lik[i] = normal_lpdf(log(y_obs_train[i]) | f[x_idx[y_train_idx[i]]], sigma);  // log likelihood
    for (i in 1:N_pred) log_lik[y_test_idx[i]] = normal_lpdf(log(y_obs_test[i])| f[x_idx[y_test_idx[i]]], sigma);

    // f_pred = gp_pred_rng(x[x_idx[y_test_idx]], log_y, x[x_idx[y_train_idx]], alpha, rho, sigma, delta);
    // y_tot[y_test_idx] = lognormal_rng(f_pred, sigma);

    g_pred = gp_deriv_rng(x, log_y, x[x_idx[y_train_idx]], alpha, rho, sigma, delta);  // predicts derivative
    growth_rate = max(g_pred);  // the growth rate is the maximum derivative of GP
  }
  y_rep = y_tot[y_train_idx];
  y_pred = y_tot[y_test_idx];
}
