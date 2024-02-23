//
// Latent GP model that models blank data and fungal growth in same file
//
functions {
  matrix generate_K_deriv(real[] x1, real[] x2, real alpha, real rho) {
    //
    // Generates kernel for GP and derivative,
    // where kernel can be seen in https://doi.org/10.1038/ncomms13766
    //
    int N1 = size(x1);
    int N2  = size(x2);
    matrix[N1+N2, N1+N2] K;
    matrix[N2, N2] k_x2_x2 = cov_exp_quad(x2, alpha, rho);
    K[1:N1, 1:N1] = cov_exp_quad(x1, alpha, rho); // exp. quadratic kernel evaluated at observed times x1: K(x1, x1)
    K[1:N1, N1+1:N1+N2] = (1/square(rho))*(rep_matrix(to_vector(x1), N2) - rep_matrix(to_row_vector(x2), N1)).*cov_exp_quad(x1, x2, alpha, rho);  // Partial derivative of K(x1, x2) w.r.t. x2: delta_2(K(x1, x2))
    K[N1+1:N1+N2, 1:N1] = K[1:N1, N1+1:N1+N2]';  // delta_1(K(x2, x1))
    for (i in 1:N2)
      for (j in 1:N2)
        K[N1+i, N1+j] = (k_x2_x2[i, j]*inv_square(rho))*(1 - (square(x2[i] - x2[j])*inv_square(rho)));  // delta_1(delta_2(K(x2, x2)))
    return K;
  }
}

data {
  int<lower=1> N_obs;  // number of obs
  int<lower=0> N_miss;  // number of missing obs
  int<lower=0> N_pred;  // number of predictions
  int<lower=0> N_c;  // number of blanks
  int<lower=0> N_exp;  // number of exps (not blanks)
  int<lower=0> c_idx[N_c];  // blank indices
  int<lower=0> exp_idx[N_exp];
  int<lower=1> N_exp_obs;  // number of obs
  int<lower=0> N_exp_miss;  // number of missing obs
  int<lower=0> N_exp_pred;  // number of predictions
  int<lower=1> R_t; // number of hours experiment observed for
  int<lower=1, upper=N_obs+N_miss+N_pred> y_train_idx[N_obs];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_missing_idx[N_miss];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_test_idx[N_pred];
  int<lower=1, upper=N_exp_obs+N_exp_miss+N_exp_pred> exp_train_idx[N_exp_obs];
  int<lower=1, upper=N_exp_obs+N_exp_miss+N_exp_pred> exp_missing_idx[N_exp_miss];
  int<lower=1, upper=N_exp_obs+N_exp_miss+N_exp_pred> exp_test_idx[N_exp_pred];
  int<lower=1, upper=R_t> x_idx[N_exp+N_c];
  real x[R_t];   // unique time series
  real x_d[R_t];  // unique times that derivative is wanted for
  vector[N_obs] y_obs_train; // input data of OD values
  vector[N_pred] y_obs_test; // input data of OD test values
  real mu_0;  // sample mean
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N = N_obs+N_miss+N_pred;
  int<lower=0> N_log_lik = N*include_likelihood;  // transformed data that sets the data size to 0 if the likelihood is not in the model
  real<lower=0> delta=1e-9;  // small offset to ensure matrices are pos. def.
  int<lower=0> N_miss_l = N_miss*include_likelihood;
  int<lower=0> N_pred_l = N_pred*include_likelihood;
  real<lower=0> sol_volume = 200;  // volume init inocula of spores are in
  real<lower=0> max_ic = 2e5/sol_volume; // max init. conds scaled to [N/ul]
}

parameters {
  vector[2*R_t] f_g_tilde; // GP for fungus with derivative
  real<lower=0> rho;  // lengthscale (GP hyperparam)
  real<lower=0> alpha;  // marginal var (GP hyperparam)
  real<lower=0> sigma; // obs. scale of noise
  real<lower=0> basal;  // linear transform parameter (offset)
  real<lower=0> delta_tilde;  // linear transform parameter (slope)
  vector[N_miss_l] log_y_missing;
  vector[N_pred_l] log_y_test;
}

transformed parameters {
  real<lower=1> delta_prop = pow(10, delta_tilde);
  vector[2*R_t] f_g;  // vector where 1:R_t is GP for logged fungal growth and R_t+1:2*R_t is the derivative
  vector[N] log_y;
  log_y[y_train_idx] = log(y_obs_train);
  if (include_likelihood) {
    log_y[y_missing_idx] = log_y_missing;
    log_y[y_test_idx] = log_y_test;
  }
  {
    matrix[2*R_t, 2*R_t] K = generate_K_deriv(x, x_d, alpha, rho) + diag_matrix(rep_vector(delta, 2*R_t));   // uses kernel for GP + derivative (generate_K_deriv)
    matrix[2*R_t, 2*R_t] L_K = cholesky_decompose(K);
    f_g = L_K * f_g_tilde;
  }
}

model {
  rho ~ gamma(10, 1);
  alpha ~ normal(5, 2);
  sigma ~ normal(0, 0.5);
  f_g_tilde ~ std_normal();
  basal ~ lognormal(log(mu_0), 1);
  delta_tilde ~ cauchy(log10(max_ic), 1);

  if (include_likelihood) {
    log_y[exp_idx] ~ normal(log(basal + exp(f_g[x_idx[exp_idx]])/delta_prop), sigma);  // exp transform of GP needed as GP is for logged fungal growth
    log_y[c_idx] ~ normal(log(basal), sigma);   // for blank data
  }
}

generated quantities {
  real y_tot[N];
  real log_lik[N_log_lik];
  real growth_rate;
  vector[N] log_y_obs;
  log_y_obs[y_train_idx] = log(y_obs_train);
  log_y_obs[y_missing_idx] = log_y[y_missing_idx];
  log_y_obs[y_test_idx] = log(y_obs_test);

  y_tot[exp_idx] = lognormal_rng(log(basal + exp(f_g[x_idx[exp_idx]])/delta_prop), sigma);
  for (i in c_idx) y_tot[i] = lognormal_rng(log(basal), sigma);
  growth_rate = max(f_g[R_t+1:(2*R_t)]);  // growth rate is max of the derivative (derivative values are from R_t+1:(2*R_t))
  if (include_likelihood) {  // calc. log likelihood
    for (i in c_idx) log_lik[i] = normal_lpdf(log_y_obs[i] | log(basal), sigma);
    for (i in exp_idx) log_lik[i] = normal_lpdf(log_y_obs[i] | log(basal + exp(f_g[x_idx[i]])/delta_prop), sigma);
  }
}
