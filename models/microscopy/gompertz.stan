//
// Gompertz function modelled with multplicative noise
//
data {
  int<lower=1> N_obs;  // number of obs
  int<lower=0> N_miss;  // number of missing obs
  int<lower=0> N_pred;  // number of predictions
  int<lower=1, upper=N_obs+N_miss+N_pred> y_train_idx[N_obs];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_missing_idx[N_miss];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_test_idx[N_pred];
  vector[N_obs] y_obs_train; // input data of OD values
  vector[N_pred] y_obs_test; // input data of OD test values
  vector[N_obs+N_miss+N_pred] time;  // time series
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N = N_obs+N_miss+N_pred;  // number of obs
  int<lower=0> N_log_lik = N*include_likelihood;  // sets the data size to 0 if the likelihood is not in the model
  int<lower=0> N_miss_l = N_miss*include_likelihood;
  int<lower=0> N_pred_l = N_pred*include_likelihood;
}

parameters {
  real<lower=0> tau;  // delay in hrs
  real<lower=0> L;  // carrying capacity
  real<lower=0> beta;  // growth parameter
  real<lower=0> f_0;  // initial condition of OD
  real<lower=0> sigma_meas;  // scale of noise
  vector<lower=0>[N_miss_l] y_missing;
  vector<lower=0>[N_pred_l] y_test;
}

transformed parameters {
  real<lower=1> K = pow(10, L);  // carrying capacity
  vector[N] log_f;
  vector[N] y;
  y[y_train_idx] = y_obs_train;
  if (include_likelihood) {
    y[y_missing_idx] = y_missing;
    y[y_test_idx] = y_test;
  }
  for (i in 1:N) log_f[i] = log(f_0) + log(K/f_0)*(1-exp(-(beta*fmax(time[i]-tau, 0.0))));  // log gompertz function
}

model {
  tau ~ gamma(5, 1);
  sigma_meas ~ std_normal();
  L ~ normal(4, 2);
  f_0 ~ std_normal();
  beta ~ std_normal();

  if (include_likelihood) y ~ lognormal(log_f, sigma_meas);
}

generated quantities {
  real y_tot[N];   // generated hyphal length values using samples from the model block
  real y_rep[N_obs];
  real y_pred[N_pred];
  vector[N_log_lik] log_lik;  // log likelihood
  vector[N] y_obs;
  y_obs[y_train_idx] = y_obs_train;
  y_obs[y_missing_idx] = y[y_missing_idx];
  y_obs[y_test_idx] = y_obs_test;

  y_tot = lognormal_rng(log_f, sigma_meas);
  // calculate log likelihood
  if (include_likelihood) {
    for (j in 1:N)
      log_lik[j] = lognormal_lpdf(y_obs[j] | log_f[j], sigma_meas);
  }
  y_rep = y_tot[y_train_idx];
  y_pred = y_tot[y_test_idx];
}
