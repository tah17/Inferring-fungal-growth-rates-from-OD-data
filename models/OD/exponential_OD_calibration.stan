//
// Exponential function modelled with multplicative noise
//
data {
  int<lower=1> N_obs;  // number of obs
  int<lower=0> N_miss;  // number of missing obs
  int<lower=0> N_pred;  // number of predictions
  int<lower=0> N_c;  // number of blanks
  int<lower=0> N_exp;  // number of non blanks "exps"
  int<lower=0> c_idx[N_c];   // blank indices
  int<lower=0> exp_idx[N_exp];  // exp indices
  int<lower=1, upper=N_obs+N_miss+N_pred> y_train_idx[N_obs];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_missing_idx[N_miss];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_test_idx[N_pred];
  vector[N_obs+N_miss+N_pred] IC;  // init. conds
  vector[N_obs] y_obs_train; // input data of OD train values
  vector[N_pred] y_obs_test; // input data of OD test values
  vector[N_obs+N_miss+N_pred] time;  // time series
  real mu_0;  // mean value of blanks
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N = N_obs+N_miss+N_pred;
  int<lower=0> N_log_lik = N*include_likelihood;  // sets the data size to 0 if the likelihood is not in the model
  int<lower=0> N_miss_l = N_miss*include_likelihood;
  int<lower=0> N_pred_l = N_pred*include_likelihood;
  real<lower=0> sol_volume = 200;  // volume init inocula of spores are in
  vector[N] IC_scaled = IC ./ sol_volume;  // init. conds scaled (so that f [N/ul])
  vector[N] log_IC = log(IC_scaled);
}

parameters {
  real<lower=0> tau;  // delay before growth
  real<lower=0> basal;  // linear transform parameter
  real<lower=0> beta;  // growth parameter
  real<lower=log10(min(IC))> delta_tilde;  // linear transform parameter
  real<lower=0> sigma_meas;  // scale of the obs noise
  vector<lower=0>[N_miss_l] y_missing;
  vector<lower=0>[N_pred_l] y_test;
}

transformed parameters {
  real<lower=min(IC_scaled)> delta = pow(10, delta_tilde);
  vector[N_exp] log_z;  // latent log exponential growth
  vector[N] log_y;
  log_y[y_train_idx] = log(y_obs_train);
  if (include_likelihood) {
    log_y[y_missing_idx] = y_missing;
    log_y[y_test_idx] = y_test;
  }
  for (i in 1:N_exp) log_z[i] = log_IC[exp_idx[i]] + beta*fmax(time[exp_idx[i]]-tau, 0.0);
}

model {
  tau ~ gamma(5, 1);
  basal ~ lognormal(log(mu_0), 1);
  delta_tilde ~ cauchy(log10(max(IC_scaled)), 1);
  sigma_meas ~ normal(0, 0.5);
  beta ~ std_normal();
  if (include_likelihood) {
    log_y[exp_idx] ~ normal(log(basal + exp(log_z)/delta), sigma_meas);
    log_y[c_idx] ~ normal(log(basal), sigma_meas);
  }
}

generated quantities {
  real y_tot[N];
  real y_rep[N_obs];
  real y_pred[N_pred];
  vector[N_log_lik] log_lik;  // log likelihood
  vector[N] log_y_obs;
  log_y_obs[y_train_idx] = log(y_obs_train);
  log_y_obs[y_missing_idx] = log_y[y_missing_idx];
  log_y_obs[y_test_idx] = log(y_obs_test);

  y_tot[exp_idx] = lognormal_rng(log(basal + exp(log_z)/delta), sigma_meas);
  for (i in c_idx)
    y_tot[i] = lognormal_rng(log(basal), sigma_meas);
  // calculate log likelihood
  if (include_likelihood) {
    for (i in c_idx)
      log_lik[i] = normal_lpdf(log_y_obs[i] | log(basal), sigma_meas);
    for (j in 1:N_exp)
      log_lik[exp_idx[j]] = normal_lpdf(log_y_obs[exp_idx[j]] | log(basal + exp(log_z[j])/delta), sigma_meas);
  }
  y_rep = y_tot[y_train_idx];
  y_pred = y_tot[y_test_idx];
}
