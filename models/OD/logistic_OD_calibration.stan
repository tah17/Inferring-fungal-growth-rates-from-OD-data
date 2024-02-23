//
// Latent logistic function with delay in beginning growth modelled with multiplicative noise
//
data {
  int<lower=1> N_obs;  // number of obs
  int<lower=0> N_miss;  // number of missing obs
  int<lower=0> N_pred;  // number of predictions
  int<lower=0> N_c;  // number of blanks
  int<lower=0> N_exp;  // number of exps (not blanks)
  int<lower=0> c_idx[N_c];  // blank indices
  int<lower=0> exp_idx[N_exp];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_train_idx[N_obs];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_missing_idx[N_miss];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_test_idx[N_pred];
  vector[N_obs+N_miss+N_pred] IC;  // init. conds
  vector[N_obs] y_obs_train; // input data of OD values
  vector[N_pred] y_obs_test; // input data of OD test values
  vector[N_obs+N_miss+N_pred] time;  // time series
  real mu_0;  // sample mean of blanks
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N = N_obs+N_miss+N_pred;
  int<lower=0> N_log_lik = N*include_likelihood;  // sets the data size to 0 if the likelihood is not in the model
  int<lower=0> N_miss_l = N_miss*include_likelihood;
  int<lower=0> N_pred_l = N_pred*include_likelihood;
  real<lower=0> sol_volume = 200;  // volume init inocula of spores are in
  vector[N] IC_scaled = IC ./ sol_volume;  // init. conds scaled (so that f [N/ul])
}

parameters {
  real<lower=0> tau;  // growth delay
  real<lower=log10(max(IC_scaled))> L;
  real<lower=log10(min(IC_scaled))> delta_tilde;
  real<lower=0> beta;  // hyphal growth rate
  real<lower=0> basal;  // linear transform parameter (offset)
  real<lower=0> sigma_meas;  // scale of observed noise
  vector<lower=0>[N_miss_l] y_missing;
  vector<lower=0>[N_pred_l] y_test;
}

transformed parameters {
  vector[N_exp] f;
  real<lower=max(IC_scaled)> K = pow(10, L);  // carrying capacity
  real<lower=min(IC_scaled)> delta = pow(10, delta_tilde);  // linear transform parameter
  vector[N] y;
  y[y_train_idx] = y_obs_train;
  if (include_likelihood) {
    y[y_missing_idx] = y_missing;
    y[y_test_idx] = y_test;
  }
  for (i in 1:N_exp) {
    f[i] = K/(1 + ((K - IC_scaled[exp_idx[i]])/IC_scaled[exp_idx[i]])*exp(-(beta)*fmax(time[exp_idx[i]] - tau, 0.0)));  // logistic function and fmax ensures function starts at t=tau
  }
}

model {
  tau ~ gamma(5, 1);
  sigma_meas ~ normal(0, 0.5);
  L ~ normal(9, 2);
  beta ~ std_normal();
  // beta ~ cauchy(0, 1);  // for prior sensitivity
  basal ~ lognormal(log(mu_0), 1);
  delta_tilde ~ cauchy(log10(max(IC_scaled)), 1);

  if (include_likelihood) {
    y[exp_idx] ~ lognormal(log(basal + f/delta), sigma_meas);
    y[c_idx] ~ lognormal(log(basal), sigma_meas);
  }
}

generated quantities {
  real y_tot[N];
  real y_rep[N_obs];
  real y_pred[N_pred];
  vector[N_log_lik] log_lik;  // log likelihood
  vector[N] y_obs;
  y_obs[y_train_idx] = y_obs_train;
  y_obs[y_missing_idx] = y[y_missing_idx];
  y_obs[y_test_idx] = y_obs_test;

  y_tot[exp_idx] = lognormal_rng(log(basal + f/delta), sigma_meas);
  for (i in c_idx)
    y_tot[i] = lognormal_rng(log(basal), sigma_meas);
  // calculate log likelihood
  if (include_likelihood) {
    for (i in c_idx)
      log_lik[i] = lognormal_lpdf(y_obs[i] | log(basal), sigma_meas);
    for (j in 1:N_exp)
      log_lik[exp_idx[j]] = lognormal_lpdf(y_obs[exp_idx[j]] | log(basal + f[j]/delta), sigma_meas);
  }
  y_rep = y_tot[y_train_idx];
  y_pred = y_tot[y_test_idx];
}
