//
// Defines a simple Random Walk (RW) model, with a vector of values 'y' modeled as a RW
// with standard deviation 'sigma'.
//
data {
  int<lower=1> N_obs;  // number of obs
  int<lower=0> N_miss;  // number of missing obs
  int<lower=0> N_pred;  // number of predictions
  int<lower=1, upper=N_obs+N_miss+N_pred> y_train_idx[N_obs];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_missing_idx[N_miss];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_test_idx[N_pred];
  int<lower=0> time[N_obs+N_miss+N_pred]; // time of OD readings
  vector[N_obs] y_obs_train; // input data of OD values
  vector[N_pred] y_obs_test; // input data of OD values for testing - to be used in calculating LPD
  int<lower=0, upper=1> include_likelihood; // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N = N_obs+N_miss+N_pred;
  int<lower=0> N_log_lik = N*include_likelihood;  // transformed data that sets the data size to 0 if the likelihood is not in the model
  int<lower=1> R_t = max(time)-min(time)+1;  // number of hours experiment observed for
  int<lower=1> R = N/R_t;  // tot number of obs reps
  int<lower=0,upper=N> start_idxs[R];
  for (i in 1:R) start_idxs[i] = (R_t*(i-1))+1;
}

parameters {
  real<lower=0> sigma;  // s.d. of the RW
  vector[N_miss] y_missing;
  vector[N_pred] y_test;
}

transformed parameters {
  vector[N] y;
  y[y_train_idx] = y_obs_train;
  y[y_missing_idx] = y_missing;
  y[y_test_idx] = y_test;
}

model {
  sigma ~ normal(0, 0.5);
  if (include_likelihood) {
    for (i in 1:R) y[((R_t*(i-1))+2):(R_t*i)] ~ normal(y[((R_t*(i-1))+1):((R_t*i)-1)], sigma);
  }
}

generated quantities {
  real y_tot[N];
  real log_lik[N_log_lik];
  real y_rep[N_obs];
  real y_pred[N_pred];
  vector[N] y_obs;
  y_obs[y_train_idx] = y_obs_train;
  y_obs[y_missing_idx] = y[y_missing_idx];
  y_obs[y_test_idx] = y_obs_test;

  for (i in 1:R) {
    y_tot[start_idxs[i]] = y[start_idxs[i]];
    if (include_likelihood) log_lik[start_idxs[i]] = normal_lpdf(y[start_idxs[i]] | y[start_idxs[i]], sigma);  // log likelihood
    for (j in ((R_t*(i-1))+1):((R_t*i)-1)) {
      y_tot[j+1] = normal_rng(y[j], sigma);   //change to y_tot for fake data simulation
      if (include_likelihood) log_lik[j+1] = normal_lpdf(y_obs[j+1] | y[j], sigma);
    }
  }
  y_rep = y_tot[y_train_idx];
  y_pred = y_tot[y_test_idx];
}
