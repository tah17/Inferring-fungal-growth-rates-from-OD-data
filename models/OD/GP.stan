//
// GP model that models blank data and OD of fungal growth in same file
//
functions {
  vector gp_pred_rng(real[] x2,  // prediction times
                     vector y1,  // observed OD
                     real[] x1,  // observed time
                     real alpha, // hyperparameters for cov_exp_quad kernel
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
      matrix[N1, N1] K;
      matrix[N1, N1] L_K;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;

      K = cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(square(sigma), N1)); // kernel K evaluated at the observed times x1 + observed noise of scale sigma
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);  // kernel K evaluated at observed times x1 and predicted times x2, K(x1, x2)
      f2_mu = (k_x1_x2' * K_div_y1);  
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      cov_f2 = cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred;  
      diag_delta = diag_matrix(rep_vector(delta, N2)); 

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
  vector gp_deriv_rng(real[] x2,  // derivative times 
                      vector y1,  // observed OD
                      real[] x1,  // observed time
                      real alpha,
                      real rho,
                      real sigma,
                      real delta) {
    //
    // Calculates the analytical GP posterior of derivative
    // where mean and covariance of derivative are detailed in https://doi.org/10.1038/ncomms13766
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
          dx2_k_x1_x2[i,j] = (x1[i]-x2[j])*inv_square(rho)*k_x1_x2[i,j];  // partial derivative of K(x1, x2) w.r.t to variable x2, delta_2(K(x1, x2))
        }
      }
      g_mu = (dx2_k_x1_x2' * K_div_y1);   // transpose of delta_2(K(x1, x2)) is partial derivative of K(x2, x1) w.r.t to variable x2, delta_1(K(x2, x1)).
      v_pred = mdivide_left_tri_low(L_K, dx2_k_x1_x2);
      vT_pred_v_pred = v_pred' * v_pred;
      k_x2_x2 = cov_exp_quad(x2, alpha, rho);
      for (i in 1:N2) {
        for (j in 1:N2) {
          // delta_1(delta_2(K(x2, x2))) - ...
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
  int<lower=0> G;   // # of groups which is 2, one for fungus and one for blanks (needs to be in this order)
  int<lower=1> N_obs;  // number of obs
  int<lower=0> N_miss;  // number of missing obs
  int<lower=0> N_pred;  // number of predictions
  int<lower=1> R_t; // number of hours experiment observed for
  int<lower=1, upper=N_obs+N_miss+N_pred> y_train_idx[N_obs];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_missing_idx[N_miss];
  int<lower=1, upper=N_obs+N_miss+N_pred> y_test_idx[N_pred];
  int<lower=1, upper=R_t> x_idx[N_obs+N_miss+N_pred];  // time index
  int<lower=1, upper=G> y_g_idx[N_obs+N_miss+N_pred];  // inidices of which group the data belongs to
  real x[R_t];  // unique time course
  vector[N_obs] y_obs_train; // input data of OD values
  vector[N_pred] y_obs_test; // input data of OD test values
  int s_train[G];   // list of the number of training points in each group
  int s_test[G];   // list of the number of testing points in each group
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N = N_obs+N_miss+N_pred;
  int<lower=0> N_log_lik = N*include_likelihood;  // transformed data that sets the data size to 0 if the likelihood is not in the model
  real<lower=0> delta=1e-9;  // small offset to ensure matrices are pos. def
}

parameters {
  matrix[R_t, G] f_tilde; // GP for blank and GP for fungus
  vector<lower=0>[G] rho;  // Hyperparams of GP kernels for each GP, rho is lengthscale and alpha is the marginal var. 
  vector<lower=0>[G] alpha;  
  real<lower=0> sigma;  // scale of observed noise on log of OD
}

transformed parameters {
  vector[N_obs] log_y = log(y_obs_train);
  matrix[R_t, G] f;
  for (g in 1:G) {
    // GP with exponentiated quadratic kernel
    matrix[R_t, R_t] K = cov_exp_quad(x, alpha[g], rho[g]) + diag_matrix(rep_vector(delta, R_t));
    matrix[R_t, R_t] L_K = cholesky_decompose(K);  // transformation w/ cholesky decomp. for efficiency as recommended by Stan manual.
    f[, g] = L_K * f_tilde[, g];
  }
}

model {
  rho ~ gamma(10, 1);
  alpha ~ normal(0, 0.5);
  sigma ~ normal(0, 0.5);
  for (g in 1:G) f_tilde[, g] ~ std_normal();
  if (include_likelihood) {
    for (i in 1:N_obs) log_y[i] ~ normal(f[x_idx[y_train_idx[i]], y_g_idx[y_train_idx[i]]], sigma);
  }
}

generated quantities {
  real y_tot[N];
  real log_lik[N_log_lik];
  real y_rep[N_obs];
  real y_pred[N_pred];
  vector[N_pred] f_pred;
  matrix[R_t, G] g_pred;
  vector[G] growth_rate_w_control; // growth rate and growth rate of the blank data
  vector[G-1] growth_rate;  // growth rate of interest only
  int pos1;
  int pos2;
  pos1 = 1;  // init. counter of training data points 
  pos2 = 1;  // init. counter of testing data points

  for (i in y_train_idx) y_tot[i] = lognormal_rng(f[x_idx[i], y_g_idx[i]], sigma);
  for (i in y_missing_idx) y_tot[i] = lognormal_rng(f[x_idx[i], y_g_idx[i]], sigma);
  if (include_likelihood) {
    for (i in 1:N_obs) // log likelihood
      log_lik[y_train_idx[i]] = normal_lpdf(log(y_obs_train[i]) | f[x_idx[y_train_idx[i]], y_g_idx[y_train_idx[i]]], sigma);  
    for (g in 1:G) {
      //
      // checks that the test set of a group doesn't have a size of 0 and 
      // then predicts an output at the prediction times segmented by
      // number of test points in group g (x[x_idx[y_test_idx[pos2:(pos2+s_test[g]-1)]]]) 
      // using a GP trained on training points in group g (segment(log_y, pos1, s_train[g])).
      //
      if (s_test[g]!=0) f_pred[pos2:(pos2+s_test[g]-1)] = gp_pred_rng(x[x_idx[y_test_idx[pos2:(pos2+s_test[g]-1)]]], segment(log_y, pos1, s_train[g]), x[x_idx[segment(y_train_idx, pos1, s_train[g])]], alpha[g], rho[g], sigma, delta);
      //
      // predicts the derivative for each group at all the training points in that group g
      //
      g_pred[, g] = gp_deriv_rng(x, segment(log_y, pos1, s_train[g]), x[x_idx[segment(y_train_idx, pos1, s_train[g])]], alpha[g], rho[g], sigma, delta);
      growth_rate_w_control[g] = max(g_pred[, g]);  // the growth rate is the maximum derivative of GP 
      pos1 = pos1 + s_train[g]; // adds to the training points counter
      pos2 = pos2 + s_test[g];  // adds to the testing points counter
    }
    growth_rate = growth_rate_w_control[1:(G-1)];  // drop the growth rate calculated for the blanks
    for (i in 1:N_pred) {
      log_lik[y_test_idx[i]] = normal_lpdf(log(y_obs_test[i])| f_pred[i], sigma);  // log likelihood 
      y_tot[y_test_idx[i]] = lognormal_rng(f_pred[i], sigma);
    }
  }
  y_rep = y_tot[y_train_idx];
  y_pred = y_tot[y_test_idx];
}
