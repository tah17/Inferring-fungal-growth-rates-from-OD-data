#
# Script that defines the Stan model meta data for the GP-OD-calibration model.
#
get_stan_data <- function(dat, run_likelihood = NULL){
  #
  # takes tibble of data and returns list that will be used as the Stan data
  # data must match the order of the data block in the corresponding model
  #
  if (is.null(run_likelihood)) {  #defaults to including the likelihood, and so the data, in the model.
    run_likelihood = 1
  }
  #
  # arranges data in format stan model is expecting
  #
  dat %>% 
    arrange(ID) -> df
  
  stan_data <- list(N_obs = length(which(df$training_data)),
                    N_miss = length(which(df$missing_data)),
                    N_pred = length(which(df$testing_data)),
                    N_c = length(which(df$blanks)),
                    N_exp = length(which(!df$blanks)),
                    c_idx = which(df$blanks),
                    exp_idx = which(!df$blanks),
                    y_train_idx = which(df$training_data),
                    y_missing_idx = which(df$missing_data),
                    y_test_idx = which(df$testing_data),
                    N_exp_obs = length(which(filter(df, !blanks)$training_data)),
                    N_exp_miss = length(which(filter(df, !blanks)$missing_data)),
                    N_exp_pred = length(which(filter(df, !blanks)$testing_data)),
                    R_t = length(unique(df$time)),
                    exp_train_idx = which(filter(df, !blanks)$training_data),
                    exp_missing_idx = which(filter(df, !blanks)$missing_data),
                    exp_test_idx = which(filter(df, !blanks)$testing_data),
                    x_idx = df$time + 1,
                    x = unique(df$time),
                    x_d = unique(df$time),
                    y_obs_train = filter(df, training_data)$OD,
                    y_obs_test = filter(df, testing_data)$OD,
                    mu_0 = mean(filter(df, training_data & blanks)$OD),
                    include_likelihood = run_likelihood)
  return(stan_data)
}
#
# arranges data in format stan model is expecting
#
get_processed_df <- function(dat){
  dat %>% 
    arrange(ID) -> df
  return(df)
}

iter <- 4000  #number of iterations the sampler will run for
no_of_chains <- 4  #number of chains the sampler will run for
params <- c("rho", "alpha", "sigma")  #parameters in the model - needs to be assigned as params
adapt_delta <- 0.99
#
## for high IC only
# iter <- 5000
# warmup <- 2000
# adapt_delta <- 0.99
#
# population parameters of the model
# to be called by tidybayes::gather_draws
#
pop_params2lang <- "cbind(sigma, rho, alpha, growth_rate, basal, delta_tilde)"
