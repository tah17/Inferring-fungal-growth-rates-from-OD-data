#
# Script that defines the Stan model meta data for the exponential model and function which creates a Stan data list.
#
get_stan_data <- function(dat, run_likelihood = NULL){
  #
  # takes tibble of data and returns list that will be used as the Stan data
  # data must match the order of the data block in the corresponding model
  #
  if (is.null(run_likelihood)) {  #defaults to including the likelihood, and so the data, in the model.
    run_likelihood = 1
  }

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
                    y_obs_train = filter(df, training_data)$OD,
                    y_obs_test = filter(df, testing_data)$OD,
                    time = df$time,
                    mu_0 = mean(filter(df, training_data & blanks)$OD),
                    include_likelihood = run_likelihood)
  return(stan_data)
}

get_processed_df <- function(dat){
  dat %>% arrange(ID) -> df
  return(df)
}

iter <- 2000  #number of iterations the sampler will run for
no_of_chains <- 4  #number of chains the sampler will run for
params <- c("tau", "basal", "f_0", "beta", "sigma_meas")  #parameters in the model - needs to be assigned as params
#
# population parameters of the model
# to be called by tidybayes::gather_draws
#
pop_params2lang <- "cbind(basal, tau, f_0, beta, sigma_meas)"
