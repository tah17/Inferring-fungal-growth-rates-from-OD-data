#
# Script that defines the Stan model meta data for the gompertz model and function which creates a Stan data list.
#
get_stan_data <- function(dat, run_likelihood = NULL){
  #
  # takes tibble of data and returns list that will be used as the Stan data
  # data must match the order of the data block in the corresponding model
  #
  if (is.null(run_likelihood)) {  #defaults to including the likelihood, and so the data, in the model.
    run_likelihood = 1
  }

  stan_data <- list(N_obs = length(which(dat$training_data)),
                    N_miss = length(which(dat$missing_data)),
                    N_pred = length(which(dat$testing_data)),
                    y_train_idx = which(dat$training_data),
                    y_missing_idx = which(dat$missing_data),
                    y_test_idx = which(dat$testing_data),
                    y_obs_train = filter(dat, training_data)$length,
                    y_obs_test = filter(dat, testing_data)$length,
                    time = dat$time,
                    include_likelihood = run_likelihood)
  return(stan_data)
}
#
# empty function because data was not changed when fed into the stan model
#
get_processed_df <- function(dat){
  return(dat)
}

iter <- 2000  #number of iterations the sampler will run for
no_of_chains <- 4  #number of chains the sampler will run for
params <- c("f_0", "beta", "tau", "sigma_meas", "L")  #parameters in the model - needs to be assigned as params
#
# population parameters of the model
# to be called by tidybayes::gather_draws
#
pop_params2lang <- "cbind(f_0, beta, tau, sigma_meas, L)"
