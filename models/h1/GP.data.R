#
# Script that defines the Stan model meta data for a simple GP.
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
  dat %>% arrange(ID) -> df
  
  stan_data <- list(N_obs = length(which(df$training_data)),
                    N_miss = length(which(df$missing_data)),
                    N_pred = length(which(df$testing_data)),
                    R_t = length(unique(df$time)),
                    y_train_idx = which(df$training_data),
                    y_missing_idx = which(df$missing_data),
                    y_test_idx = which(df$testing_data),
                    x_idx = df$time + 1,
                    x = unique(df$time),
                    y_obs_train = filter(df, training_data)$nuclei,
                    y_obs_test = filter(df, testing_data)$nuclei,
                    include_likelihood = run_likelihood)
  return(stan_data)
}
#
# function that arranges data to the same format that stan model received the data 
#
get_processed_df <- function(dat){
  dat %>% arrange(ID) -> df
  return(df)
}

iter <- 2000  #number of iterations the sampler will run for
no_of_chains <- 4  #number of chains the sampler will run for
params <- c("rho", "alpha", "sigma")  #parameters in the model - needs to be assigned as params
#
# population parameters of the model
# to be called by tidybayes::gather_draws
#
pop_params2lang <- "cbind(rho, alpha, sigma)"