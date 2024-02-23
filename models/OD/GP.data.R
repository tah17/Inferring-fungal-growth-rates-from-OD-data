#
# Script that defines the Stan model meta data for the GP.
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
    arrange(ID) %>% 
    mutate(groups=case_when(blanks ~ 2, !blanks ~ 1)) %>%  ## arrange by blanks or fungus so they can be indexed in Stan file (we fit a GP to each for fair prediction comparison)
    arrange(groups) -> df
  
  stan_data <- list(G = length(unique(df$groups)),
                    N_obs = length(which(df$training_data)),
                    N_miss = length(which(df$missing_data)),
                    N_pred = length(which(df$testing_data)),
                    R_t = length(unique(df$time)),
                    y_train_idx = which(df$training_data),
                    y_missing_idx = which(df$missing_data),
                    y_test_idx = which(df$testing_data),
                    x_idx = df$time + 1,
                    y_g_idx = as.numeric(df$groups),
                    x = unique(df$time),
                    y_obs_train = filter(df, training_data)$OD,
                    y_obs_test = filter(df, testing_data)$OD,
                    s_train = unlist(lapply(unique(df$groups), function(x) length(filter(filter(df, training_data), groups==x)$OD))),
                    s_test = unlist(lapply(unique(df$groups), function(x) length(filter(filter(df, testing_data), groups==x)$OD))),
                    include_likelihood = run_likelihood)
  return(stan_data)
}
#
# arranges data in format stan model got data
#
get_processed_df <- function(dat){
  dat %>% 
    arrange(ID) %>% 
    mutate(groups=case_when(blanks ~ 2, !blanks ~ 1)) %>%
    arrange(groups) -> df
  return(df)
}

iter <- 2000  #number of iterations the sampler will run for
no_of_chains <- 4  #number of chains the sampler will run for
params <- c("rho", "alpha", "sigma")  #parameters in the model - needs to be assigned as params
#
# population parameters of the model
# to be called by tidybayes::gather_draws
#
pop_params2lang <- "cbind(sigma)"
ind_params2lang <- "cbind(rho, alpha, growth_rate_w_control)[rep]"  #individual level parameters of the model
