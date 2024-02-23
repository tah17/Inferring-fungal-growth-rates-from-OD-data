rm(list = ls())
library(rstan)
library(reshape2)
library(tidyverse)
library(forecast)
seed <- 683861
set.seed(seed)
source("R/03_validate/functions.R")
model_list <- c("RW", "mixed_logistic_OD_calibration", "exponential", "gompertz", "gompertz_OD_calibration", "GP", "logistic", "exponential_OD_calibration", "logistic_OD_calibration_no_delay", "logistic_OD_calibration", "GP_OD_calibration")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
#
# Read in command line arguments
#
myargs = commandArgs(trailingOnly=TRUE)
IC_used <- myargs[1] # ICs used - all or high
if (length(myargs) < 1) { 
  stop("Error: missing command line argument 'IC_used'.")
} else if (length(myargs) < 2) {
  fold_no <- 5
  warning("Using default Number of Folds")
} else {
  fold_no <- as.integer(myargs[2])
}
#
# Read in data set
#
df <- readRDS(file=paste("data/OD_", IC_used, ".Rda", sep=""))
cv_idxs <- readRDS(file = paste("data/split_idxs_", fold_no, "_fold_", IC_used, ".Rda", sep=""))  #Gets cross validation indices
job_idx_list <- expand.grid(Fold = 1:length(cv_idxs$train_idxs), Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
#
# Get Job Idx
#
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
#
# Specify Stan Model
#
model_idx <- job_idx_list$Model[job_idx]
model <- model_list[model_idx]
stanfile <- paste('models/OD/', model, '.stan', sep = "")
stanmodel <- stan_model(stanfile)
source(paste('models/OD/', model, '.data.R', sep = ""))
#
# Read in data of job
#
chosen_fold <- job_idx_list$Fold[job_idx]
chosen_chain <- job_idx_list$Chains[job_idx]
fold_df <- get_labelled_data(chosen_fold, df, cv_idxs)
# Fit Model --------------------------------------------------
if (exists("warmup")) {
  stan_fit <- sampling(stanmodel, data = get_stan_data(fold_df), chains = 1, warmup = warmup, iter = iter, seed = seed, chain_id=chosen_chain, cores = 1)
} else if (exists("adapt_delta")) {
  stan_fit <- sampling(stanmodel, data = get_stan_data(fold_df), chains = 1, iter = iter, seed = seed, chain_id=chosen_chain, cores = 1, control=list(adapt_delta=adapt_delta))
} else {
  stan_fit <- sampling(stanmodel, data = get_stan_data(fold_df), chains = 1, iter = iter, seed = seed, chain_id=chosen_chain, cores = 1)
}
if (!dir.exists(file.path(storage_loc, paste("/OD_output/", IC_used, sep="")))) {
  dir.create(file.path(storage_loc, paste("/OD_output/", IC_used, sep="")))
}
saveRDS(stan_fit, file=paste(storage_loc, "/OD_output/", IC_used, "/", fold_no, "_fold_", model, "_", job_idx, "_fit.Rda", sep=""))  #stores model
