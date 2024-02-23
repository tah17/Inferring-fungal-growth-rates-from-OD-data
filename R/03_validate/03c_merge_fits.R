rm(list = ls())
library(tidyverse)
library(rstan)
library(forecast)
library(data.table)
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
IC_used <-  myargs[1]  #ICs used - all or high
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
no_of_folds <- length(cv_idxs$train_idxs)
full_job_idx_list <- expand.grid(Fold = 1:length(cv_idxs$train_idxs), Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets prev up job list used in 03b_
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
full_job_idx_list %>%
  filter(Model == job_idx) %>%
  select (-Model) -> job_idx_list
#
# get model
#
model <- model_list[job_idx]
source(paste('models/OD/', model, '.data.R', sep = ""))
#
# read in model training fits
#
fit_idx <- which(full_job_idx_list$Model==job_idx)
stan_list <- lapply(fit_idx, function(x) readRDS(file = paste(storage_loc, "/OD_output/", IC_used, "/", fold_no, "_fold_", model, "_", x, "_fit.Rda", sep="")))
#
# merge different chains from same fold
#
stan_kfold <- lapply(1:no_of_folds, function(x) sflist2stanfit(stan_list[job_idx_list$Fold==x])) # combine stan fits from chains but using the same cv fold (aka fitting to the same data) into a single stan fit using sflist2stanfit
#
# save merged stan fit objects
#
lapply(1:no_of_folds, function(i) saveRDS(stan_kfold[[i]], paste(storage_loc, '/OD_output/', IC_used, "/", fold_no, "_fold_", model, '_train', i, '.rds', sep = "")))  #saves fits
#
# calculate training + testing errors and save
#
dat <- lapply(1:no_of_folds, function(x) get_labelled_data(x, df, cv_idxs))
lapply(1:no_of_folds, function(x) get_metrics(get_processed_df(dat[[x]]), stan_kfold[[x]]) %>% mutate(Fold = paste("Fold", x, sep=""))) %>%
  rbindlist() -> fit_metrics
saveRDS(fit_metrics, paste('output/', IC_used, "/", fold_no, "_fold_", model, '_fit_stats.rds', sep = ""))
