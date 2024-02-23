library(tidyverse)
library(loo)
#' get_data_from_idxs
#'
#' Takes a tibble where the data is assigned a CV training fold and returns a stan data list, by calling the user defined get_stan_data.
#'
#' @param i Int. The cross validation fold label.
#' @param df Tibble. Tibble that contains the original data set.
#' @param cv_idxs List. List of training (and maybe testing) indices of the data set.
#'
#' @return Tibble of data for training fold i in the k-fold cross validation
#'
get_labelled_data <- function(i, df, cv_idxs){
  fold_idxs <- cv_idxs$train_idxs
  df %>%
    mutate(training_data = ID %in% fold_idxs[[i]]) %>%  #splits data by replicate ID
    mutate(testing_data  = !training_data) %>%
    select(c(time, OD, ID, fungal_ic, technical_reps, biological_reps, blanks, training_data, testing_data, missing_data)) %>%
    unite("fold_ID", fungal_ic:biological_reps, remove = FALSE) %>%
    mutate(fold_ID = as.numeric(factor(fold_ID))) %>%
    arrange(ID)  %>%
    mutate(training_data = training_data & !missing_data) %>%
    mutate(testing_data = testing_data & !missing_data) -> res
  return(res)
}
#' get_metrics
#'
#' Takes a stan fit and data tibble and returns tibble of mean and SEM accuracy over the replicates.
#'
#' @param dat_i Tibble. Tibble that contains the training + testing data set of fold i.
#' @param stan_i Stan object. Stan fit object trained on training data of fold i.
#'
#' @return Tibble of Mean and SEM of metrics (ME, RMSE, MAE, MPE & MAPE) between the replicate posterior predictive distribution (y_tot) and the data provided, for fold i.
#'
get_metrics <- function(dat_i, stan_i) {
  summary_fit <- as.data.frame(summary(stan_i)$summary)
  y_tot <- summary_fit[grep("y_tot", row.names(summary_fit)), ]
  log_lik <- extract_log_lik(stan_i)
  lpds <- colMeans(exp(log_lik))
  dat_i %>%
    mutate(y_rep = y_tot$mean) %>%
    mutate(lpd = log(lpds)) %>%
    filter(!missing_data) %>%  #do not want to eval. performance on missing data
    group_by(fold_ID, training_data, testing_data) %>%
    summarise(accs = c(as.list(accuracy(OD, y_rep)), sum(lpd))) %>%
    mutate(accs = as.numeric(accs)) %>%
    mutate(Metric = c('ME', 'RMSE', 'MAE', 'MPE', 'MAPE', 'LPD')) %>%
    group_by(Metric, training_data, testing_data) %>%
    summarise(Mean = mean(accs), SEM = sd(accs)/sqrt(n())) %>%
    distinct(Metric, training_data, testing_data, Mean, SEM) %>%
    ungroup() %>%
    mutate(Training = training_data) %>%
    select(-c(training_data, testing_data)) -> metrics
  return(metrics)
}
