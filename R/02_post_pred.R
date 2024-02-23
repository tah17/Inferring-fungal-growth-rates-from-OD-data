#
# Script that runs a full posterior inference and posterior check (to be visually assessed) for a model in models/
#
rm(list = ls())
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(bayesplot)
library(loo)
library(tidybayes)

seed <- 683861
set.seed(seed)
save_fit <- TRUE

# Read in Data ------------------------------------------------------------
data_type <- "OD"  ## pick type of data out of "OD", "microscopy" (hyphal length) and "h1" (nuclear count)
if (identical(data_type, "OD")) {
  datafile <- "data/OD_all_ICs.Rda"
  readRDS(file = datafile) %>%
    mutate(missing_data = is.na(OD)) -> df_no_labels
} else if (identical(data_type, "microscopy")) {
  datafile <- "data/microscopy.Rda"
  readRDS(file = datafile) %>%
    mutate(missing_data = is.na(length)) -> df_no_labels
} else if (identical(data_type, "h1")) {
  datafile <- "data/h1_data.Rda"
  readRDS(file = datafile) %>%
    mutate(missing_data = is.na(nuclei)) -> df_no_labels
}
## adds training and testing labels
df_no_labels %>%
  mutate(training_data = TRUE) %>%
  mutate(testing_data = FALSE) %>%
  filter(training_data | testing_data) %>%
  mutate(training_data = training_data & !missing_data) %>%
  mutate(testing_data = testing_data & !missing_data) %>%
  # filter(fungal_ic%in%c(0, 2e5)) %>%  # when fitting to only high IC
  arrange(ID) -> df

# Specify Model -----------------------------------------------------------
model <- "logistic_OD_calibration"
stanfile <- paste('models/', data_type, '/', model, '.stan', sep = "")  #read in stan code
#
# read in associated model.data.R file which specified the stan data, iterations and number of chains.
# also defines the parameters in the model
#
source(paste('models/', data_type, '/', model, '.data.R', sep = ""))

# Fit Model ---------------------------------------------------------------
#
# i.e. get samples from the posterior
#
fit <- stan(file = stanfile,
            data = get_stan_data(df),
            seed = seed,
            chains = no_of_chains,
            iter = iter)
#
# stan diagnostics and plots
#
traceplot(fit)
plot(fit, pars=params)
#
# save fit and parameter draws for later plots
#
if (save_fit) {
  saveRDS(fit, file=paste("output/", data_type, '/', model, "_full_fit.Rda", sep=""))  #stores model
  saveRDS(gather_draws(fit, cbind(beta)), file=paste("output/", data_type, '/', model, "_params.Rda", sep="")) #stores growth parameter
  # saveRDS(gather_draws(fit, cbind(growth_rate)), file=paste("output/", data_type, '/', model, "_params.Rda", sep="")) #stores growth parameter for GP
}
# Posterior Predictive Check ----------------------------------------------
#
# visualize samples from the posterior pred. distribution
#
# real data vs the samples from the posterior
#
summary_fit <- as.data.frame(summary(fit)$summary)  #extracts a dataframe of the samples from the stan object
y_reps <- summary_fit[grep("^y_tot", row.names(summary_fit)), ]  #extracts the y_reps (which are samples from the posterior pred) from summary_fit
get_processed_df(df) %>%
  cbind(y_rep = y_reps$`50%`, y_rep_upper = y_reps$`97.5%`, y_rep_lower = y_reps$`2.5%`) %>%
  mutate(mean_y_upper = mean(y_rep_upper)) %>%  #gets the mean of the quantiles of the y_reps
  mutate(mean_y_lower = mean(y_rep_lower)) %>%
  ggplot(aes(x = time, y = {if (identical(data_type, "OD")) OD else if (identical(data_type, "microscopy")) length else if (identical(data_type, "h1")) nuclei}, group = ID)) +
  geom_point(size=0.5) +
  geom_ribbon(aes(ymin = y_rep_lower, ymax = y_rep_upper), fill = "lightgray", alpha = .1) +
  geom_line(aes(y = y_rep), color="red", size = 0.5) +
  facet_grid(. ~ fungal_ic) +
  scale_y_log10() +
  theme_bw(base_size = 20)
