#
# Script that runs a prior predictive check (to be visually assessed) and a
# fake data check for a model in models/
#
rm(list = ls())
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(bayesplot)
library(data.table)
library(tidybayes)
library(posterior)
library(RColorBrewer)

seed <- 683861
set.seed(seed)

save_prior <- FALSE  ## set to TRUE to save the prior and fake data check results of a model
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
  # filter(fungal_ic%in%c(0, 2e5)) %>%  # uncomment when fitting to only high IC
  arrange(ID) -> df

# Specify Model -----------------------------------------------------------
model <- "logistic_OD_calibration"
stanfile <- paste('models/', data_type, '/', model, '.stan', sep = "")  #read in stan code
#
# read in associated model.data.R file which specified the stan data, iterations and number of chains.
# also defines the parameters in the model
#
source(paste('models/', data_type, '/', model, '.data.R', sep = ""))
# Sample from the Prior Predictive Distribution -------------------------------------------------
prior_pred <- stan(file = stanfile,
                   data = get_stan_data(df, run_likelihood = FALSE),
                   seed = seed,
                   chains = no_of_chains,
                   iter = iter)
#
# check hmc diagnostics
#
plot(prior_pred, pars=params)
pairs(prior_pred, pars = c(params, "lp__"))
# Plot Samples  ---------------------------------
#
# visualize samples from the prior pred. distribution
#
# Real data vs the samples from the prior
#
summary_prior_pred <- as.data.frame(summary(prior_pred)$summary) #extracts a dataframe of the samples from the stan object
y_reps <- summary_prior_pred[grep("^y_tot", row.names(summary_prior_pred)), ] #extracts the y replicates (which are samples from the prior pred) from summary_prior_pred
#
# plots the mean and 95% quantiles of the samples
#
cbind(df, y_rep = y_reps$mean, y_rep_upper = y_reps$`97.5%`, y_rep_lower = y_reps$`2.5%`) %>%
  group_by(time) %>%
  mutate(mean_y_upper = mean(y_rep_upper)) %>%   #gets the mean of the quantiles of the y_reps
  mutate(mean_y_lower = mean(y_rep_lower)) %>%
  ggplot(aes(x = time, y = {if (identical(data_type, "OD")) OD else if (identical(data_type, "microscopy")) length else if (identical(data_type, "h1")) nuclei}, group = interaction(fungal_ic, technical_reps, biological_reps))) +
  geom_ribbon(aes(ymin = y_rep_lower, ymax = y_rep_upper), fill = "lightgrey", alpha = .5) +
  geom_point(size=0.5) +
  facet_grid(. ~ fungal_ic) +
  scale_y_log10() +
  geom_line(alpha=.5) +
  theme_bw(base_size = 25) +
  xlab("Time [hrs]") +
  ylab("Values")
#
# save prior
#
if (save_prior) {
  saveRDS(prior_pred, file=paste("output/", data_type, '/', model, "_full_prior.Rda", sep=""))  #stores model prior
  saveRDS(summary_prior_pred, file=paste("output/", data_type, '/', model, "_prior.Rda", sep=""))  #stores model prior
}
# Sample Fake Data ------------------------------------------
#
# draws one parameter sample from the prior predictive distribution
#
draw <- sample(1:(iter/2), 1)
gather_draws(prior_pred, y_tot[1]) %>%
  filter(.draw == draw) -> y_rep_draws
y_rep_draw <- y_rep_draws$.value[[1]]
#
# plots the sample
#
get_processed_df(df) %>%
  cbind(y_rep_draw = y_rep_draw) %>%
  ggplot(aes(x = time, y = {if (identical(data_type, "OD")) OD else if (identical(data_type, "microscopy")) length else if (identical(data_type, "h1")) nuclei}, group = ID)) +
  geom_point(aes(y = y_rep_draw), size = 0.5) +
  facet_grid(. ~ fungal_ic) +
  scale_y_log10() +
  theme_bw(base_size = 25)
#
# generates fake data with this sample
#
if (identical(data_type, "OD")) {
  get_processed_df(df) %>%
    mutate(OD = y_rep_draw) -> fake_df
} else if (identical(data_type, "microscopy")) {
  get_processed_df(df) %>%
    mutate(length = y_rep_draw) -> fake_df
} else if (identical(data_type, "h1")) {
  get_processed_df(df) %>%
    mutate(nuclei = y_rep_draw) -> fake_df
}
#
# gets the population level parameters that generated this sample, for example if the basal OD level for each replicate follows
# a Normal distribution around the mean basal OD with a scale sigma_mu - pop_params_draws contains the value of sigma_mu
#
# pop_params2lang should be specified as a string of a list of variables in model.data.R
# e.g. pop_params2lang <- "cbind(sigma_mu)"
#
gather_draws(prior_pred, !!str2lang(pop_params2lang)) %>%
  filter(.draw == draw) -> pop_params_draws
if (save_prior) {
  saveRDS(pop_params_draws, file=paste("output/", data_type, '/', model, "_prior_draw.Rda", sep=""))  #stores "true" model pop params
}
#
# gets the individual level parameters that generated this sample, for example, if a replicate, j, had a basal OD level B_j
# ind_params_draws contains B_j
#
# ind_params_draws should be specified as a string of a list of variables in model.data.R
# ind_params2lang <- "cbind(B)[rep]"
#
gather_draws(prior_pred, !!str2lang(ind_params2lang)) %>%  #!! - evaluate as literal
  filter(.draw == draw) %>%
  mutate(".truevalue" = .value) %>%
  select(c(.variable, rep, .truevalue)) -> ind_params_draws
# Fit to Fake Data --------------------------------------------------------
#
# fit the model to its own fake data
#
fake_data_check <- stan(file = stanfile,
                        data = get_stan_data(fake_df, run_likelihood = TRUE),
                        seed = seed,
                        chains = no_of_chains,
                        iter = iter)
#
# check hmc diagnostics e.g. look for divergent transitions etc
#
summary(fake_data_check)
traceplot(fake_data_check, pars = params)
pairs(fake_data_check, pars = params)
mcmc_areas(fake_data_check, pars = params, prob = 0.8)

# Evaluate Fake Data Fit -----------------------------------------------------
#
# gets samples from the posteriors of the fake data check for both the
# population and individual level parameters
#
gather_draws(fake_data_check, !!str2lang(pop_params2lang)) -> pop_fake_data_check
if (save_prior) {
  saveRDS(pop_fake_data_check, file=paste("output/", data_type, '/', model, "_prior_draw_fit.Rda", sep=""))  #stores estimated model params from fake data
}
gather_draws(fake_data_check, !!str2lang(ind_params2lang)) -> ind_fake_data_check
#
# plots the population parameters' credible intervals against the true values (black dots)
#
ggplot() +
  stat_pointinterval(data = pop_fake_data_check, aes(x =.variable, y = .value, colour = .variable)) +
  geom_point(data = pop_params_draws, aes(y = .value, x = .variable), shape=4, size=2, stroke=1) +
  theme_bw(base_size = 20) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="none") +
  xlab("Parameter") +
  ylab("Value")
#
# plots a individual parameters' credible intervals against the true values (black dots)
#
ind_fake_data_check %>%
  full_join(ind_params_draws) %>%
  ggplot(aes(y = factor(rep), x = .value, colour = .variable)) +
    stat_pointinterval() +
    facet_grid(. ~ .variable) +
    geom_point(aes(x = .truevalue, y = factor(rep)), colour = "black") +
    theme_bw(base_size = 15)
#
# extracts the posterior of the fake data check and plots the posterior predictive samples (red)
# and their 95% quantiles (shaded grey) against the true fake data - fake_df (black)
#
summary_f_data_check <- as.data.frame(summary(fake_data_check)$summary)
y_reps <- summary_f_data_check[grep("^y_tot", row.names(summary_f_data_check)), ]
fake_df %>%
  cbind(y_rep = y_reps$mean, y_rep_upper = y_reps$`97.5%`, y_rep_lower = y_reps$`2.5%`) %>%
  ggplot(aes(x = time, y = {if (identical(data_type, "OD")) OD else if (identical(data_type, "microscopy")) length else if (identical(data_type, "h1")) nuclei}, group=ID)) +
  geom_ribbon(aes(ymin = y_rep_lower, ymax = y_rep_upper), fill = "lightgrey", alpha = .3) +
  geom_line(aes(x = time, y = y_rep), colour="red") +
  geom_point(size=0.5) +
  facet_grid(. ~ fungal_ic, labeller="label_parsed") +
  scale_y_log10() +
  theme_bw(base_size = 25)
