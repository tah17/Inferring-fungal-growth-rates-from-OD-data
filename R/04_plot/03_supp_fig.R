#
# Script that plots figure S3: CV per fold fit vis for GP and logistic models
#
rm(list = ls())
library(ggplot2)
library(rstan)
library(data.table)
library(tidyverse)
library(cowplot)
library(forecast)
seed <- 683861
set.seed(seed)
source("R/03_validate/functions.R")

get_fold_fit <- function(fold, model, df, cv_idxs){
  source(paste('models/OD/', model, '.data.R', sep = ""))
  get_labelled_data(fold, df, cv_idxs) %>% add_column(Fold=fold) -> fold_df
  fit <- readRDS(file=paste("OD_output/", IC_used, "/5_fold_", model, "_train", fold, ".rds", sep=""))
  extracted_fit <- rstan::extract(fit)
  y_replicates <- extracted_fit$y_tot
  summary_fit <- as.data.frame(summary(fit)$summary)  #extracts a dataframe of the samples from the stan object
  get_processed_df(fold_df) %>%
    cbind(y_rep_median = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.5)))),
          y_rep_97_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.975)))),
          y_rep_2_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.025)))),
          y_rep_90 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.90)))),
          y_rep_10 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.10)))),
          y_rep_80 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.80)))),
          y_rep_20 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.20))))) -> fit_fold_df
  if (model=="GP") {
    f_reps <- summary_fit[grep("^f\\[\\d+,1\\]", row.names(summary_fit)), ] ## for GP-OD
    f_df <- tibble(f_reps, time=0:25)
    fit_fold_df %>%
      filter(!testing_data) %>%
      full_join(f_df) %>%
      select(-c(`25%`, `75%`, n_eff, Rhat, mean, se_mean, sd)) %>% 
      rename(f_rep = `50%`, f_rep_upper=`97.5%`, f_rep_lower=`2.5%`) -> train_df
    f_preds <- summary_fit[grep("^f_pred\\[\\d+\\]", row.names(summary_fit)), ] ## for GP-OD
    fit_fold_df %>%
      group_by(ID) %>%
      filter(testing_data) %>%
      cbind(f_rep = f_preds$`50%`, f_rep_upper = f_preds$`97.5%`, f_rep_lower = f_preds$`2.5%`) -> test_df
    rbind(train_df, test_df) %>%
      select(-groups) %>%
      filter(fungal_ic!=0) -> tot_df
  } else if (model=="GP_OD_calibration") {
    f_g_reps <- summary_fit[grep("^f_g\\[\\d+\\]", row.names(summary_fit)), ]   ## for GP-OD-calibration
    f_g_df <- tibble(f_g_reps, time=rep(0:25, 2), deriv=rep(c(FALSE, TRUE), each=length(unique(time))))
    fit_fold_df %>%
      full_join(f_g_df) %>%
      filter(fungal_ic!=0, !deriv) %>% 
      select(-c(deriv, `25%`, `75%`, n_eff, Rhat, mean, se_mean, sd)) %>% 
      rename(f_rep = `50%`, f_rep_upper=`97.5%`, f_rep_lower=`2.5%`) -> tot_df
  } else {
    f_reps <- summary_fit[grep("^f\\[\\d+\\]", row.names(summary_fit)), ] 
    fit_fold_df %>%
      group_by(ID) %>%
      filter(!blanks) %>%
      cbind(f_rep = f_reps$`50%`, f_rep_upper = f_reps$`97.5%`, f_rep_lower = f_reps$`2.5%`) -> tot_df
  }
  return(tot_df)
}

# Read in Data ------------------------------------------------------------
fungal_ic_labels <- c(expression(1~'[N/'~mu~'l]'), expression(10^1~'[N/'~mu~'l]'), expression(10^2~'[N/'~mu~'l]'), expression(10^3~'[N/'~mu~'l]'))
IC_used <- "all_ICs"
fold_no <- 5
df <- readRDS(file=paste("data/OD_", IC_used, ".Rda", sep=""))
cv_idxs <- readRDS(file = paste("data/split_idxs_", fold_no, "_fold_", IC_used, ".Rda", sep=""))  #Gets cross validation indices

model_1 <- "GP"
lapply(1:5, function(x) get_fold_fit(x, model_1, df, cv_idxs)) %>% 
  rbindlist() %>%
  add_column(Model = model_1) %>%
  mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> fit_df_1

model_2 <- "GP_OD_calibration"
lapply(1:5, function(x) get_fold_fit(x, model_2, df, cv_idxs)) %>% 
  rbindlist() %>%
  add_column(Model = model_2) %>%
  mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> fit_df_2

rbind(fit_df_1, fit_df_2) %>%
  mutate(Model = factor(Model, levels = c(model_1, model_2), labels = c("GP-OD", "GP-OD\n-calibration"))) -> gp_fit_df

model_3 <- "logistic_OD_calibration"
lapply(1:5, function(x) get_fold_fit(x, model_3, df, cv_idxs)) %>% 
  rbindlist() %>%
  add_column(Model = model_3) %>%
  mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> fit_df_3

model_4 <- "logistic_OD_calibration_no_delay"
lapply(1:5, function(x) get_fold_fit(x, model_4, df, cv_idxs)) %>% 
  rbindlist() %>%
  add_column(Model = model_4) %>%
  mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> fit_df_4

rbind(fit_df_3, fit_df_4) %>%
  mutate(Model = factor(Model, levels = c(model_3, model_4), labels = c("Logistic-OD\n-calibration", "No-delay logistic-\nOD-calibration"))) -> logistic_fit_df

chosen_fold <- 3
#
# plotting underlying f of GPs
#
gp_fit_df %>%
  filter(Fold==chosen_fold) %>%
  ggplot(aes(x = time, y = f_rep, group = interaction(ID, Model))) +
    geom_ribbon(aes(ymin = f_rep_lower, ymax = f_rep_upper), fill="#969696", alpha = 0.05) +
    geom_line(aes(y = f_rep), size = 1) +
    facet_grid(Model ~  fungal_ic_new, labeller=labeller(.cols=label_parsed), scales="free_y") +
    theme_bw(base_size = 11) +
    xlab("Time [hrs]") +
    ylab(expression(paste(g(t), " [N/", mu, "l]"))) +
    theme(legend.position="none", text = element_text(size=11), strip.text.y.right = element_text(angle = 0), strip.background.y = element_blank()) -> p1
#
# plotting fits of GPs
#
gp_fit_df %>%
  filter(Fold==chosen_fold) %>%
  ggplot(aes(x = time, y = OD, group = interaction(ID, Model))) +
  geom_ribbon(aes(ymin = y_rep_20, ymax = y_rep_80), fill="#525252", alpha = 0.1) +
  geom_ribbon(aes(ymin = y_rep_10, ymax = y_rep_90), fill="#969696", alpha = 0.05) +
  geom_ribbon(aes(ymin = y_rep_2_5, ymax = y_rep_97_5), fill="#bdbdbd", alpha = .01) +
  geom_line(aes(y = y_rep_median), size = 1) +
  geom_point(aes(color=training_data), size=0.3) +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.08, 0.5)) +
  facet_grid(Model ~ fungal_ic_new, labeller=labeller(.cols=label_parsed)) +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("OD") +
  scale_color_manual(name = "Data Split", breaks = c(TRUE, FALSE), labels = c("Train", "Test"), values = c("#018571", "#d01c8b")) +
  theme(text = element_text(size=11), strip.text.y.right = element_text(angle = 0), strip.background.y = element_blank()) -> p2

legend <- get_legend(
  # create some space to the left of the legend
  p2 + theme(legend.box.margin = margin(0, 0, 0, 12), legend.key.size = unit(2,"line"), legend.position = "bottom") + guides(colour = guide_legend(override.aes = list(size=3)))
)

gp_plot <- plot_grid(p1, p2 + theme(legend.position="none"), ncol=1)

#
# plotting underlying f of logistics
#
logistic_fit_df %>%
  filter(Fold==chosen_fold) %>%
  ggplot(aes(x = time, y = f_rep, group = interaction(ID, Model))) +
  geom_ribbon(aes(ymin = f_rep_lower, ymax = f_rep_upper), fill="#969696", alpha = 0.05) +
  geom_line(aes(y = f_rep), size = 1) +
  facet_grid(Model ~  fungal_ic_new, labeller=labeller(.cols=label_parsed)) +
  theme_bw(base_size = 11) +
  scale_y_log10(breaks=c(1e1, 1e2, 1e3, 1e4, 1e5), labels = 1:5) +
  xlab("Time [hrs]") +
  ylab(expression("log"[10]~paste(f(t), " [N/", mu, "l]"))) +
  theme(legend.position="none", text = element_text(size=11), strip.text.y.right = element_text(angle = 0), strip.background.y = element_blank()) -> p3
#
# plotting fits
#
logistic_fit_df %>%
  filter(Fold==chosen_fold) %>%
  ggplot(aes(x = time, y = OD, group = interaction(ID, Model))) +
  geom_ribbon(aes(ymin = y_rep_20, ymax = y_rep_80), fill="#525252", alpha = 0.1) +
  geom_ribbon(aes(ymin = y_rep_10, ymax = y_rep_90), fill="#969696", alpha = 0.05) +
  geom_ribbon(aes(ymin = y_rep_2_5, ymax = y_rep_97_5), fill="#bdbdbd", alpha = .01) +
  geom_line(aes(y = y_rep_median), size = 1) +
  geom_point(aes(color=training_data), size=0.3) +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.08, 0.5)) +
  facet_grid(Model ~ fungal_ic_new, labeller=labeller(.cols=label_parsed)) +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("OD") +
  scale_color_manual(name = "Data Split", breaks = c(TRUE, FALSE), labels = c("Train", "Test"), values = c("#018571", "#d01c8b")) +
  theme(legend.position="none", text = element_text(size=11), strip.text.y.right = element_text(angle = 0), strip.background.y = element_blank()) -> p4

logistic_plot <- plot_grid(p3, p4, ncol=1)
tot_plot <- plot_grid(gp_plot, logistic_plot, labels="auto", nrow=1, label_size=12)

tiff("FigS3.tif", width=29.5, height=18, units = "cm", res=300)
plot_grid(tot_plot, legend, rel_heights = c(1, .1),  ncol=1)
dev.off()
png("FigS3.png", width=29.5, height=18, units = "cm", res=300)
plot_grid(tot_plot, legend, rel_heights = c(1, .1),  ncol=1)
dev.off()