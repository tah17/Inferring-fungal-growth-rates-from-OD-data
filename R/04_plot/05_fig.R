#
# Script that plots figure 5: modelling calibration overcomes issues
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(tidybayes)
library(rstan)
seed <- 683861
set.seed(seed)
source("R/04_plot/functions.R")

# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD_all_ICs.Rda") -> df
fungal_ic_labels <- c(expression(0~'[N/'~mu~'l]'), expression(1~'[N/'~mu~'l]'), expression(10^1~'[N/'~mu~'l]'), expression(10^2~'[N/'~mu~'l]'), expression(10^3~'[N/'~mu~'l]'))
df %>% mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> df
#
# get posterior samples
#
model <- "logistic_OD_calibration"
source(paste('models/OD/', model, '.data.R', sep = ""))
fit <- readRDS(file=paste("output/OD/", model, "_full_fit.Rda", sep=""))
extracted_fit <- rstan::extract(fit)
y_replicates <- extracted_fit$y_tot
summary_fit <- as.data.frame(summary(fit)$summary)  #extracts a dataframe of the samples from the stan object
f_reps <- summary_fit[grep("^f\\[\\d+\\]", row.names(summary_fit)), ] 
cbind(get_processed_df(df), 
      y_rep_median = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.5)))),
      y_rep_97_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.975)))),
      y_rep_2_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.025)))),
      y_rep_90 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.90)))),
      y_rep_10 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.10)))),
      y_rep_80 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.80)))),
      y_rep_20 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.20))))) -> logistic_post
#
# plot fit to OD
#
logistic_post %>%
  group_by(ID) %>%
  filter(!blanks) %>%
  cbind(f_rep = f_reps$`50%`, f_rep_upper = f_reps$`97.5%`, f_rep_lower = f_reps$`2.5%`) %>%
  ggplot(aes(x = time, y = OD, group = ID)) +
  geom_ribbon(aes(ymin = y_rep_20, ymax = y_rep_80), fill="#525252", alpha = 0.1) +
  geom_ribbon(aes(ymin = y_rep_10, ymax = y_rep_90), fill="#969696", alpha = 0.05) +
  geom_ribbon(aes(ymin = y_rep_2_5, ymax = y_rep_97_5), fill="#bdbdbd", alpha = .01) +
  geom_line(aes(y = y_rep_median), colour="#0091D4", size = 1) +
  geom_point(size=0.1) +
  facet_grid(. ~  fungal_ic_new, labeller="label_parsed") +
  theme_bw(base_size = 11) +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.1, 0.5)) +
  xlab("Time [hrs]") +
  ylab("OD") +
  theme(legend.position="none", text = element_text(size=11)) -> p1
#
# plot underlying f
#
logistic_post %>%
  group_by(ID) %>%
  filter(!blanks) %>%
  cbind(f_rep = f_reps$`50%`, f_rep_upper = f_reps$`97.5%`, f_rep_lower = f_reps$`2.5%`) %>%
  ggplot(aes(x = time, y = f_rep, group = ID)) +
  geom_ribbon(aes(ymin = f_rep_lower, ymax = f_rep_upper), fill="#969696", alpha = 0.05) +
  geom_line(aes(y = f_rep), colour="#0091D4", size = 1) +
  facet_grid(. ~  fungal_ic_new, labeller="label_parsed") +
  theme_bw(base_size = 11) +
  scale_y_log10(breaks=c(1e1, 1e2, 1e3, 1e4, 1e5), labels = 1:5) +
  xlab("Time [hrs]") +
  ylab(expression("log"[10]~paste(f(t), " [N/", mu, "l]"))) +
  theme(legend.position="none", text = element_text(size=11)) -> p2
#
# read in logistic model parameters and plot
#
logistic_models <- c("OD/logistic", "OD/logistic_OD_calibration", "microscopy/logistic", "h1/logistic")
lapply(logistic_models, function(x) readRDS(paste("output/", x, "_params.Rda", sep="")) %>% 
         mutate(Model=x)) %>% 
  rbindlist()  %>%
  add_column(model_type="Logistic") -> parameters_logistic

parameters_logistic %>%
  mutate(data = case_when(Model%in%logistic_models[1:2] ~ "OD", Model==logistic_models[3] ~ "HL", Model==logistic_models[4] ~ "NC")) %>%
  mutate(data_type = case_when(Model%in%logistic_models[1:2] ~ "Indirect", Model%in%logistic_models[3:4] ~ "Direct")) %>%
  mutate(data = factor(data)) %>%
  mutate(data_type = factor(data_type, levels = c("Direct", "Indirect"))) %>%
  mutate(Model = factor(Model, levels=c("OD/logistic_OD_calibration", "OD/logistic", "h1/logistic", "microscopy/logistic"))) %>%
  ggplot(aes(fill = data_type, color = Model, x = .value, y = Model)) +
  stat_pointinterval(point_size=1.7) +
  facet_grid(data_type ~ ., scales="free", space = "free_y", labeller = "label_parsed") +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks=c(1e-2, 1e-1, 1e-0), trans="log10",  labels = c(-2, -1, 0), limits = c(1e-2,1e-0)) +
  scale_y_discrete(breaks=c("microscopy/logistic", "h1/logistic", "OD/logistic", "OD/logistic_OD_calibration"), labels = c("Logistic-HL", "Logistic-NC", "Logistic-OD", "Logistic-OD\n-calibration")) +
  xlab(expression("log"[10]~"Inferred Hyphal Growth Rate")) +
  ylab("Model") +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.right = element_text(angle = 0), axis.text.x = element_text(angle=0), legend.position="none", text=element_text(size=12)) +
  scale_color_manual(name = "Data Type", breaks = logistic_models, values = c("#373A36", "#0091D4", "#31a354", "#31a354")) +
  scale_fill_manual(name = "Data Type", breaks = logistic_models, values = c("#373A36", "#0091D4", "#31a354", "#31a354")) -> p3

p3_space <- plot_grid(NULL, p3, NULL, ncol=1, rel_heights = c(0.1, 1, 0.1))
p3_space2 <- plot_grid(p3_space, NULL, nrow=1, rel_widths = c(0.7, 0.29))

tiff("Fig5.tif", width = 14.5, height = 16, units = "cm", res=300)
plot_grid(p3_space2, p2, p1, labels=c('a', 'b', 'c'), ncol=1, rel_heights = c(1, 0.9, 0.9), label_size=12)
dev.off()

  