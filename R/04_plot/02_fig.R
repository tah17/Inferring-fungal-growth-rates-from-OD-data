#
# Script that plots figure 2: showing issues with directly fitting logistic model to OD
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(tidybayes)
seed <- 683861
set.seed(seed)
source("R/04_plot/functions.R")

# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD_all_ICs.Rda") -> df
fungal_ic_labels <- c(expression(0~'[N/'~mu~'l]'), expression(1~'[N/'~mu~'l]'), expression(10^1~'[N/'~mu~'l]'), expression(10^2~'[N/'~mu~'l]'), expression(10^3~'[N/'~mu~'l]'))
df %>% mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> df
#
# get logistic samples
#
model <- "logistic"
source(paste('models/OD/', model, '.data.R', sep = ""))
fit <- readRDS(file=paste("output/OD/", model, "_full_fit.Rda", sep=""))
extracted_fit <- rstan::extract(fit)
y_replicates <- extracted_fit$y_tot
cbind(get_processed_df(df), y_rep_median = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.5)))),
      y_rep_97_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.975)))),
      y_rep_2_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.025)))),
      y_rep_90 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.90)))),
      y_rep_10 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.10)))),
      y_rep_80 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.80)))),
      y_rep_20 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.20))))) -> logistic_post
#
# plot fit to all ICs
#
logistic_post %>%
  group_by(interaction(ID)) %>%
  ggplot(aes(x = time, y = OD, group = ID)) +
  geom_ribbon(aes(ymin = y_rep_20, ymax = y_rep_80), fill="#525252", alpha = 0.1) +
  geom_ribbon(aes(ymin = y_rep_10, ymax = y_rep_90), fill="#969696", alpha = 0.05) +
  geom_ribbon(aes(ymin = y_rep_2_5, ymax = y_rep_97_5), fill="#bdbdbd", alpha = .01) +
  geom_line(aes(y = y_rep_median), colour="#373A36", size = 1) +
  geom_point(size=0.1) +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.08, 0.5)) +
  facet_grid(. ~ fungal_ic_new, labeller="label_parsed") +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("OD") +
  theme(legend.position="none", text = element_text(size=11)) -> p1
#
# get high IC posterior samples
#
model <- "logistic_high_IC"
fit_high_IC <- readRDS(file=paste("output/OD/", model, "_full_fit.Rda", sep=""))
extracted_fit <- rstan::extract(fit_high_IC)
y_replicates <- extracted_fit$y_tot
cbind(get_processed_df(df) %>% filter(fungal_ic%in%c(0, 2e5)), y_rep_median = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.5)))),
      y_rep_97_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.975)))),
      y_rep_2_5 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.025)))),
      y_rep_90 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.90)))),
      y_rep_10 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.10)))),
      y_rep_80 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.80)))),
      y_rep_20 = sapply(1:dim(y_replicates)[2], function(x) unname(quantile(y_replicates[,x], probs = c(0.20))))) -> logistic_post_high_IC
#
# plot fit to only IC=2e5
#
logistic_post_high_IC %>%
  group_by(ID) %>%
  filter(fungal_ic==2e5) %>%
  ggplot(aes(x = time, y = OD, group = ID)) +
  geom_ribbon(aes(ymin = y_rep_20, ymax = y_rep_80), fill="#525252", alpha = 0.1) +
  geom_ribbon(aes(ymin = y_rep_10, ymax = y_rep_90), fill="#969696", alpha = 0.05) +
  geom_ribbon(aes(ymin = y_rep_2_5, ymax = y_rep_97_5), fill="#bdbdbd", alpha = .01) +
  geom_line(aes(y = y_rep_median), colour="#373A36", size = 1) +
  geom_point(size=0.1) +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.1, 0.5)) +
  facet_grid(. ~ fungal_ic_new, labeller="label_parsed") +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("OD") +
  theme(legend.position="none", text = element_text(size=11)) -> p2
#
# read in logistic model parameters and plot
#
logistic_models <- c("OD/logistic", "OD/logistic_high_IC", "microscopy/logistic", "h1/logistic")
lapply(logistic_models, function(x) readRDS(paste("output/", x, "_params.Rda", sep="")) %>% 
         mutate(Model=x)) %>% 
  rbindlist()  %>%
  add_column(model_type="Logistic") -> parameters_logistic

parameters_logistic %>%
  filter(.variable=="beta") %>%
  mutate(data = case_when(Model==logistic_models[1:2] ~ "OD", Model==logistic_models[3] ~ "HL", Model==logistic_models[4] ~ "NC")) %>%
  mutate(data_type = case_when(Model%in%logistic_models[1:2] ~ "Indirect", Model%in%logistic_models[3:4] ~ "Direct")) %>%
  mutate(data = factor(data)) %>%
  mutate(data_type = factor(data_type, levels = c("Direct", "Indirect"))) %>%
  mutate(Model = factor(Model, levels=c("OD/logistic_high_IC", "OD/logistic", "h1/logistic", "microscopy/logistic"))) %>%
  ggplot(aes(fill = data_type, color = Model, x = .value, y = Model)) +
  stat_pointinterval(point_size=1.5) +
  facet_grid(data_type ~ ., scales="free", space = "free_y", labeller = "label_parsed") +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks=c(1e-2, 1e-1, 1e-0), trans="log10",  labels = c(-2, -1, 0), limits = c(1e-2,1e-0)) +
  scale_y_discrete(breaks=c("microscopy/logistic", "h1/logistic", "OD/logistic", "OD/logistic_high_IC"), labels = c("Logistic-HL", "Logistic-NC", "Logistic-OD",  "Logistic-OD (High IC)")) +
  xlab(expression("log"[10]~"Inferred Hyphal Growth Rate")) +
  ylab("Model") +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.right = element_text(angle = 0), axis.text.x = element_text(angle=0), legend.position="none", text=element_text(size=12)) +
  scale_color_manual(name = "Data Type", breaks = logistic_models, values = c("#373A36", "#373A36", "#31a354", "#31a354")) +
  scale_fill_manual(name = "Data Type", breaks = logistic_models, values = c("#373A36", "#373A36", "#31a354", "#31a354")) -> p3

p3_space <- plot_grid(NULL, p3, NULL, ncol=1, rel_heights = c(0.1, 1, 0.1))
p2_3 <- plot_grid(p2, p3_space, NULL, labels=c('b', 'c', ''), nrow=1, rel_widths = c(0.25, 0.6, 0.1), label_size=12)

tiff("Fig2.tif", width = 19, height = 11, units = "cm", res=300)
plot_grid(p1, p2_3, ncol=1, labels=c('a', ''), label_size=12)
dev.off()

