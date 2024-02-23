#
# Script that plots figure 3: posterior samples and data
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
seed <- 683861
set.seed(seed)
source("R/04_plot/functions.R")

# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD_all_ICs.Rda") -> df
fungal_ic_labels <- c(expression(0~'[N/'~mu~'l]'), expression(1~'[N/'~mu~'l]'), expression(10^1~'[N/'~mu~'l]'), expression(10^2~'[N/'~mu~'l]'), expression(10^3~'[N/'~mu~'l]'))
df %>% mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> df
#
# get latent var logistic samples
#
model <- "logistic_OD_calibration"
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

logistic_post %>%
  group_by(interaction(ID)) %>%
  ggplot(aes(x = time, y = OD, group = ID)) +
  geom_ribbon(aes(ymin = y_rep_20, ymax = y_rep_80), fill="#525252", alpha = 0.1) +
  geom_ribbon(aes(ymin = y_rep_10, ymax = y_rep_90), fill="#969696", alpha = 0.05) +
  geom_ribbon(aes(ymin = y_rep_2_5, ymax = y_rep_97_5), fill="#bdbdbd", alpha = .01) +
  geom_line(aes(y = y_rep_median), colour="#0091D4", size = 1) +
  geom_point(size=0.1) +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.1, 0.5)) +
  facet_grid(. ~ fungal_ic_new, labeller="label_parsed") +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("OD") +
  theme(legend.position="none", text = element_text(size=11)) -> p

tiff("Fig3.tif", width = 19, height = 5.5, units = "cm", res=300)
p
dev.off()
