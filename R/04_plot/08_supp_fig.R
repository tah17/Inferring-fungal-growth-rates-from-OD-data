#
# Script that plots supp figure S8: prior samples and data
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
# read_prior
#
model <- "logistic_OD_calibration"
source(paste('models/OD/', model, '.data.R', sep = ""))
summary_prior_pred <- readRDS(file=paste("output/OD/", model, "_prior.Rda", sep=""))
y_reps <- summary_prior_pred[grep("^y_tot", row.names(summary_prior_pred)), ] #extracts the y replicates (which are samples from the prior pred) from summary_prior_pred
#
# plots the mean and 95% quantiles of the samples
#
cbind(get_processed_df(df), y_rep = y_reps$mean, y_rep_upper = y_reps$`97.5%`, y_rep_lower = y_reps$`2.5%`) %>%
  mutate(mean_y_upper = mean(y_rep_upper)) %>%   #gets the mean of the quantiles of the y_reps
  mutate(mean_y_lower = mean(y_rep_lower)) %>%
  ggplot(aes(x = time, y = OD, group = interaction(fungal_ic, technical_reps, biological_reps))) +
  geom_ribbon(aes(ymin = y_rep_lower, ymax = y_rep_upper), fill = "lightgrey", alpha = .5) +
  geom_point(size=0.1) +
  facet_grid(. ~ fungal_ic_new, labeller="label_parsed") +
  scale_y_log10() +
  geom_line(alpha=.5) +
  theme_bw(base_size = 12) +
  xlab("Time [hrs]") +
  ylab("OD") -> p

tiff("FigS8.tif", width = 19, height = 5.5, units = "cm", res=300)
p
dev.off()
png("FigS8.png", width = 19, height = 5.5, units = "cm", res=300)
p
dev.off()