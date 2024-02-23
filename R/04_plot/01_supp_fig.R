#
# Script that plots supp figure 1: Exponential-OD and No delay exponential function posterior samples and data
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
seed <- 683861
set.seed(seed)

# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD_all_ICs.Rda") -> df
fungal_ic_labels <- c(expression(0~'[N/'~mu~'l]'), expression(1~'[N/'~mu~'l]'), expression(10^1~'[N/'~mu~'l]'), expression(10^2~'[N/'~mu~'l]'), expression(10^3~'[N/'~mu~'l]'))
df %>% mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> df
#
# get exponential models' samples
#
fit1 <- readRDS(file="output/OD/exponential_full_fit.Rda")
extracted_fit1 <- rstan::extract(fit1)
y_replicates1 <- extracted_fit1$y_tot
cbind(df, y_rep_median = sapply(1:dim(y_replicates1)[2], function(x) unname(quantile(y_replicates1[,x], probs = c(0.5)))),
      y_rep_97_5 = sapply(1:dim(y_replicates1)[2], function(x) unname(quantile(y_replicates1[,x], probs = c(0.975)))),
      y_rep_2_5 = sapply(1:dim(y_replicates1)[2], function(x) unname(quantile(y_replicates1[,x], probs = c(0.025)))),
      y_rep_90 = sapply(1:dim(y_replicates1)[2], function(x) unname(quantile(y_replicates1[,x], probs = c(0.90)))),
      y_rep_10 = sapply(1:dim(y_replicates1)[2], function(x) unname(quantile(y_replicates1[,x], probs = c(0.10)))),
      y_rep_80 = sapply(1:dim(y_replicates1)[2], function(x) unname(quantile(y_replicates1[,x], probs = c(0.80)))),
      y_rep_20 = sapply(1:dim(y_replicates1)[2], function(x) unname(quantile(y_replicates1[,x], probs = c(0.20)))),
      Model="Exponential-OD") -> exp_post

fit2 <- readRDS(file="output/OD/exponential_no_delay_full_fit.Rda")
extracted_fit2 <- rstan::extract(fit2)
y_replicates2 <- extracted_fit2$y_tot
cbind(df, y_rep_median = sapply(1:dim(y_replicates2)[2], function(x) unname(quantile(y_replicates2[,x], probs = c(0.5)))),
      y_rep_97_5 = sapply(1:dim(y_replicates2)[2], function(x) unname(quantile(y_replicates2[,x], probs = c(0.975)))),
      y_rep_2_5 = sapply(1:dim(y_replicates2)[2], function(x) unname(quantile(y_replicates2[,x], probs = c(0.025)))),
      y_rep_90 = sapply(1:dim(y_replicates2)[2], function(x) unname(quantile(y_replicates2[,x], probs = c(0.90)))),
      y_rep_10 = sapply(1:dim(y_replicates2)[2], function(x) unname(quantile(y_replicates2[,x], probs = c(0.10)))),
      y_rep_80 = sapply(1:dim(y_replicates2)[2], function(x) unname(quantile(y_replicates2[,x], probs = c(0.80)))),
      y_rep_20 = sapply(1:dim(y_replicates2)[2], function(x) unname(quantile(y_replicates2[,x], probs = c(0.20)))),
      Model="No-delay-exponential") -> exp_post_no_delay

posteriors <- full_join(exp_post, exp_post_no_delay)

posteriors %>%
  mutate(Model = factor(Model, levels=c("Exponential-OD", "No-delay-exponential"), labels=c("Exponential-OD", "No-delay \nexponential model"))) %>%
  group_by(interaction(ID, Model)) %>%
  ggplot(aes(x = time, y = OD, group = ID)) +
  geom_ribbon(aes(ymin = y_rep_20, ymax = y_rep_80), fill="#525252", alpha = 0.1) +
  geom_ribbon(aes(ymin = y_rep_10, ymax = y_rep_90), fill="#969696", alpha = 0.05) +
  geom_ribbon(aes(ymin = y_rep_2_5, ymax = y_rep_97_5), fill="#bdbdbd", alpha = .01) +
  geom_line(aes(y = y_rep_median), colour="#373A36", size = 1) +
  geom_point(size=0.1) +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.07, 0.5)) +
  facet_grid(Model ~ fungal_ic_new, labeller=labeller(.cols=label_parsed)) +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("OD") +
  theme(legend.position="none", text = element_text(size=11), strip.placement = "outside", strip.background.y = element_blank(), strip.text.y.right = element_text(angle = 0)) -> p

tiff("FigS1.tif", width = 18.5, height = 8, units = "cm", res=300)
p
dev.off()
png("FigS1.png", width = 18.5, height = 8, units = "cm", res=300)
p
dev.off()
