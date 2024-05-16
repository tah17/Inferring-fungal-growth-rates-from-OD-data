#
# Script that plots supp figure S10: prior samples and data
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(tidybayes)
seed <- 683861
set.seed(seed)

# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD_all_ICs.Rda") -> df
fungal_ic_labels <- c(expression(0~'[N/'~mu~'l]'), expression(1~'[N/'~mu~'l]'), expression(10^1~'[N/'~mu~'l]'), expression(10^2~'[N/'~mu~'l]'), expression(10^3~'[N/'~mu~'l]'))
df %>% mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels)) -> df
#
# read_prior
#
model <- "logistic_OD_calibration"
source(paste('models/OD/', model, '.data.R', sep = ""))
prior_pred <- readRDS(file=paste("output/OD/", model, "_full_prior.Rda", sep=""))
gather_draws(prior_pred, !!str2lang(pop_params2lang)) %>%
  add_column(distribution = "prior") -> priors

readRDS(file=paste('output/OD/', model, "_prior_draw.Rda", sep="")) %>%
  add_column(distribution = "fake_data_draw") %>%
  mutate(Distribution = factor(distribution, levels = c("fake_data_draw"), labels = c("Sample from prior \nthat generates fake data"))) -> fake_data_check_draw


readRDS(file=paste('output/OD/', model, "_prior_draw_fit.Rda", sep="")) %>%
  add_column(distribution = "posterior") -> fake_data_check_post

rbind(priors, fake_data_check_post) %>%
  mutate(Distribution = factor(distribution, levels = c("prior", "posterior"), labels = c("Prior", "Posterior \n(fake data check)"))) -> draws
#
# plots
#
ggplot() +
  stat_pointinterval(data = draws, aes(x =.variable, y = .value, color=Distribution), position="dodge", .width = c(0.80, 0.95)) +
  geom_point(data = fake_data_check_draw, aes(y = .value, x = .variable, shape=Distribution), size=2, stroke=1) +
  theme_bw(base_size = 11) +
  scale_color_brewer(palette="Dark2") +
  scale_x_discrete(breaks = c("delta_tilde", "beta", "sigma_meas", "L", "basal", "tau"),
                  labels = c(expression(tilde(delta)), expression(beta), expression(sigma), "L", "B",  expression(tau))) +
  scale_shape_manual(name = "", values=4) +
  xlab("Parameter") +
  ylab("Value") -> p

tiff("FigS10.tif", width = 19, height = 7, units = "cm", res=300)
p
dev.off()
png("FigS10.png", width = 19, height = 7, units = "cm", res=300)
p
dev.off()
#
# get specific CIs for logistic-OD-calibration model
#
draws %>%
  filter(distribution=="posterior") %>%
  group_by(.variable) %>%
  summarise(median_qi(.value, .width = c(.95))) -> post_cis
