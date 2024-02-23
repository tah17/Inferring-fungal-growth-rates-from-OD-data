#
# Script that plots figure S7: plot of prior growth rates sensitivity
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(data.table)
library(tidybayes)
library(grid)
seed <- 683861
set.seed(seed)

# Read in Data ------------------------------------------------------------
#
# read in logistic model parameters and plot
#
logistic_models <- c("OD/logistic", "OD/logistic_OD_calibration", "microscopy/logistic", "h1/logistic")
lapply(logistic_models, function(x) readRDS(paste("output/", x, "_params.Rda", sep="")) %>% 
         mutate(Model=x)) %>% 
  rbindlist()  %>%
  add_column(prior_type="Normal") -> parameters_normal

lapply(logistic_models, function(x) readRDS(paste("output/", x, "_cauchy_beta_params.Rda", sep="")) %>% 
         mutate(Model=x)) %>% 
  rbindlist()  %>%
  add_column(prior_type="Cauchy") -> parameters_cauchy

parameters_logistic <- rbind(parameters_cauchy, parameters_normal)

parameters_logistic %>%
  mutate(data = case_when(Model%in%logistic_models[1:2] ~ "OD", Model==logistic_models[3] ~ "HL", Model==logistic_models[4] ~ "NC")) %>%
  mutate(data_type = case_when(Model%in%logistic_models[1:2] ~ "Indirect", Model%in%logistic_models[3:4] ~ "Direct")) %>%
  mutate(data = factor(data)) %>%
  mutate(data_type = factor(data_type, levels = c("Direct", "Indirect"))) %>%
  mutate(Model = factor(Model, levels=c("OD/logistic", "OD/logistic_OD_calibration", "h1/logistic", "microscopy/logistic"))) %>%
  ggplot(aes(fill = data_type, color = prior_type, group=prior_type, x = .value, y = Model)) +
  stat_pointinterval(position = position_dodge(0.3), point_size=1.7) +
  facet_grid(data_type ~ ., scales="free", space = "free_y", labeller = "label_parsed") +
  theme_bw(base_size = 11) +
  scale_x_continuous(breaks=c(1e-2, 1e-1, 1e-0), trans="log10",  labels = c(-2, -1, 0), limits = c(1e-2,1e-0)) +
  scale_y_discrete(breaks=c("microscopy/logistic", "h1/logistic", "OD/logistic", "OD/logistic_OD_calibration"), labels = c("Logistic-HL", "Logistic-NC", "Logistic-OD", "Logistic-OD\n-calibration")) +
  xlab(expression("log"[10]~"Inferred Hyphal Growth Rate")) +
  ylab("Model") +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.right = element_text(angle = 0), axis.text.x = element_text(angle=0)) +
  scale_color_manual(name = "Prior", breaks = c("Normal", "Cauchy"), values = c("#7b3294", "#008837")) +
  scale_fill_manual(name = "Prior", breaks = c("Normal", "Cauchy"), values = c("#7b3294", "#008837")) -> p

tiff("FigS7.tif", width = 14, height = 6.5, units = "cm", res=300)
p
dev.off()
png("FigS7.png", width = 14, height = 6.5, units = "cm", res=300)
p
dev.off()