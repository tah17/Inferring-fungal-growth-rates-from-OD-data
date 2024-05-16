#
# Script that plots figure S4: plot of growth rates on high IC
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
# read in logistic model parameters from fitting to high IC and plot
#
logistic_models <- c("OD/logistic_high_IC", "OD/logistic_OD_calibration_high_IC", "microscopy/logistic", "h1/logistic")
lapply(logistic_models, function(x) readRDS(paste("output/", x, "_params.Rda", sep="")) %>% 
         mutate(Model=x)) %>% 
  rbindlist()  %>%
  add_column(model_type="Logistic") -> parameters_logistic

parameters_logistic %>%
  mutate(data = case_when(Model%in%logistic_models[1:2] ~ "OD", Model==logistic_models[3] ~ "Microscopy", Model==logistic_models[4] ~ "Nuclear Count")) %>%
  mutate(data_type = case_when(Model%in%logistic_models[1:2] ~ "Indirect", Model%in%logistic_models[3:4] ~ "Direct")) %>%
  mutate(data = factor(data)) %>%
  mutate(data_type = factor(data_type, levels = c("Direct", "Indirect"))) %>%
  mutate(Model = factor(Model, levels=c("OD/logistic_high_IC", "OD/logistic_OD_calibration_high_IC", "h1/logistic", "microscopy/logistic"))) %>%
  ggplot(aes(fill = data_type, color = Model, x = .value, y = Model)) +
  stat_pointinterval(point_size=1.7, .width = c(0.80, 0.95)) +
  facet_grid(data_type ~ ., scales="free", space = "free_y", labeller = "label_parsed") +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks=c(1e-2, 1e-1, 1e-0), trans="log10",  labels = c(-2, -1, 0), limits = c(1e-2,1e-0)) +
  scale_y_discrete(breaks=c("microscopy/logistic", "h1/logistic", "OD/logistic_high_IC", "OD/logistic_OD_calibration_high_IC"), labels = c("Logistic-HL", "Logistic-NC", "Logistic-OD", "Logistic-OD\n-calibration")) +
  xlab(expression("log"[10]~"Inferred Hyphal Growth Rate")) +
  ylab("Model") +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.right = element_text(angle = 0), axis.text.x = element_text(angle=0), legend.position="none") +
  scale_color_manual(name = "Data Type", breaks = logistic_models, values = c("#373A36", "#0091D4", "#31a354", "#31a354")) +
  scale_fill_manual(name = "Data Type", breaks = logistic_models, values = c("#373A36", "#0091D4", "#31a354", "#31a354")) -> p

tiff("FigS4.tif", width = 12, height = 5, units = "cm", res=300)
p
dev.off()
png("FigS4.png", width = 12, height = 5, units = "cm", res=300)
p
dev.off()

#
# get specific CIs for logistic-OD-calibration model
# 
# parameters_logistic %>%
#   filter(Model==logistic_models[2]) %>%
#   pull(.value) %>%
#   median_qi(.width = c(.80, .95))
