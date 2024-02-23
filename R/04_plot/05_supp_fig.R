#
# Script that plots figure S5: plot of initial growth rates of gompertz
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
# read in gompertz model parameters from fitting to all and high ICs and plot
#
gompertz_models <- c("OD/gompertz", "OD/gompertz_high_IC", "OD/gompertz_OD_calibration", "OD/gompertz_OD_calibration_high_IC", "microscopy/gompertz", "h1/gompertz")
lapply(gompertz_models, function(x) readRDS(paste("output/", x, "_params.Rda", sep="")) %>%
         mutate(Model=x)) %>%
  rbindlist() -> parameters_gompertz

parameters_gompertz %>%
  mutate(data = case_when(Model%in%gompertz_models[1:4] ~ "OD", Model==gompertz_models[5] ~ "Microscopy", Model==gompertz_models[6] ~ "Nuclear Count")) %>%
  mutate(data_type = case_when(Model%in%gompertz_models[1:4] ~ "Indirect", Model%in%gompertz_models[5:6] ~ "Direct")) %>%
  mutate(data = factor(data)) %>%
  mutate(data_type = factor(data_type, levels = c("Direct", "Indirect"))) %>%
  mutate(Model = factor(Model, levels=c( "OD/gompertz_OD_calibration_high_IC", "OD/gompertz_OD_calibration", "OD/gompertz_high_IC", "OD/gompertz", "h1/gompertz", "microscopy/gompertz"))) %>%
  ggplot(aes(fill = data_type, color = data_type, x = .value, y = Model)) +
  stat_pointinterval(point_size=1.7) +
  facet_grid(data_type ~ ., scales="free", space = "free_y", labeller = "label_parsed") +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks=c(1e-2, 1e-1, 1e-0), trans="log10",  labels = c(-2, -1, 0), limits = c(0.008,1e-0)) +
  scale_y_discrete(breaks=c("microscopy/gompertz", "h1/gompertz", "OD/gompertz", "OD/gompertz_high_IC", "OD/gompertz_OD_calibration", "OD/gompertz_OD_calibration_high_IC"), labels = c("Gompertz-HL", "Gompertz-NC", "Gompertz-OD", "Gompertz-OD (High IC)",  "Gompertz-OD\n-calibration", "Gompertz-OD-calibration \n(High IC)")) +
  xlab(expression("log"[10]~"Inferred Initial Hyphal Growth Rate")) +
  ylab("Model") +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.right = element_text(angle = 0), axis.text.x = element_text(angle=0), legend.position="none") +
  scale_color_manual(name = "Data Type", breaks = c("Indirect", "Direct"), labels = c("Indirect", "Direct"), values = c("#373A36", "#31a354")) +
  scale_fill_manual(name = "Data Type", breaks = c("Indirect", "Direct"),  labels = c("Indirect", "Direct"), values = c("#373A36", "#31a354")) -> p

tiff("FigS5.tif", width = 14, height = 6.5, units = "cm", res=300)
p
dev.off()
png("FigS5.png", width = 14, height = 6.5, units = "cm", res=300)
p
dev.off()
