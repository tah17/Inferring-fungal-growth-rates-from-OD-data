#
# Script that plots figure S6: plot of GP growth rates
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(data.table)
library(tidybayes)
seed <- 683861
set.seed(seed)

# Read in Data ------------------------------------------------------------
#
# read in GP model growth rates and plot
#
GP_models <- c("OD/GP_OD_calibration_high_IC", "OD/GP_OD_calibration", "OD/GP", "OD/GP_high_IC", "microscopy/GP", "h1/GP")
#
# read in growth rate from GP fit to high IC OD only
#
readRDS("output/OD/GP_high_IC_params.Rda") %>%  
  mutate(Model="OD/GP_high_IC") %>%
  ungroup() %>%
  select(-rep) %>%
  add_column(model_type="GP") -> parameters_GP_OD_high_IC
#
# read in growth rate from GP fit to OD 
#
readRDS("output/OD/GP_params.Rda") %>%  
  mutate(Model="OD/GP") %>%
  ungroup() %>%
  select(-rep) %>%
  add_column(model_type="GP") %>%
  rbind(parameters_GP_OD_high_IC) -> parameters_GP_OD
#
# read in growth rate from GP-OD-calibration model fit to OD on high IC
#
readRDS("output/OD/GP_OD_calibration_high_IC_params.Rda") %>%  
  mutate(Model="OD/GP_OD_calibration_high_IC") %>%
  ungroup() %>%
  add_column(model_type="GP") -> parameters_GP_OD_calibration_high_IC
#
# read in growth rate from GP-OD-calibration model fit to OD and combine with growth rate from high IC
#
readRDS("output/OD/GP_OD_calibration_params.Rda") %>%  
  mutate(Model="OD/GP_OD_calibration") %>%
  ungroup() %>%
  add_column(model_type="GP") %>%
  rbind(parameters_GP_OD_calibration_high_IC) -> parameters_GP_OD_calibration
#
# read in growth rate from GP fit to hyphal length data
#
readRDS("output/microscopy/GP_params.Rda") %>%
  mutate(Model="microscopy/GP") %>%
  ungroup() %>%
  add_column(model_type="GP")-> parameters_GP_microscopy
#
# read in growth rate from GP fit to nuclear count data
#
readRDS("output/h1/GP_params.Rda") %>%
  mutate(Model="h1/GP") %>%
  add_column(model_type="GP") -> parameters_GP_h1
#
# combine all growth rates and plot
#
rbind(parameters_GP_microscopy, parameters_GP_OD, parameters_GP_h1, parameters_GP_OD_calibration) -> parameters_GP
# Plot ------------------------------------------------------------
parameters_GP %>%
  mutate(Model = factor(Model, levels=c(GP_models[1], GP_models[2], GP_models[4], GP_models[3], GP_models[6], GP_models[5]))) %>%
  mutate(data = case_when(Model%in%GP_models[1:4] ~ "Indirect", Model%in%GP_models[5:6] ~ "Direct")) %>%
  ggplot(aes(fill = data, color = data, group = interaction(.variable, Model), x = .value, y = Model)) +
  stat_pointinterval(point_size=1.7) +
  facet_grid(~data ~ ., scales="free", space = "free", labeller = "label_parsed") +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks=c(1e-2, 1e-1, 1e-0), trans="log10",  labels = c(-2, -1, 0), limits = c(1e-2,1e-0)) +
  scale_y_discrete(breaks=c("microscopy/GP", "h1/GP", "OD/GP", "OD/GP_high_IC", "OD/GP_OD_calibration", "OD/GP_OD_calibration_high_IC"), labels = c("GP-HL", "GP-NC", "GP-OD", "GP-OD (High IC)",  "GP-OD-calibration", "GP-OD-calibration \n(High IC)")) +
  ylab("Model") +
  xlab(expression("log"[10]~"Inferred Growth Rate")) +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.right = element_text(angle = 0), legend.position = "none") +
  scale_color_manual(name = "Data Type", breaks = c("Indirect", "Direct"), labels = c("Indirect", "Direct"), values = c("#373A36", "#31a354")) +
  scale_fill_manual(name = "Data Type", breaks = c("Indirect", "Direct"),  labels = c("Indirect", "Direct"), values = c("#373A36", "#31a354")) -> p_GP

tiff("FigS6.tif", width = 14, height = 8, units = "cm", res=300)
p_GP
dev.off()
png("FigS6.png", width = 14, height = 8, units = "cm", res=300)
p_GP
dev.off()