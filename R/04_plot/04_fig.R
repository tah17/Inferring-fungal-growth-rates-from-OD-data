#
# Script that plots figure 4: model validation, and supp fig 2: model validation on high_IC
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(loo)
library(rlang)
library(data.table)
library(forecast)
library(cowplot)
source("R/04_plot/functions.R")
seed <- 683861
set.seed(seed)

fold_no <- 5
IC_used <- "all_ICs"
readRDS(file=paste("data/OD_", IC_used, ".Rda", sep="")) %>% drop_na() -> df

model_list <- c("RW",
                "mixed_logistic_OD_calibration",
                "gompertz",
                "exponential", 
                "gompertz_OD_calibration",
                "GP",
                "logistic",
                "exponential_OD_calibration",
                "logistic_OD_calibration_no_delay",
                "logistic_OD_calibration",
                "GP_OD_calibration")

lapply(model_list, function(x) readRDS(paste('output/', IC_used, "/", fold_no, "_fold_", x, '_fit_stats.rds', sep = "")) %>% mutate(Model = x)) %>% rbindlist() %>%
  mutate(reference = case_when(Model=="logistic_OD_calibration" ~ "Our Model", 
                               Model %in% c("mixed_logistic_OD_calibration", "logistic_OD_calibration_no_delay") ~ "Extensions",
                               Model %in% c("gompertz_OD_calibration", "exponential_OD_calibration", "GP_OD_calibration") ~ "Reference Models \nwith Calibration",
                               Model == "logistic" ~ "Reference Logistic Model \nwithout Calibration",
                               Model %in% c("RW", "gompertz", "exponential", "GP") ~ "Reference Models \nwithout Calibration")) %>%
  mutate(reference = factor(reference, levels=c("Our Model", "Reference Logistic Model \nwithout Calibration", "Reference Models \nwithout Calibration", "Reference Models \nwith Calibration", "Extensions"))) %>%
  mutate(Model = factor(Model, levels=c("GP", "gompertz", "exponential", "RW", "logistic_OD_calibration", "logistic", "GP_OD_calibration", "gompertz_OD_calibration", "exponential_OD_calibration", "logistic_OD_calibration_no_delay", "mixed_logistic_OD_calibration"))) -> res_list

if (IC_used =="high_IC") {
  plot_RMSE_LPD(res_list %>% filter(!Model%in%c("RW", "gompertz_OD_calibration", "exponential_OD_calibration", "GP_OD_calibration", "mixed_logistic_OD_calibration", "logistic_OD_calibration_no_delay")), "#084594", "#1d91c0", df) -> p
} else {
  plot_RMSE_LPD(res_list %>% filter(Model!="RW"), "#084594", "#1d91c0", df) -> p
}
legend <- get_legend(
  # create some space to the left of the legend
  p + theme(legend.box.margin = margin(0, 0, 0, 12), legend.key.size = unit(2,"line"), legend.position = "bottom")
)
#
# save fig
#
if (IC_used =="high_IC") {
  tiff("FigS2.tif", width = 16, height = 8.8, units = "cm", res=300)
  plot_grid(p + theme(legend.position="none") + scale_x_log10(breaks=c(0.01, 0.03, 0.1, 1, 10, 30, 100), labels = c(-2, -1.52, -1, 0, 1, 1.48, 2)), rel_heights = c(1, .1), legend, ncol=1)
  dev.off()
  png("FigS2.png", width = 16, height = 8.8, units = "cm", res=300)
  plot_grid(p + theme(legend.position="none") + scale_x_log10(breaks=c(0.01, 0.03, 0.1, 1, 10, 30, 100), labels = c(-2, -1.52, -1, 0, 1, 1.48,  2)), rel_heights = c(1, .1), legend, ncol=1)
  dev.off()
} else { 
  tiff("Fig4.tif", width = 16, height = 13, units = "cm", res=300)
  plot_grid(p + theme(legend.position="none") + scale_x_log10(breaks=c(0.03, 0.1, 0.3, 1, 3, 10, 30, 100), labels = c(-1.52, -1, -0.52, 0, 0.48, 1, 1.48,  2)), rel_heights = c(1, .1), legend, ncol=1)
  dev.off()
}

