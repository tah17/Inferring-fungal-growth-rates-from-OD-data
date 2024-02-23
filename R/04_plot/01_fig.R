#
# Script that plots figure 1: plot of experimental data
#
rm(list = ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)

seed <- 683861
set.seed(seed)
# Read in Data ------------------------------------------------------------
#
# read and process OD data
#
readRDS(file = "data/OD_all_ICs.Rda") -> df_OD
fungal_ic_labels_OD <- c(expression(0~'[N/'~mu~'l]'), expression(1~'[N/'~mu~'l]'), expression(10^1~'[N/'~mu~'l]'), expression(10^2~'[N/'~mu~'l]'), expression(10^3~'[N/'~mu~'l]'))
df_OD %>% mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels_OD)) -> df_OD
#
# read and process hyphal length data
#
# 
readRDS(file = "data/microscopy.Rda") -> df_micro
fungal_ic_labels_micro <- c(expression(1~'[N/'~mu~'l]'))
df_micro %>% mutate(fungal_ic_new = factor(fungal_ic, labels = fungal_ic_labels_micro)) -> df_micro
#
# read in nuclear count data
#
readRDS(file = "data/h1_data.Rda") -> df_h1
# Plot --------------------------------------------------------------------
point_size <- 1.5
point_colour <- "#000000"
line_colour <- "#525252"
#
# plot the OD data by fungal IC
#
df_OD %>%
  ggplot(aes(x = time, y = OD, group = interaction(fungal_ic, technical_reps, biological_reps))) +
  geom_point(shape=20, size=point_size, colour=point_colour) +
  geom_line(alpha=.5, colour=line_colour) +
  facet_grid(. ~ fungal_ic_new, labeller = "label_parsed") +
  scale_y_log10(breaks=c(0.1, 0.3, 0.5), limits=c(0.1, 0.5)) +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("OD") +
  theme(legend.position="none", text = element_text(size=11)) -> p_OD
#
# plot the hyphal length data
#
df_micro %>%
  ggplot(aes(x = time, y = length, group = interaction(ID))) +
  geom_point(shape=20, size=point_size, colour=point_colour) +
  geom_line(alpha=.5, colour=line_colour) +
  facet_grid(. ~ fungal_ic_new, labeller = "label_parsed") +
  scale_y_log10() +
  xlim(0, 25) +
  theme_bw(base_size = 11)+
  xlab("Time [hrs]") +
  ylab(expression(paste("Hyphal length (", mu, "m)"))) +
  theme(legend.position="none", text = element_text(size=11)) -> p_micro
#
# plot nuclear count data
#
df_h1 %>%
  mutate(fungal_ic = factor(fungal_ic, labels=c(expression(10^1~'[N/'~mu~'l]')))) %>%
  ggplot(aes(x = time, y = nuclei, colour= fungal_ic, group = interaction(ID))) +
  geom_point(shape=20, size=point_size, colour=point_colour) +
  geom_line(alpha=.5, colour=line_colour) +
  facet_grid(. ~ fungal_ic, labeller="label_parsed") +
  scale_y_log10() +
  xlim(0, 25) +
  theme_bw(base_size = 11) +
  xlab("Time [hrs]") +
  ylab("Nuclei") +
  theme(legend.position="none", text = element_text(size=11)) -> p_h1

micro_h1_data <- plot_grid(p_micro, p_h1, NULL, labels=c('b', 'c', ''), rel_widths = c(1.2, 1.2, 2), nrow=1, label_size=12)

tiff("Fig1.tif", width = 19, height = 11, units = "cm", res=300)
plot_grid(p_OD, micro_h1_data, ncol=1, labels=c('a', ''), label_size=12)
dev.off()

