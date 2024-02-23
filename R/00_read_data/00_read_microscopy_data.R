#
# Script that reads in the hyphal length data of the microscopy experiment, performed by Natasha Motsi.
#
rm(list = ls())
library(tidyverse)
library(readxl)
library(ggplot2)
#'
#' read_data
#'
#' Reads the hyphal length data 
#'
#' @param excel_file String. The name of the excel spreadsheet the data holds.
#' @param sheet_name String. Sheet name of excel sheet that contains raw data.
#' @param data_range String. Range of cells in the excel sheet that contain the data.
#' @param fungal_ic. List of doubles. The initial concentration of spores.
#' @param rep_labels. List of strings. Technical replicates of data in the excel sheet. Must be in the same order as the data.
#'
#' @return A data frame of hyphal length readings.
#'
read_data <- function(excel_file, sheet_name, data_range, fungal_ic, rep_labels){
  read_excel(excel_file, sheet = sheet_name, range = data_range) %>% 
    select_if(~!all(is.na(.))) -> df_raw
  df_raw %>%
    rename_with(.cols = 2:dim(df_raw)[2], ~ rep_labels) %>%
    rename(time=Time) %>%
    pivot_longer(-time, names_to="replicate", values_to="length") %>%
    add_column(fungal_ic = fungal_ic) -> df  
  return(df)
}
#
# read in raw data
#
excel_file <- "data/Hyphal\ Length.xlsx"
sheet_name <- "2X10^4"
data_range <- "A1:D20"
fungal_ic <- 1 
rep_labels <- c(paste("R", 1:3, sep=""))
df_raw <- read_data(excel_file, sheet_name, data_range, fungal_ic, rep_labels)
#
# join and save data
#
df_raw %>%
  unite("ID", c(fungal_ic, replicate), remove = FALSE) %>%
  mutate(ID = as.numeric(factor(ID))) %>%
  arrange(ID) -> df
saveRDS(df, file="data/microscopy.Rda")
#
# quick plot 
#
df %>%
  ggplot(aes(x = time, y = length, colour = factor(replicate), group = interaction(fungal_ic, replicate))) +
  geom_point(size=0.5) +
  geom_line(alpha=.5) +
  facet_grid(. ~ fungal_ic) +
  scale_y_log10() +
  theme_bw(base_size = 20) -> p

p <- p + xlab("Time [hrs]") +
  ylab("Length") +
  scale_color_manual(name = "Replicate", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))

p

