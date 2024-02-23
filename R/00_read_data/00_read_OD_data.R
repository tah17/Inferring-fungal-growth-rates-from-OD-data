#
# Script that reads in the spreadsheet data from the OD reader performed by Natasha Motsi.
#
#
rm(list = ls())
library(tidyverse)
library(readxl)
library(ggplot2)
library(reshape2)
#'
#' read_data
#'
#' Reads the OD data 
#'
#' @param excel_file String. The name of the excel spreadsheet the data holds.
#' @param sheet_name String. Sheet name of excel sheet that contains raw data from OD reader of a particular biological replicate.
#' @param data_range String. Range of cells in the excel sheet that contain the OD readings.
#' @param fungal_ic. List of doubles. The initial concentration of conidia of each set of OD readings over time. Must be in the same order as the data.
#' @param rep_labels. List of strings. Technical replicates of OD readings in the excel sheet. Must be in the same order as the data.
#'
#' @return A data frame of OD readings for different initial fungal concentrations.
#'
read_data <- function(excel_file, sheet_name, data_range, fungal_ic, rep_labels){
  #
  # create tibble that holds the fungal_ic and technical reps in the same order
  # as the OD data in the sheet.
  #
  tibble(fungal_ic) %>%
    add_column(technical_reps = rep_labels) %>%
    mutate(experiment = as.character(row_number())) -> col_names
  #
  # read in the OD values
  #
  df1 <- read_excel(excel_file, sheet = sheet_name, range = data_range)
  names(df1) <- c("time", "temp", col_names$experiment)  #rename the data to match the order as that in col_names
  #
  # match the  fungal_ic and technical reps to the correct experiment number
  # & place in long format
  #
  pivot_longer(df1, !c(time, temp), names_to = "experiment", values_to = "OD") %>%
    full_join(col_names) %>%
    select(-experiment) %>%
    drop_na() -> df2
  return(df2)
}

# Read in Bio Rep 1 -------------------------------------------------
data_range <- "B23:R49"
fungal_ic <- c(rep(c(2*10^c(5, 4, 3, 2), 0), each = 2), c(2*10^c(5, 2, 4), 0, 2*10^3))  ## in spores [N] (has to be converted to [N/ul] in models)
rep_labels <- c(rep(c("R1", "R2"), each = 1, times = 5), rep("R3", times = 5))

excel_file1 <- "data/OD.xlsx"
sheet_name1 <- "Plate\ 1\ -\ Data\ from\ OD\ reader"

read_data(excel_file1, sheet_name1, data_range, fungal_ic, rep_labels) %>%
  add_column(biological_reps = "R1") -> bio_rep1
#
# quick plot of the 1st biological rep
#
bio_rep1 %>%
  ggplot(aes(x = time, y = OD, group = technical_reps)) +
  geom_point(size=0.5) +
  geom_line(alpha=.5) +
  facet_grid(. ~ fungal_ic) +
  scale_y_log10() +
  theme_bw(base_size = 20)

# Read in Bio Rep 2 -------------------------------------------------
excel_file2 <- "data/OD\ 2\ 25.09.xlsx"
sheet_name2 <- "Plate\ 1\ -\ Sheet1"
read_data(excel_file2, sheet_name2, data_range, fungal_ic, rep_labels) %>%
  add_column(biological_reps = "R2") -> bio_rep2
#
# quick plot of the 2nd biological rep
#
bio_rep2 %>%
  ggplot(aes(x = time, y = OD, group = technical_reps)) +
  geom_point(size=0.5) +
  geom_line(alpha=.5) +
  facet_grid(. ~ fungal_ic) +
  scale_y_log10() +
  theme_bw(base_size = 20)

# Read in Bio Rep 3 -------------------------------------------------
sheet_name3 <- "Bio\ rep\ 3\ 28.09"
read_data(excel_file1, sheet_name3, data_range, fungal_ic, rep_labels) %>%
  add_column(biological_reps = "R3") -> bio_rep3
#
# quick plot of the 3rd biological rep
#
bio_rep3 %>%
  ggplot(aes(x = time, y = OD, group = technical_reps)) +
  geom_point(size=0.5) +
  geom_line(alpha=.5) +
  facet_grid(. ~ fungal_ic) +
  theme_bw(base_size = 20)

# Create and Save Entire Data Set -----------------------------------------
rbind(bio_rep1, bio_rep2, bio_rep3) %>%
  unite("ID", fungal_ic:biological_reps, remove = FALSE) %>%
  mutate(biological_reps = as.numeric(factor(biological_reps))) %>%
  mutate(blanks = fungal_ic == 0) %>%
  mutate(ID = as.numeric(factor(ID))) -> df
#
# add in NAs for missing data
#
df %>%
  add_row(time=rep(25, dim(filter(df, biological_reps==2, time==24))[1]),
          ID=filter(df, biological_reps==2, time==24)$ID,
          fungal_ic=filter(df, biological_reps==2, time==24)$fungal_ic,
          technical_reps=filter(df, biological_reps==2, time==24)$technical_reps,
          biological_reps=filter(df, biological_reps==2, time==24)$biological_reps,
          blanks=filter(df, biological_reps==2, time==24)$blanks) %>%
  mutate(missing_data = is.na(OD)) %>%
  arrange(ID) -> df_w_missing
#
# save data
#
saveRDS(df_w_missing, file="data/OD_all_ICs.Rda")
#
# plot all the data
#
df_w_missing %>%
  ggplot(aes(x = time, y = OD, colour = factor(biological_reps), group = interaction(technical_reps, biological_reps))) +
  geom_point(size=0.5) +
  geom_line(alpha=.5) +
  facet_grid(. ~ fungal_ic) +
  scale_y_log10() +
  theme_bw(base_size = 20)

#
# save only high IC version of data
#
df_w_missing %>%
  filter(fungal_ic%in%c(0, 2e5)) %>%
  select(-ID) %>%
  unite("ID", fungal_ic:biological_reps, remove = FALSE) %>%
  mutate(ID = as.numeric(factor(ID))) %>%
  saveRDS(file="data/OD_high_IC.Rda")

