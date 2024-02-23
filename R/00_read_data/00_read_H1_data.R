rm(list = ls())
library(tidyverse)
library(readxl)
library(ggplot2)
library(data.table)
#'
#' read_data
#'
#' Reads the nuclear count data 
#'
#' @param idx Int. Index of ID of data.
#' @param excel_file String. The name of the excel spreadsheet the data holds.
#' @param data_range String. Range of cells in the excel sheet that contain the ID's data.
#'
#' @return A data frame of nuclear count readings
#'
read_data <- function(x, excel_file, data_range){
  read_excel(excel_file, range=data_range) %>%
    drop_na() %>%
    rename(time = `Time (hrs`) %>%
    add_column(fungal_ic=1e1) %>%   ## IC in [N/ul]
    add_column(ID = x) -> df
  return(df)
}

excel_file <- "data/H1.xlsx"
data_ranges <- c("A1:B10", "A13:B25", "A29:B39", "A45:B57", "A61:B73", "A78:B91", "A96:B110")
lapply(1:length(data_ranges), function(x) read_data(x, excel_file, data_ranges[x])) %>% 
  rbindlist() %>%
  complete(time, ID) %>% 
  fill(fungal_ic) %>%
  arrange(ID) -> df

saveRDS(df, file ="data/h1_data.Rda")

df %>%
  ggplot(aes(x = time, y = nuclei, group = ID)) +
  geom_point(size=0.5) +
  geom_line(alpha=.5) +
  scale_y_log10() +
  theme_bw(base_size = 20) -> p

p <- p + xlab("Time [hrs]") +
  ylab("Nuclei")

p
