rm(list = ls())
library(reshape2)
library(data.table)
library(tidyverse)
library(caret)
seed <- 683861
set.seed(seed)
#
# Read in command line arguments
#
myargs = commandArgs(trailingOnly=TRUE)
IC_used <- myargs[1] # ICs used - all_ICs or high_IC
if (length(myargs) < 1) { 
  stop("Error: missing command line argument 'IC_used'.")
} else if (length(myargs) < 2) {
  fold_no <- 5
  warning("Using default Number of Folds")
} else {
  fold_no <- as.integer(myargs[2])
}
if ((!identical(IC_used, "all_ICs"))&&(!identical(IC_used, "high_IC"))) {
  stop(paste("Requested ", IC_used, ", which is not a valid option \nPlease choose either: 'all_ICs' or 'high_IC'"))
}
# Read in Data ------------------------------------------------------------
#
df <- readRDS(file=paste("data/OD_", IC_used, ".Rda", sep=""))
#
# Create idx lists
#
no_of_folds <- fold_no
train_idxs <- createFolds(unique(df$ID), k = no_of_folds, returnTrain = TRUE)  #createFolds from caret library
cv_idxs <- list(train_idxs = unname(train_idxs))
#
# Save lists
#
filename <- paste("data/split_idxs_", fold_no, "_fold_", IC_used, ".Rda", sep="")
saveRDS(cv_idxs, file=filename)
