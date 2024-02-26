# Inferring-fungal-growth-rates-from-OD-data

## Table of contents
* [Introduction](#introduction)
* [Technologies](#technologies)
* [Setup](#setup)
* [Usage](#usage)

## Introduction
The repository holds the code to infer growth rates from optical density (OD) fungal growth data.

## Technologies
The code is written in R (v4.2.0) and Stan. Details of the packages and their versions can be found in the [renv.lock](renv.lock) file.

## Setup
First either clone or download the repository to your machine. The package renv can be used to download the correct versions of the R packages used in this project. As detailed in the [Introduction to renv](https://rstudio.github.io/renv/articles/renv.html) vignette, once the repository is opened the appropriate version of renv will be automatically installed. To then download and install the required packages run:

```
renv::restore()
```

in the R console when prompted.

## Usage

The code is intended to be run by using the .R scripts in the [`R/`](R) folder:
* [`00_read_data/`](R/00_read_data/): Each script in this folder takes an excel spreadsheet of data from an experiment and stores the processed data as a tibble in a "data" folder.
* [`01_prior_pred.R`](R/01_prior_pred.R): For prior predictive and fake data checks of the models.  
* [`02_post_pred.R`](R/02_post_pred.R): Samples from the posterior of the model when using all the data available. Performs a posterior predictive check.
* [`03_validate/`](R/03_validate/): Scripts to conduct k-fold cross validation.
* [`04_plot/`](R/04_plot/): Holds the plotting scripts used to generate each of the plots in the manuscript.

The above scripts are run for each of the Stan models in the [`models/`](models) folder. The [`models/`](models) folder is structured by having a file that details the probabilistic model in Stan, `model.stan`, and a corresponding meta-file `model.data.R`, which is a script that loads in the arguments and data corresponding to the model in `model.stan`.
