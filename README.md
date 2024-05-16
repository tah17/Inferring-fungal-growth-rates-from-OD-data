# Inferring-fungal-growth-rates-from-OD-data

## Table of contents
* [Introduction](#introduction)
* [Technologies](#technologies)
* [Setup](#setup)
* [Usage](#usage)
* [License](#license)

## Introduction
The repository holds the code for the manuscript [Hameed _et al._ (2024), "Inferring fungal growth rates from optical density data"](https://doi.org/10.1371/journal.pcbi.1012105).

## Technologies
The code is written in R (v4.2.0) and Stan. Details of the packages and their versions can be found in the [renv.lock](renv.lock) file.

## Setup
First either clone or download the repository to your machine. The package renv can be used to download the correct versions of the R packages used in this project. As detailed in the [Introduction to renv](https://rstudio.github.io/renv/articles/renv.html) vignette, once the repository is opened the appropriate version of renv will be automatically installed. To then download and install the required packages run:

```
renv::restore()
```

in the R console when prompted.

The scripts in [`03_validate/`](R/03_validate/) are intended to be run using Imperial College London's High Performance Computing (HPC) services. To run these scripts on Imperial's HPC facility, a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) environment with R at v4.2.0 should be created:

```
module load anaconda3/personal
conda create -n r420 r-base=4.2.0 -c conda-forge
source activate r420
```

and the correct versions of the R packages can then downloaded in the same manner as detailed above.

## Usage

The code is intended to be run by using the .R scripts in the [`R/`](R) folder:
* [`00_read_data/`](R/00_read_data/): Each script in this folder takes an excel spreadsheet of data from an experiment and stores the processed data as a tibble in a "data" folder.
* [`01_prior_pred.R`](R/01_prior_pred.R): For prior predictive and fake data checks of the models.  
* [`02_post_pred.R`](R/02_post_pred.R): Samples from the posterior of the model when using all the data available. Performs a posterior predictive check.
* [`03_validate/`](R/03_validate/): Scripts to conduct k-fold cross validation.
* [`04_plot/`](R/04_plot/): Holds the plotting scripts used to generate each of the plots in the manuscript.

The above scripts are run for each of the Stan models in the [`models/`](models) folder. The [`models/`](models) folder is structured by having a file that details the probabilistic model in Stan, `model.stan`, and a corresponding meta-file `model.data.R`, which is a script that loads in the arguments and data corresponding to the model in `model.stan`.

## License
Licensed under the GPLv3 license. See [LICENSE](LICENSE) for more information.
