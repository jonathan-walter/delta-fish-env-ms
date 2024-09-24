# delta-fish-env-ms
Derived data products and code accompanying Walter et al. (in review).

This repository reproduces analyses from the study ``Quantifying drivers of fish population change using multiple long-term surveys.''
It contains the directories code, data, and outputs.
Data used in this study have been archived in the Enviromental Data Initiative repository, and detailed metadata are provided there:

Walter, J.A. 2024. Age-0 fish abundances and putative environmental drivers for the Sacramento-San Joaquin Delta, California, 1980 to 2020. ver 1. Environmental Data Initiative. https://doi.org/10.6073/pasta/d2f8c04a81efc0cd3238dbc503d07036.

The code directory contains the files:
1) fn_modsel_regionwide.R: a helper function for running MARSS models and performing model selection based on AICc.
2) metaanalysis.R: this script synthesizes modeling results across many species and reproduces manuscript figures.
3) modelFitting.R: this script runs input data processing and statistical modeling over many species and outputs .RDS objects that are then used by metanalysis.R. It is computationally intensive to run due to the number of models and species and will take several days to run if run serially.

The data directory contains the following files:
1) Age0_thresh_fall_comb.csv: Age-0 maximum length thresholds for fall months for candidate species.
2) env_drivers_combined.csv: Annualized (water-year basis), centered, and scaled environmental driver time series.
3) fish_cleaned_fallSurveys.csv: Fish abundance data.
4) SpeciesDetails.csv: Basic ecological characteristics and taxonomy of candidate species.

The outputs directory is designated to hold output .RDS files from modelFitting.R and figures produced by metaanalysis.R. Before running these scripts, it contains only a placeholder text file so that the directory is recognized by Git.

The code was written in R version 4.4.0/R Studio version 2024.04.1+748 and uses the packages MARSS v3.11.9, trend v1.1.6, viridis v0.6.5, tidyverse v2.0.0, lubridate v1.9.3, mgcv v1.9.1, fields v15.2.
