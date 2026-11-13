# Interplay between climate, childhood mixing, and population-level susceptibility explains a sudden shift in RSV seasonality in Japan

[![DOI](https://zenodo.org/badge/742662924.svg)](https://doi.org/10.5281/zenodo.17396084)

------

This repository contains all code and data used for the analysis.

* `R`: contains an R script for simulating the SIRS model
* `analysis`: contains an R script for comparing changes in transmission rate before and after 2017
* `data`: contains time series data for RSV cases and number of children attending childcare facilities as well as data on population sizes by prefecture
* `doc`: contains latex files for the manuscript
* `figure`: contains R scripts for plotting all figures
  * `figure1.R`: R script for generating figure 1 and figure S1
  * `figure_comb_sirs_npi.R`: R script for generating figure 2
  * `figure_comb_sirs_change.R`: R script for generating figure 3
  * `figure_childcare.R`: R script for generating figure 4
  * `figure_map.R`: R script for generating figure S2
  * `figure_joint_climate.R`: R script for generating figure S3
  * `figure_ryukyu_sirs_change.R`: R script for generating figure S4
* `script`: contains an R script for processing data
* `simulate_sirs`: contains R scripts for simulating the fitted SIRS model across different islands to explore the relationship between seasonal forcing and center of gravity
  * `simulate_sirs_honshu.R`: R script for simulating epidemic dynamics in Honshu island
  * `simulate_sirs_kyushu.R`: R script for simulating epidemic dynamics in Kyushu island
  * `simulate_sirs_ryukyu.R`: R script for simulating epidemic dynamics in Ryukyu island
  * `simulate_sirs_shikoku.R`: R script for simulating epidemic dynamics in Shikoku island
  * `simulate_sirs_interpolate.R`: R script for simulating epidemic dynamics using interpolated transmission rates
* `stanfit_sirs`: contains R scripts for fitting SIRS model using Stan
* `stanfit_sirs2`: contains R scripts for fitting SIRS model that allows for changes in transmission using Stan
* `stanmodel`: constrains stan scripts for deterministic models 

------

* R scripts in `stanfit_sirs` and  `stanfit_sirs2` folders can be run independently as standalone files; these files need to be run first to generate rda files for model fits.
* R scripts in `simulate_sirs` folder can be run after all stan models have been fitted. These scripts will generate rda files that contain a summary of analyses of fitted models.
* R scripts in `figure` folder need to be run after all models have been fitted and analyzed.

------

All code was run on M2 MacBook Pro, 2023. Each model fit takes <10 minutes to run. All other code will take much less time than model fitting.
