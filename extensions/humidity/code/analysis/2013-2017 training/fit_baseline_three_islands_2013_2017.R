project_root <- Sys.getenv("PERTURBATION_ROOT")
support_dir <- paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/npieff_support")
islands <- "Hokkaido"
models <- "baseline"
validation_only <- FALSE
force_refit <- FALSE

sampling_settings <- list(
  chains = 4L, cores = 4L,
  iter = 16000L,
  warmup = 8000L,
  seed = 101L,
  adapt_delta = 0.99,
  max_treedepth = 15L
)

source(file.path(support_dir, "npieff_prediction_helpers.R"), encoding = "UTF-8")
source(file.path(support_dir, "npieff_workflow.R"), encoding = "UTF-8")
stan_file <- file.path(project_root, "stanmodel", "sirs_npieff_holdout.stan")
output_root <- file.path(project_root, "analysis",
                         "npieff_baseline_humidity_train2013_2017_forecast2018_2019")
npieff_run <- np_run(
  project_root, stan_file, output_root, islands, models, sampling_settings,
  validation_only = validation_only, force_refit = force_refit
)

