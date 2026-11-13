# Baseline keeps npieff; humidity_seasonal excludes npieff.
# Source runs each fitting script; matching cached posteriors are reused.
# Training 2013-2017 / forecast 2018-2019. No automatic diagnostic retries.
# For checks without fitting/writing: options(mixed_no_npieff.validation_only = TRUE)
mixed_validation_only <- isTRUE(getOption("mixed_no_npieff.validation_only", FALSE))
mixed_old_options <- options(
  npieff.validation_only = mixed_validation_only,
  humidity_no_npieff.validation_only = mixed_validation_only
)
tryCatch({
  source(paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training/fit_baseline_three_islands_2013_2017.R"), encoding = "UTF-8")
  for (humidity_script_name in c("humidity_seasonal_hokkaido.R",
                                "humidity_seasonal_sikoku.R", "humidity_seasonal_ryukyu.R")) {
    tryCatch(
      source(file.path(paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training"), humidity_script_name),
             encoding = "UTF-8"),
      error = function(e) message(humidity_script_name, ": ", conditionMessage(e))
    )
  }
  if (!mixed_validation_only)
    source(paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training/Baseline vs Humidity anomaly.R"), encoding = "UTF-8")
}, finally = options(mixed_old_options))

