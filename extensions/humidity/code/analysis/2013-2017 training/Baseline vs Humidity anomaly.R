# Baseline retains npieff; the seasonal-humidity model excludes npieff.
# This file does NOT run MCMC or relax any diagnostic threshold.
g_prediction <- NULL
comparison_metrics <- comparison_diagnostics <- prediction_plot_data <- comparison_status <- NULL
comparison_metrics_by_year <- NULL
project_root <- Sys.getenv("PERTURBATION_ROOT")
support_dir <- paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/npieff_support")
.libPaths(c(file.path(project_root, ".R-library"), .libPaths()))
library(ggplot2)
source(file.path(support_dir, "npieff_prediction_helpers.R"), encoding = "UTF-8")
source(file.path(support_dir, "npieff_workflow.R"), encoding = "UTF-8")
source(file.path(support_dir, "humidity_no_npieff_support.R"), encoding = "UTF-8")
humidity_script_dir <- paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training")
stan_file <- file.path(project_root, "stanmodel", "sirs_npieff_holdout.stan")
output_root <- file.path(project_root, "analysis",
                         "npieff_baseline_humidity_train2013_2017_forecast2018_2019")
islands <- c("Hokkaido", "Shikoku", "Ryukyu")
models <- c("baseline", "humidity_seasonal")
model_labels <- c(baseline = "Baseline + npieff",
                  humidity_seasonal = "Humidity seasonal (no npieff)")
valid <- statuses <- list()
for (island in islands) for (model in models) {
  key <- paste(island, model, sep = "/")
  checked <- tryCatch(
    if (model == "baseline") {
      np_collect(project_root, stan_file, output_root, island, model)
    } else hn_collect(project_root, island, humidity_script_dir),
    error = function(e) e
  )
  ready <- !inherits(checked, "error")
  reason <- if (ready) "" else conditionMessage(checked)
  state <- if (ready) "ready" else if (grepl("diagnostics_failed|MCMC diagnostics failed", reason)) {
    "diagnostics_failed"
  } else if (grepl("Missing .*fit/predictions", reason)) {
    "missing_result"
  } else "validation_error"
  statuses[[key]] <- data.frame(island = island, model_key = model,
    model = unname(model_labels[model]), status = state, reason = reason)
  if (ready) {
    checked$curves$model <- checked$metrics$model <- checked$diagnostics$model <- unname(model_labels[model])
    valid[[key]] <- checked
  }
}
comparison_status <- do.call(rbind, statuses)
rownames(comparison_status) <- NULL
print(comparison_status[c("island", "model_key", "status")], row.names = FALSE)
if (!length(valid)) {
  stop("No validated predictions are available. Read comparison_status$reason. No MCMC was run.",
       call. = FALSE)
}
prediction_plot_data <- do.call(rbind, lapply(valid, function(x) x$curves))

# Score the same saved weekly predictive distributions, without fitting again.
# MAE uses the predictive MEDIAN (the blue line); RMSE_mean uses the predictive MEAN.
# log_score is the log posterior-predictive probability mass of the observed count.
# Coverage is a fraction, with nominal target 0.95. Interpret width WITH coverage.
forecast_metrics_summary <- function(d, period) {
  required <- c("observed", "pred_median", "pred_mean", "lower", "upper", "log_score")
  stopifnot(nrow(d) > 0L, all(required %in% names(d)),
            all(is.finite(as.matrix(d[setdiff(required, "log_score")]))),
            all(d$lower <= d$upper), !anyNA(d$log_score),
            !any(d$log_score == Inf))
  data.frame(
    island = as.character(d$island[1]), model = as.character(d$model[1]),
    period = period, n_weeks = nrow(d),
    mean_log_score = mean(d$log_score),
    MAE = mean(abs(d$observed - d$pred_median)),
    RMSE_mean = sqrt(mean((d$observed - d$pred_mean)^2)),
    coverage_95 = mean(d$observed >= d$lower & d$observed <= d$upper),
    mean_interval_width = mean(d$upper - d$lower)
  )
}
overall_metrics <- yearly_metrics <- list()
observations_by_island <- list()
for (key in names(valid)) {
  d <- valid[[key]]$curves
  d <- d[order(d$year, d$week), , drop = FALSE]
  stopifnot(nrow(d) == 104L,
            all(d$year == rep(2018:2019, each = 52L)),
            all(d$week == rep(1:52, 2L)))
  island_key <- as.character(d$island[1])
  obs <- d[c("year", "week", "observed")]
  rownames(obs) <- NULL
  if (!is.null(observations_by_island[[island_key]])) {
    stopifnot(isTRUE(all.equal(obs, observations_by_island[[island_key]],
                              check.attributes = FALSE)))
  }
  observations_by_island[[island_key]] <- obs
  overall_metrics[[key]] <- forecast_metrics_summary(d, "2018-2019")
  for (yr in 2018:2019) {
    yearly_metrics[[paste(key, yr, sep = "/")]] <-
      forecast_metrics_summary(d[d$year == yr, , drop = FALSE], as.character(yr))
  }
}
comparison_metrics <- do.call(rbind, overall_metrics)
comparison_metrics_by_year <- do.call(rbind, yearly_metrics)
rownames(comparison_metrics) <- rownames(comparison_metrics_by_year) <- NULL
comparison_diagnostics <- do.call(rbind, lapply(valid, function(x) x$diagnostics))
comparison <- list(curves = prediction_plot_data, metrics = comparison_metrics,
                   metrics_by_year = comparison_metrics_by_year,
                   diagnostics = comparison_diagnostics, status = comparison_status)
g_prediction <- hn_plot_comparison(prediction_plot_data)
unavailable <- comparison_status[comparison_status$status != "ready", , drop = FALSE]
if (nrow(unavailable)) {
  # These grey lines are observed counts ONLY; no failed or unvalidated predictions are used.
  observed <- np_read_data(project_root)$cases
  missing_observed <- do.call(rbind, lapply(seq_len(nrow(unavailable)), function(j) {
    u <- unavailable[j, ]
    d <- observed[observed$island == u$island & observed$year %in% 2018:2019, ]
    d <- d[order(d$year, d$week), ]
    stopifnot(nrow(d) == 104L)
    data.frame(island = u$island, model = u$model,
               step = seq_len(nrow(d)), observed = d$cases)
  }))
  unavailable$panel_label <- ifelse(unavailable$status == "diagnostics_failed",
    "No prediction: MCMC diagnostics failed",
    ifelse(unavailable$status == "missing_result",
           "No prediction: new fit/result not saved",
           "No prediction: saved-input validation failed"))
  missing_observed$island <- factor(missing_observed$island, levels = islands)
  missing_observed$model <- factor(missing_observed$model, levels = unname(model_labels))
  unavailable$island <- factor(unavailable$island, levels = islands)
  unavailable$model <- factor(unavailable$model, levels = unname(model_labels))
  g_prediction <- g_prediction +
    geom_line(data = missing_observed, aes(x = step, y = observed),
              inherit.aes = FALSE, colour = "grey60", linewidth = 0.4) +
    geom_label(data = unavailable, aes(x = 52.5, y = Inf, label = panel_label),
               inherit.aes = FALSE, vjust = 1.2, size = 3.2,
               colour = "#9C3B10", fill = "white") +
    labs(caption = paste(
      "Black/grey: observed; blue: predictive median; shading: 95% predictive interval.",
      "Grey-only panels have no validated prediction. No failed fits were plotted.", sep = "\n"))
}
g_prediction <- g_prediction +
  labs(subtitle = sprintf("Training: 2013-2017; %d/6 prediction panels ready; 2018-2019 conditional hindcast",
                          length(valid)))
print(comparison_metrics, row.names = FALSE, digits = 4)
message("Year-specific scores: print(comparison_metrics_by_year, row.names = FALSE, digits = 4)")
print(g_prediction)

# An option allows a local preview without writing to the Desktop.
output_dir <- getOption("npieff.comparison_output_dir",
  paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "baseline_npieff_vs_humidity_no_npieff_2018_2019.png"),
       plot = g_prediction, width = 14, height = 10, units = "in", dpi = 200)
write.csv(comparison_metrics, file.path(output_dir, "comparison_metrics.csv"), row.names = FALSE)
write.csv(comparison_metrics_by_year, file.path(output_dir, "comparison_metrics_by_year.csv"), row.names = FALSE)
write.csv(comparison_diagnostics, file.path(output_dir, "comparison_diagnostics.csv"), row.names = FALSE)
write.csv(prediction_plot_data, file.path(output_dir, "prediction_plot_data.csv"), row.names = FALSE)
write.csv(comparison_status, file.path(output_dir, "comparison_status.csv"), row.names = FALSE)
message("Saved available prediction panels and status table in: ", output_dir)

