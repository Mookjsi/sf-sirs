# Plot saved 2018-2019 forecast metrics only. No MCMC or model fitting.
project_root <- Sys.getenv("PERTURBATION_ROOT")
.libPaths(c(file.path(project_root, ".R-library"), .libPaths()))
library(ggplot2)

metrics_file <- paste0(
  paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training/"),
  "plots_baseline_npieff_vs_humidity_no_npieff/comparison_metrics.csv"
)
metrics <- read.csv(metrics_file, check.names = FALSE)
metric_keys <- c("mean_log_score", "MAE", "RMSE_mean", "coverage_95", "mean_interval_width")
stopifnot(all(c("island", "model", "period", "n_weeks", metric_keys) %in% names(metrics)),
          nrow(metrics) == 6L, all(metrics$period == "2018-2019"),
          all(metrics$n_weeks == 104L), !anyDuplicated(metrics[c("island", "model")]),
          all(is.finite(as.matrix(metrics[metric_keys]))),
          all(metrics$coverage_95 >= 0 & metrics$coverage_95 <= 1))

island_levels <- c("Hokkaido", "Shikoku", "Ryukyu")
model_levels <- c("Baseline + npieff", "Humidity seasonal (no npieff)")
metric_labels <- c(
  "\ud3c9\uade0 \ub85c\uadf8\uc608\uce21\uc810\uc218\n\ud074\uc218\ub85d \uc88b\uc74c (0\uc5d0 \uac00\uae4c\uc6b8\uc218\ub85d)",
  "MAE \u00b7 \uc2e0\uace0 \uc218\n\uc791\uc744\uc218\ub85d \uc88b\uc74c",
  "RMSE \u00b7 \uc2e0\uace0 \uc218\n\uc791\uc744\uc218\ub85d \uc88b\uc74c",
  "95% \uc608\uce21\uad6c\uac04 \ud3ec\ud568\ub960 \u00b7 %\n\uc810\uc120\uc740 \ubaa9\ud45c \ud3ec\ud568\ub960 95%",
  "\ud3c9\uade0 \uc608\uce21\uad6c\uac04 \ud3ed \u00b7 \uc2e0\uace0 \uc218\n\ud3ec\ud568\ub960\uacfc \ud568\uaed8 \ud574\uc11d"
)
bar_data <- do.call(rbind, lapply(seq_along(metric_keys), function(j) {
  v <- metrics[[metric_keys[j]]]
  if (metric_keys[j] == "coverage_95") v <- 100 * v
  data.frame(
    island = metrics$island, model = metrics$model,
    metric = metric_labels[j], value = v,
    value_label = if (metric_keys[j] == "mean_log_score") sprintf("%.2f", v) else sprintf("%.1f", v),
    label_vjust = ifelse(metric_keys[j] == "coverage_95" & v >= 88, 1.5,
                         ifelse(v < 0, 1.3, -0.45)),
    label_colour = ifelse(metric_keys[j] == "coverage_95" & v >= 88, "white", "#111111")
  )
}))
bar_data$island <- factor(bar_data$island, levels = island_levels)
bar_data$model <- factor(bar_data$model, levels = model_levels)
bar_data$metric <- factor(bar_data$metric, levels = metric_labels)
stopifnot(!anyNA(bar_data))
coverage_reference <- data.frame(metric = factor(metric_labels[4], levels = metric_labels), target = 95)
coverage_limit <- data.frame(
  metric = factor(metric_labels[4], levels = metric_labels),
  island = factor("Hokkaido", levels = island_levels), value = 100
)

g_metrics <- ggplot(bar_data, aes(x = island, y = value, fill = model)) +
  geom_hline(yintercept = 0, colour = "#657183", linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(data = coverage_reference, aes(yintercept = target),
             inherit.aes = FALSE, linetype = "dashed", colour = "#505968", linewidth = 0.6) +
  geom_blank(data = coverage_limit, aes(x = island, y = value), inherit.aes = FALSE) +
  geom_text(aes(label = value_label, vjust = label_vjust, colour = label_colour, group = model),
            position = position_dodge(width = 0.75), size = 3.8, family = "Malgun Gothic") +
  facet_wrap(~ metric, ncol = 2, scales = "free") +
  scale_colour_identity() +
  scale_x_discrete(labels = c(Hokkaido = "\ud64b\uce74\uc774\ub3c4", Shikoku = "\uc2dc\ucf54\ucfe0", Ryukyu = "\ub958\ud050")) +
  scale_y_continuous(expand = expansion(mult = c(0.13, 0.14))) +
  scale_fill_manual(
    values = c("Baseline + npieff" = "#2873AD", "Humidity seasonal (no npieff)" = "#D77A29"),
    labels = c("Baseline (npieff \ud3ec\ud568)", "\ube44\uc2b5 \ubaa8\ud615 (npieff \uc81c\uc678)")
  ) +
  labs(
    title = "Baseline\uacfc \ube44\uc2b5 \ubaa8\ud615\uc758 \uc608\uce21 \uc131\ub2a5",
    subtitle = "\ud559\uc2b5: 2013\u20132017\ub144  |  \ud3c9\uac00: 2018\u20132019\ub144, \uc12c\ubcc4 104\uc8fc",
    x = NULL, y = NULL, fill = NULL,
    caption = paste(
      "MAE\ub294 \uc608\uce21\uc911\uc559\uac12, RMSE\ub294 \uc608\uce21\ud3c9\uade0 \uae30\uc900. \uac19\uc740 \uc12c \uc548\uc5d0\uc11c \ub450 \ubaa8\ud615\uc744 \ube44\uad50\ud569\ub2c8\ub2e4.",
      "\uad6c\uac04 \ud3ed\uc774 \uc881\uc544\ub3c4 \ud3ec\ud568\ub960\uc774 \ub0ae\uc73c\uba74 \uc88b\uc740 \uc608\uce21\uc774\ub77c\uace0 \ubcfc \uc218 \uc5c6\uc2b5\ub2c8\ub2e4. \ub450 \ubaa8\ud615\uc740 npieff \ud3ec\ud568 \uc5ec\ubd80\ub3c4 \ub2e4\ub985\ub2c8\ub2e4.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 13, base_family = "Malgun Gothic") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 23, face = "bold", colour = "#192B3D", margin = margin(b = 8)),
    plot.subtitle = element_text(size = 13, colour = "#526171", margin = margin(b = 12)),
    plot.caption = element_text(hjust = 0, size = 11, colour = "#526171", lineheight = 1.3, margin = margin(t = 18)),
    legend.position = "top", legend.justification = "left", legend.text = element_text(size = 12),
    strip.text = element_text(hjust = 0, face = "bold", size = 13, lineheight = 1.25, margin = margin(10, 8, 10, 8)),
    strip.background = element_rect(fill = "#EEF2F6", colour = NA),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E5EAF0", linewidth = 0.3),
    axis.text = element_text(colour = "#354555", size = 11),
    panel.spacing = grid::unit(1.1, "lines"),
    plot.margin = margin(20, 22, 18, 20)
  )

output_dir <- getOption("forecast.metrics.output_dir", dirname(metrics_file))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "forecast_metrics_bars_2018_2019.png"),
       plot = g_metrics, width = 12, height = 12, units = "in", dpi = 160, type = "cairo")
ggsave(file.path(output_dir, "forecast_metrics_bars_2018_2019.pdf"),
       plot = g_metrics, width = 12, height = 12, units = "in", device = cairo_pdf)
if (interactive()) print(g_metrics)
message("Saved metric bar charts in: ", output_dir)
