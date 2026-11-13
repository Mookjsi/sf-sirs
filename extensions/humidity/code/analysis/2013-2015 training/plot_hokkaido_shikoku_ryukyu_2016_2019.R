# Plot saved predictions for Hokkaido, Shikoku, and Ryukyu.
# Training: 2013-2015. Prediction: 2016-2019. Open and click Source.
# No Stan compilation or MCMC is run.
project_root <- Sys.getenv("PERTURBATION_ROOT")
if (dir.exists(file.path(project_root, ".R-library"))) {
  .libPaths(c(file.path(project_root, ".R-library"), .libPaths()))
}
library(ggplot2)

train_years <- 2013:2015
test_years <- 2016:2019
island_folders <- c(
  Hokkaido = "humidity_seasonal_hokkaido_train2013_2015",
  Shikoku = "humidity_seasonal_shikoku_train2013_2015",
  Ryukyu = "humidity_seasonal_ryukyu_train2013_2015"
)
output_dir <- paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2015 training/plots_hokkaido_shikoku_ryukyu")

# Check that the saved fits and predictions use the intended dates and inputs.
weather <- read.csv(file.path(project_root, "data_processed", "data_processed_humidity.csv"),
                    encoding = "UTF-8")
population <- read.csv(file.path(project_root, "data", "japan_pop.csv"), encoding = "UTF-8")
population$prefecture <- sub("-.*", "", population$prefecture)
expected_prefecture_count <- c(Hokkaido = 1L, Shikoku = 4L, Ryukyu = 1L)
curves <- vector("list", length(island_folders))
diagnostics <- vector("list", length(island_folders))

for (i in seq_along(island_folders)) {
  island_id <- names(island_folders)[i]
  result_dir <- file.path(project_root, "analysis", island_folders[[i]])
  result_file <- file.path(result_dir, "humidity_seasonal_result.rds")
  fit_file <- file.path(result_dir, "humidity_seasonal_fit_cache.rds")
  if (!file.exists(result_file) || !file.exists(fit_file)) {
    stop("Saved result/fit missing: ", island_id, ". Checked: ", result_dir)
  }
  res <- readRDS(result_file)
  fit <- readRDS(fit_file)
  if (!isTRUE(res$diagnostics$summary$diagnostics_ok) ||
      !isTRUE(fit$diagnostics$summary$diagnostics_ok)) {
    stop("MCMC diagnostics failed: ", island_id, ". No new plot was generated.")
  }
  if (!isTRUE(all.equal(fit$standata, res$standata)) ||
      !isTRUE(all.equal(fit$diagnostics, res$diagnostics))) {
    stop("Fit cache and saved prediction do not match: ", island_id)
  }
  curve <- as.data.frame(res$result$curve)
  expected_year <- rep(test_years, each = 52L)
  expected_week <- rep(1:52, length(test_years))
  stopifnot(
    res$standata$N == 52L * length(train_years),
    identical(as.integer(res$train$year), rep(train_years, each = 52L)),
    identical(as.integer(res$train$week), rep(1:52, length(train_years))),
    identical(as.integer(res$test$year), expected_year),
    identical(as.integer(res$test$week), expected_week),
    identical(as.integer(curve$year), expected_year),
    identical(as.integer(curve$week), expected_week),
    identical(as.integer(curve$step), seq_along(expected_year)),
    isTRUE(all.equal(as.numeric(curve$observed), as.numeric(res$test$cases))),
    all(is.finite(as.matrix(curve[c("observed", "pred_median", "lower", "upper")]))),
    all(curve$lower <= curve$pred_median & curve$pred_median <= curve$upper)
  )
  saved <- rbind(as.data.frame(res$train), as.data.frame(res$test))
  prefectures <- population$prefecture[population$island == island_id]
  stopifnot(length(prefectures) == expected_prefecture_count[[island_id]])
  x <- weather[weather$key %in% prefectures &
                 weather$year %in% c(train_years, test_years), ]
  stopifnot(nrow(x) == nrow(saved) * length(prefectures),
            !anyNA(x[c("key", "year", "week", "humidity")]),
            !anyDuplicated(x[c("key", "year", "week")]),
            all(is.finite(x$humidity)))
  counts <- aggregate(humidity ~ year + week, x, length)
  stopifnot(nrow(counts) == nrow(saved), all(counts$humidity == length(prefectures)))
  x <- aggregate(humidity ~ year + week, x, mean)
  index <- match(paste(saved$year, saved$week), paste(x$year, x$week))
  if (!isTRUE(all.equal(as.numeric(saved$q_gkg),
                        as.numeric(x$humidity[index] * 1000), tolerance = 1e-8))) {
    stop("Saved humidity does not match the complete prefecture means: ", island_id)
  }
  curve$island <- island_id
  curves[[i]] <- curve
  diagnostics[[i]] <- cbind(island = island_id, fit$diagnostics$summary)
}

prediction_three_islands <- do.call(rbind, curves)
prediction_three_islands$island <- factor(prediction_three_islands$island,
                                          levels = names(island_folders))
diagnostics_three_islands <- do.call(rbind, diagnostics)
print(diagnostics_three_islands, row.names = FALSE, digits = 4)
n_test <- 52L * length(test_years)
if (length(test_years) == 2L) {
  x_breaks <- c(1, 27, 53, 79, n_test)
  x_labels <- c(paste(test_years[1], "W1"), paste(test_years[1], "W27"),
                paste(test_years[2], "W1"), paste(test_years[2], "W27"), paste(test_years[2], "W52"))
} else {
  x_breaks <- c(1L + 52L * seq.int(0L, length(test_years) - 1L), n_test)
  x_labels <- c(paste(test_years, "W1"), paste(max(test_years), "W52"))
}

g_three_islands_2013_2015 <- ggplot(prediction_three_islands, aes(x = step)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#90CAF9", alpha = 0.45) +
  geom_line(aes(y = observed), colour = "black", linewidth = 0.45) +
  geom_line(aes(y = pred_median), colour = "#0072B2", linewidth = 0.7) +
  geom_vline(xintercept = 52L * seq_len(length(test_years) - 1L) + 0.5,
             linetype = "dotted", linewidth = 0.4) +
  facet_wrap(~island, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels,
                     expand = expansion(mult = c(0.025, 0.055))) +
  labs(title = "humidity_seasonal: Hokkaido, Shikoku and Ryukyu",
       subtitle = sprintf("Conditional hindcast, %d-%d (trained on %d-%d)",
                          min(test_years), max(test_years), min(train_years), max(train_years)),
       x = NULL, y = "Weekly reported cases",
       caption = "Black: observed; blue: predictive median; shading: 95% predictive interval.") +
  theme_bw(base_size = 14) +
  theme(strip.text = element_text(size = 15, face = "bold"),
        panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0, size = 11))

print(g_three_islands_2013_2015)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png_file <- file.path(output_dir, "hokkaido_shikoku_ryukyu_prediction_2016_2019.png")
ggsave(png_file, g_three_islands_2013_2015, width = 12, height = 10,
       units = "in", dpi = 220, bg = "white")
message("Saved prediction figure: ", png_file)
