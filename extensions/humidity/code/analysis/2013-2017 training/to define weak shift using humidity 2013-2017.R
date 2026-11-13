# Hokkaido, Shikoku and Ryukyu peak shifts: train 2013-2017, compare 2018-2019.
# Run this file directly. It reads saved posteriors; it never runs MCMC.
# A: actual humidity. B: each island's 2013-2017 same-week mean humidity.
library(dplyr)
library(ggplot2)

project_root <- Sys.getenv("PERTURBATION_ROOT")
combined_output_dir <- getOption("humidity.peak_shift.output_dir", file.path(
  project_root, "analysis", "humidity_seasonal_three_islands_train2013_2017"
))
train_years <- 2013:2017
test_years <- 2018:2019
island_folders <- c(
  Hokkaido = "humidity_seasonal_hokkaido",
  Shikoku = "humidity_seasonal_shikoku",
  Ryukyu = "humidity_seasonal_ryukyu"
)

# Same deterministic SIRS prediction function used by the three fitted models.


hs_predict <- function(draws, week, standata, X_test) {
  stopifnot(is.matrix(draws), nrow(draws) > 0L,
            is.matrix(X_test), is.numeric(X_test),
            length(week) == nrow(X_test), length(week) > 0L,
            ncol(X_test) == standata$K,
            all(is.finite(X_test)), all(week %in% 1:52),
            all(week == as.integer(week)),
            length(standata$N) == 1L, standata$N >= 2L,
            length(standata$pop) == 1L, is.finite(standata$pop),
            standata$pop > 0,
            length(standata$mu) == 1L, is.finite(standata$mu),
            standata$mu > 0,
            length(standata$gamma) == 1L, is.finite(standata$gamma),
            standata$gamma > 0)

  beta_names <- paste0("beta[", 1:52, "]")
  theta_names <- paste0("theta[", seq_len(standata$K), "]")
  state_names <- paste0(c("S", "I", "R"), "[", standata$N, "]")
  required <- c(beta_names, theta_names, state_names,
                "rho", "tau", "omega")
  missing <- setdiff(required, colnames(draws))
  if (length(missing)) {
    stop("Missing prediction columns: ", paste(missing, collapse = ", "))
  }
  stopifnot(all(is.finite(draws[, required, drop = FALSE])))

  beta <- draws[, beta_names, drop = FALSE]
  theta <- draws[, theta_names, drop = FALSE]
  S <- draws[, state_names[1]]
  I <- draws[, state_names[2]]
  R <- draws[, state_names[3]]
  rho <- draws[, "rho"]
  tau <- draws[, "tau"]
  omega <- draws[, "omega"]
  stopifnot(all(beta > 0), all(S >= 0), all(I >= 0), all(R >= 0),
            all(rho >= 0 & rho <= 1), all(tau > 0), all(omega >= 0))

  delta <- 1 / tau
  mu <- standata$mu
  gamma <- standata$gamma
  pop <- standata$pop
  C <- matrix(NA_real_, nrow(draws), length(week))

  for (t in seq_along(week)) {
    log_multiplier <- as.vector(theta %*% X_test[t, ])
    beta_t <- beta[, week[t]] * exp(log_multiplier)
    foi <- beta_t * (I + omega) / pop
    if (any(!is.finite(foi)) || any(foi < 0)) {
      stop("Non-finite/negative force of infection at holdout step ", t,
           ". No invalid draws were removed or clamped.")
    }

    Sout <- -expm1(-(foi + mu)) * S
    StoI <- foi / (foi + mu) * Sout
    Iout <- -expm1(-(gamma + mu)) * I
    ItoR <- gamma / (gamma + mu) * Iout
    Rout <- -expm1(-(delta + mu)) * R
    RtoS <- delta / (delta + mu) * Rout
    S <- S + RtoS - Sout + mu * pop
    I <- I + StoI - Iout
    R <- R + ItoR - Rout
    C[, t] <- StoI * rho
    if (any(!is.finite(c(S, I, R, C[, t]))) ||
        any(c(S, I, R, C[, t]) < 0)) {
      stop("Invalid state or expected cases at holdout step ", t,
           ". No invalid draws were removed or clamped.")
    }
  }
  C
}


peak_info <- function(M, weeks) {
  stopifnot(is.matrix(M), ncol(M) == length(weeks),
            all(is.finite(M)), all(M >= 0))
  k <- apply(M, 1, which.max)
  maximum <- M[cbind(seq_len(nrow(M)), k)]
  data.frame(
    peak_week = weeks[k],
    boundary = k == 1L | k == ncol(M),
    tied = rowSums(M == maximum) > 1L
  )
}

make_humidity_X <- function(z, week) {
  phase <- 2 * pi * week / 52
  cbind(humidity = as.numeric(z),
        humidity_sin = as.numeric(z) * sin(phase),
        humidity_cos = as.numeric(z) * cos(phase))
}

# Check saved inputs against the supplied data, including the two name fixes.
pop_data <- read.csv(file.path(project_root, "data", "japan_pop.csv"),
                     encoding = "UTF-8") %>%
  mutate(prefecture = sub("-.*", "", prefecture),
         pop = as.numeric(gsub(",", "", pop)))
humidity_data <- read.csv(
  file.path(project_root, "data_processed", "data_processed_humidity.csv"),
  encoding = "UTF-8"
) %>%
  transmute(prefecture = case_when(key == "Hy\u014dgo" ~ "Hyogo",
                                   key == "Naoasaki" ~ "Nagasaki", TRUE ~ key),
            year, week, humidity) %>%
  filter(year %in% c(train_years, test_years))
stopifnot(
  setequal(unique(humidity_data$prefecture), pop_data$prefecture),
  !anyDuplicated(humidity_data[c("prefecture", "year", "week")]),
  all(is.finite(humidity_data$humidity))
)
humidity_by_island <- humidity_data %>%
  left_join(select(pop_data, prefecture, island), by = "prefecture") %>%
  group_by(island, year, week) %>%
  summarise(q_gkg = mean(humidity) * 1000, n_prefectures = n(), .groups = "drop")

result_dirs <- file.path(project_root, "analysis", island_folders)
names(result_dirs) <- names(island_folders)
ready <- file.exists(file.path(result_dirs, "humidity_seasonal_fit_cache.rds")) &
  file.exists(file.path(result_dirs, "humidity_seasonal_result.rds"))
if (!all(ready)) {
  stop("Missing saved fits/results: ", paste(names(result_dirs)[!ready], collapse = ", "))
}

all_results <- lapply(names(result_dirs), function(id) tryCatch({
  message("Comparing A and B: ", id)
  fit <- readRDS(file.path(result_dirs[[id]], "humidity_seasonal_fit_cache.rds"))
  res <- readRDS(file.path(result_dirs[[id]], "humidity_seasonal_result.rds"))
  tr <- as.data.frame(res$train)
  te <- as.data.frame(res$test)
  if (!isTRUE(fit$diagnostics$summary$diagnostics_ok) ||
      !isTRUE(res$diagnostics$summary$diagnostics_ok)) {
    stop("MCMC diagnostics failed: ", id)
  }
  stopifnot(
    fit$standata$N == 260L, fit$standata$K == 3L,
    identical(as.integer(tr$year), rep(train_years, each = 52L)),
    identical(as.integer(tr$week), rep(1:52, length(train_years))),
    identical(as.integer(te$year), rep(test_years, each = 52L)),
    identical(as.integer(te$week), rep(1:52, length(test_years))),
    isTRUE(all.equal(fit$standata, res$standata, tolerance = 1e-12)),
    isTRUE(all.equal(as.numeric(fit$standata$pop),
                     sum(pop_data$pop[pop_data$island == id]), tolerance = 1e-12)),
    identical(as.integer(fit$standata$cases), as.integer(tr$cases)),
    identical(as.integer(fit$standata$week), as.integer(tr$week))
  )

  weather <- humidity_by_island[humidity_by_island$island == id, ]
  weather <- weather[order(weather$year, weather$week), ]
  saved <- rbind(tr, te)
  stopifnot(nrow(weather) == nrow(saved),
            all(weather$n_prefectures == sum(pop_data$island == id)),
            identical(as.integer(weather$year), as.integer(saved$year)),
            identical(as.integer(weather$week), as.integer(saved$week)))
  if (!isTRUE(all.equal(as.numeric(saved$q_gkg), as.numeric(weather$q_gkg),
                        tolerance = 1e-8))) {
    stop("Saved humidity differs from the complete prefecture mean: ", id,
         ". Refit that island after correcting CSV encoding and prefecture names.")
  }

  q_reference <- tapply(tr$q_gkg, tr$week, mean)
  train_anomaly <- tr$q_gkg - as.numeric(q_reference[as.character(tr$week)])
  q_scale <- sd(train_anomaly)
  stopifnot(is.finite(q_scale), q_scale > 0)
  z_train <- train_anomaly / q_scale
  z_actual <- (te$q_gkg - as.numeric(q_reference[as.character(te$week)])) / q_scale
  X_A <- make_humidity_X(z_actual, te$week)
  X_B <- make_humidity_X(rep(0, nrow(te)), te$week)
  stopifnot(
    isTRUE(all.equal(as.numeric(tr$q_anomaly_std), as.numeric(z_train), tolerance = 1e-8)),
    isTRUE(all.equal(as.numeric(te$q_anomaly_std), as.numeric(z_actual), tolerance = 1e-8)),
    isTRUE(all.equal(unname(fit$standata$X),
                     unname(make_humidity_X(z_train, tr$week)), tolerance = 1e-8)),
    isTRUE(all.equal(unname(res$X_test), unname(X_A), tolerance = 1e-8))
  )

  # A and B share posterior draws and the fitted end-of-2017 S/I/R states.
  # Simulate all 104 weeks continuously; do not reset states at January 2019.
  C_A <- hs_predict(fit$draws, te$week, fit$standata, X_A)
  C_B <- hs_predict(fit$draws, te$week, fit$standata, X_B)
  stopifnot(identical(dim(C_A), dim(C_B)))
  if (!isTRUE(all.equal(unname(C_A), unname(res$C), tolerance = 1e-10))) {
    stop("Saved forecast and fit cache do not match: ", id,
         ". Refresh that island's results with force_refit <- FALSE.")
  }

  bind_rows(lapply(test_years, function(yr) {
    idx <- which(te$year == yr)
    A <- peak_info(C_A[, idx, drop = FALSE], te$week[idx])
    B <- peak_info(C_B[, idx, drop = FALSE], te$week[idx])
    data.frame(island = id, year = yr, draw = seq_len(nrow(A)),
               peak_week_A = A$peak_week, peak_week_B = B$peak_week,
               shift_weeks = A$peak_week - B$peak_week,
               boundary = A$boundary | B$boundary,
               tied = A$tied | B$tied)
  }))
}, error = function(e) {
  structure(list(island = id, message = conditionMessage(e)), class = "peak_input_error")
}))
failed <- vapply(all_results, inherits, logical(1), what = "peak_input_error")
if (any(failed)) {
  stop("Peak comparison stopped. Resolve the following saved-input problems:\n",
       paste(vapply(all_results[failed], function(x) {
         paste0(x$island, ": ", x$message)
       }, character(1)), collapse = "\n"), call. = FALSE)
}

peak_comparison_all <- bind_rows(all_results) %>%
  mutate(island = factor(island, levels = names(island_folders)))
peak_summary_all <- peak_comparison_all %>%
  group_by(island, year) %>%
  summarise(
    n_draws = n(),
    median_shift_weeks = median(shift_weeks),
    lower_95 = quantile(shift_weeks, 0.025),
    upper_95 = quantile(shift_weeks, 0.975),
    prob_earlier = mean(shift_weeks < 0),
    prob_same = mean(shift_weeks == 0),
    prob_later = mean(shift_weeks > 0),
    .groups = "drop"
  )
peak_checks_all <- peak_comparison_all %>%
  group_by(island, year) %>%
  summarise(n_draws = n(), boundary_draws = sum(boundary),
            tied_draws = sum(tied), outside_zoom_draws = sum(abs(shift_weeks) > 8),
            outside_zoom_percent = 100 * mean(abs(shift_weeks) > 8), .groups = "drop")
peak_distribution_all <- peak_comparison_all %>%
  count(island, year, shift_weeks, name = "n") %>%
  group_by(island, year) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(direction = factor(case_when(
    shift_weeks < 0 ~ "Earlier", shift_weeks == 0 ~ "Same week", TRUE ~ "Later"
  ), levels = c("Earlier", "Same week", "Later")))

g_peak_all <- ggplot(peak_distribution_all,
                     aes(x = shift_weeks, y = proportion, fill = direction)) +
  geom_col(width = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_grid(island ~ year) +
  scale_fill_manual(values = c(Earlier = "#0072B2", "Same week" = "#999999",
                                Later = "#D55E00"), drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Peak timing difference: Hokkaido, Shikoku and Ryukyu",
       subtitle = "A: actual humidity; B: each island's 2013-2017 weekly-mean humidity",
       x = "Peak timing difference, A - B (weeks)", y = "Posterior proportion", fill = NULL,
       caption = "Peaks of expected reported cases. All posterior draws are retained.") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(size = 13, face = "bold"),
        strip.text.y = element_text(angle = 0), legend.position = "bottom",
        panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))

# Zoom only the view; the percentages still use every posterior draw.
g_peak_all_zoom <- g_peak_all +
  scale_x_continuous(breaks = seq(-8, 8, 2)) +
  coord_cartesian(xlim = c(-8.5, 8.5)) +
  labs(caption = NULL)
g_peak_shift <- g_peak_all_zoom

print(peak_summary_all, n = Inf, width = Inf)
print(peak_checks_all, n = Inf, width = Inf)
print(g_peak_all_zoom)

dir.create(combined_output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(combined_output_dir, "peak_shift_three_islands_2018_2019.png"),
       g_peak_all, width = 12, height = 10, units = "in", dpi = 220, bg = "white")
ggsave(file.path(combined_output_dir, "peak_shift_three_islands_2018_2019_zoom.png"),
       g_peak_all_zoom, width = 12, height = 10, units = "in", dpi = 220, bg = "white")
write.csv(peak_summary_all, file.path(combined_output_dir, "peak_shift_summary_three_islands_2018_2019.csv"),
          row.names = FALSE)
write.csv(peak_checks_all, file.path(combined_output_dir, "peak_shift_checks_three_islands_2018_2019.csv"),
          row.names = FALSE)
write.csv(peak_distribution_all, file.path(combined_output_dir, "peak_shift_distribution_three_islands_2018_2019.csv"),
          row.names = FALSE)
message("Saved Hokkaido, Shikoku and Ryukyu peak comparisons in: ", combined_output_dir)
