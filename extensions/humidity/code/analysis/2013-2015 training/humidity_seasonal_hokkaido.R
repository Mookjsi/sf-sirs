# Hokkaido: humidity_seasonal only.
# This single R file contains its Stan model and all analysis functions.
# In RStudio, open this file and click Source.
# Defaults: reuse matching saved posterior draws; fit ONLY this model if absent.
# The original tuning scripts and results are not modified.
#
# Model:
# beta_t = beta[week] * exp(z * (theta[1] +
#                theta[2]*sin(2*pi*week/52) + theta[3]*cos(2*pi*week/52)))
# z = (specific humidity - 2013-2015 same-week mean) / training anomaly SD.
#
# Change these settings here:
project_root <- Sys.getenv("PERTURBATION_ROOT")
island_id <- "Hokkaido"

output_dir <- file.path(
  project_root,
  "analysis",
  paste0("humidity_seasonal_", tolower(island_id), "_train2013_2015")
)

force_refit <- FALSE

sampling_settings <- list(
  chains = 4L, cores = 4L, iter = 12000L, warmup = 8000L, seed = 303L,
  adapt_delta = 0.95, max_treedepth = 15L
)

# Embedded Stan model. The seasonal design below fixes K=3.
stan_code <- 'data {
  int<lower=2> N;
  array[N] int<lower=1, upper=52> week;
  array[N] int<lower=0> cases;
  real<lower=0> mu;
  real<lower=0> pop;
  real<lower=0> gamma;

  int<lower=1> K;
  matrix[N, K] X;
  vector<lower=0>[K] theta_sd;
  real<lower=0> sigma_sd;
  real<lower=0> beta_sd;
  real<lower=0> phi_sd;
  real<lower=0> tau_sd;
  real<lower=0> omega_sd;
  real<lower=0> rho_sd;
  real<lower=0, upper=1> S0_mean;
  real<lower=0> S0_sd;
}

parameters {
  real<lower=0, upper=1> S0;
  real<lower=0, upper=1> I0;
  real<lower=0, upper=1> rho;
  array[52] real<lower=0> beta;
  real<lower=0> sigma;
  real<lower=0> phi;
  real<lower=0> tau;
  real<lower=0> omega;
  vector[K] theta;
}

transformed parameters {
  vector<lower=0>[N] S;
  vector<lower=0>[N] I;
  vector<lower=0>[N] R;
  vector<lower=0>[N] C;
  real<lower=0> delta = 1 / tau;

  S[1] = S0 * pop;
  I[1] = I0 * pop;
  R[1] = (1 - S0 - I0) * pop;
  C[1] = 0;

  for (i in 2:N) {
    real foi = beta[week[i]] * exp(X[i] * theta)
               * (I[i - 1] + omega) / pop;
    real Sout = (1 - exp(-(foi + mu))) * S[i - 1];
    real StoI = foi / (foi + mu) * Sout;
    real Iout = (1 - exp(-(gamma + mu))) * I[i - 1];
    real ItoR = gamma / (gamma + mu) * Iout;
    real Rout = (1 - exp(-(delta + mu))) * R[i - 1];
    real RtoS = delta / (delta + mu) * Rout;

    S[i] = S[i - 1] + RtoS - Sout + mu * pop;
    I[i] = I[i - 1] + StoI - Iout;
    R[i] = R[i - 1] + ItoR - Rout;
    C[i] = StoI * rho;
  }
}

model {
  theta ~ normal(0, theta_sd);
  beta ~ normal(8, beta_sd);
  for (i in 2:52) {
    beta[i] ~ normal(beta[i - 1], sigma);
  }
  beta[1] ~ normal(beta[52], sigma);
  sigma ~ normal(0, sigma_sd);
  tau ~ normal(0, tau_sd);
  S0 ~ normal(S0_mean, S0_sd);
  I0 ~ normal(0, 0.001);
  rho ~ normal(0, rho_sd);
  phi ~ normal(0, phi_sd);
  omega ~ normal(0, omega_sd);
  for (i in 2:N) {
    cases[i] ~ neg_binomial_2(C[i], phi);
  }
}'

hs_saved_inputs <- function(project_root) {
  path <- file.path(project_root, "analysis", "honshu_humidity_checkpoint.rda")
  if (!file.exists(path)) stop("Saved analysis checkpoint not found: ", path)
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  needed <- c("honshu_train", "honshu_test", "standata_baseline",
              "standata_humidity", "draws_baseline", "draws_humidity")
  missing <- setdiff(needed, ls(e))
  if (length(missing)) stop("Missing checkpoint objects: ", paste(missing, collapse = ", "))
  out <- mget(needed, envir = e, inherits = FALSE)
  tr <- as.data.frame(out$honshu_train)
  te <- as.data.frame(out$honshu_test)
  stopifnot(nrow(tr) == 260L, nrow(te) == 104L,
            identical(as.integer(tr$year), rep(2013:2017, each = 52)),
            identical(as.integer(te$year), rep(2018:2019, each = 52)),
            identical(as.integer(tr$week), rep(1:52, 5)),
            identical(as.integer(te$week), rep(1:52, 2)))
  required_columns <- c("year", "week", "cases", "q_gkg", "q_anomaly_std")
  stopifnot(all(required_columns %in% names(tr)),
            all(required_columns %in% names(te)))
  q_ref <- tapply(tr$q_gkg, tr$week, mean)
  anomaly <- tr$q_gkg - q_ref[as.character(tr$week)]
  q_scale <- sd(anomaly)
  stopifnot(is.finite(q_scale), q_scale > 0)
  # tapply() returns a 1-D array. unname() retains its dimension attribute;
  # compare plain numeric values on both sides rather than array metadata.
  expected_train <- as.numeric(anomaly / q_scale)
  expected_test <- as.numeric((te$q_gkg - q_ref[as.character(te$week)]) / q_scale)
  if (!isTRUE(all.equal(as.numeric(tr$q_anomaly_std), expected_train, tolerance = 1e-8)) ||
      !isTRUE(all.equal(as.numeric(te$q_anomaly_std), expected_test, tolerance = 1e-8))) {
    stop("Stored humidity anomalies differ from the 2013-2017 reference. ",
         "Use the original saved checkpoint, not a modified counterfactual input.")
  }
  stopifnot(identical(as.integer(tr$cases), as.integer(out$standata_humidity$cases)),
            identical(as.integer(tr$week), as.integer(out$standata_humidity$week)),
            out$standata_humidity$N == 260L,
            out$standata_baseline$N == 260L)
  out$honshu_train <- tr[, required_columns]
  out$honshu_test <- te[, required_columns]
  out$q_reference <- q_ref
  out$q_scale <- q_scale
  out
}

hs_initial_values <- function(chain_id, K) {
  list(
    S0 = 0.07 + 0.005 * chain_id,
    I0 = 0.00008 + 0.00002 * chain_id,
    rho = 0.008 + 0.001 * chain_id,
    beta = rep(8, 52), sigma = 0.4, phi = 8 + chain_id,
    tau = 100 + 10 * chain_id, omega = 80 + 20 * chain_id,
    theta = if (K == 1L) array(0, dim = 1L) else rep(0, K)
  )
}


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

hs_evaluate <- function(C, draws, test_data, id, seed = 202L) {
  stopifnot(is.matrix(C), is.matrix(draws), nrow(C) > 0L,
            nrow(C) == nrow(draws), ncol(C) == nrow(test_data),
            "phi" %in% colnames(draws),
            all(c("year", "week", "cases") %in% names(test_data)),
            all(is.finite(C)), all(C >= 0),
            all(is.finite(test_data$cases)), all(test_data$cases >= 0),
            all(test_data$cases == floor(test_data$cases)),
            all(test_data$week %in% 1:52),
            all(diff(test_data$year * 52 + test_data$week) == 1))
  phi <- draws[, "phi"]
  stopifnot(all(is.finite(phi)), all(phi > 0))

  # Keep this reproducible without changing the caller's random-number stream.
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_rng <- if (had_rng) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_rng) {
      assign(".Random.seed", old_rng, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  # R matrices are column-major: repeat all draw-specific phi values per week.
  Y <- matrix(
    rnbinom(length(C), mu = as.vector(C),
            size = rep(phi, times = ncol(C))),
    nrow = nrow(C), ncol = ncol(C)
  )
  if (any(!is.finite(Y))) stop("Non-finite posterior predictive counts.")
  qs <- apply(Y, 2, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
  c_qs <- apply(C, 2, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
  log_score <- vapply(seq_len(ncol(C)), function(t) {
    lp <- dnbinom(test_data$cases[t], mu = C[, t], size = phi, log = TRUE)
    if (anyNA(lp)) stop("Invalid log predictive mass at step ", t)
    a <- max(lp)
    if (a == -Inf) return(-Inf)
    a + log(mean(exp(lp - a)))
  }, numeric(1))

  curve <- data.frame(
    id = id, step = seq_len(ncol(C)), year = test_data$year,
    week = test_data$week, observed = test_data$cases,
    lower = qs[1, ], pred_median = qs[2, ], upper = qs[3, ],
    pred_mean = colMeans(C), expected_median = c_qs[2, ],
    expected_lower = c_qs[1, ], expected_upper = c_qs[3, ],
    log_score = log_score
  )
  metrics <- data.frame(
    id = id,
    mean_log_score = mean(log_score),
    MAE = mean(abs(curve$observed - curve$pred_median)),
    coverage_95 = mean(curve$observed >= curve$lower &
                       curve$observed <= curve$upper),
    mean_interval_width = mean(curve$upper - curve$lower)
  )

  years <- sort(unique(test_data$year))
  peaks <- do.call(rbind, lapply(years, function(yr) {
    idx <- which(test_data$year == yr)
    obs <- curve$observed[idx]
    pred <- curve$pred_median[idx]
    wk <- curve$week[idx]
    io <- which.max(obs)
    ip <- which.max(pred)
    obs_peak <- obs[io]
    pred_peak <- pred[ip]
    data.frame(
      id = id, year = yr,
      observed_peak = obs_peak, predicted_curve_peak = pred_peak,
      observed_peak_week = wk[io], predicted_curve_peak_week = wk[ip],
      peak_height_error = pred_peak - obs_peak,
      peak_height_abs_error = abs(pred_peak - obs_peak),
      peak_height_relative_abs_error = if (obs_peak > 0)
        abs(pred_peak - obs_peak) / obs_peak else NA_real_,
      peak_week_error = wk[ip] - wk[io],
      peak_week_abs_error = abs(wk[ip] - wk[io]),
      prediction_at_observed_peak = pred[io],
      observed_peak_ties = sum(obs == obs_peak),
      predicted_peak_ties = sum(pred == pred_peak)
    )
  }))
  metrics$mean_peak_relative_abs_error <- mean(
    peaks$peak_height_relative_abs_error
  )
  metrics$mean_peak_week_abs_error <- mean(peaks$peak_week_abs_error)
  metrics$mean_peak_height_abs_error <- mean(peaks$peak_height_abs_error)

  annual_max <- do.call(rbind, lapply(years, function(yr) {
    idx <- which(test_data$year == yr)
    data.frame(id = id, year = yr, draw = seq_len(nrow(C)),
               expected_max = apply(C[, idx, drop = FALSE], 1, max),
               reported_max = apply(Y[, idx, drop = FALSE], 1, max))
  }))
  annual_max_summary <- do.call(rbind, lapply(years, function(yr) {
    rows <- annual_max[annual_max$year == yr, ]
    obs <- max(test_data$cases[test_data$year == yr])
    ec <- quantile(rows$expected_max, c(0.025, 0.5, 0.975), names = FALSE)
    py <- quantile(rows$reported_max, c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(
      id = id, year = yr, observed_max = obs,
      expected_lower = ec[1], expected_median = ec[2], expected_upper = ec[3],
      predictive_lower = py[1], predictive_median = py[2],
      predictive_upper = py[3],
      prob_predictive_ge_observed = mean(rows$reported_max >= obs)
    )
  }))

  list(curve = curve, metrics = metrics, peaks = peaks,
       annual_max = annual_max, annual_max_summary = annual_max_summary)
}

hs_diagnostics <- function(fit, max_treedepth = 15L) {
  if (!inherits(fit, "stanfit")) stop("Expected a stanfit object.")
  ss <- summary(fit)$summary
  stopifnot(all(c("sd", "n_eff", "Rhat") %in% colnames(ss)))
  parameters <- data.frame(
    parameter = rownames(ss), sd = ss[, "sd"],
    n_eff = ss[, "n_eff"], Rhat = ss[, "Rhat"],
    row.names = NULL
  )
  # C[1] is the sole fixed entry. Do not silently omit stuck parameters.
  checked <- parameters[parameters$parameter != "C[1]", , drop = FALSE]
  complete <- nrow(checked) > 0L && all(is.finite(checked$sd)) &&
    all(checked$sd > 0) && all(is.finite(checked$Rhat)) &&
    all(is.finite(checked$n_eff))
  max_rhat <- if (complete) max(checked$Rhat) else Inf
  min_neff <- if (complete) min(checked$n_eff) else 0

  samplers <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  if (!length(samplers)) stop("No post-warmup sampler diagnostics found.")
  chains <- do.call(rbind, lapply(seq_along(samplers), function(k) {
    s <- samplers[[k]]
    stopifnot(is.matrix(s), nrow(s) > 1L,
              all(c("divergent__", "treedepth__", "energy__") %in%
                    colnames(s)))
    energy <- s[, "energy__"]
    ev <- var(energy)
    bfmi <- if (is.finite(ev) && ev > 0 && all(is.finite(energy)))
      mean(diff(energy)^2) / ev else NA_real_
    data.frame(
      chain = k, retained_draws = nrow(s),
      divergences = sum(s[, "divergent__"]),
      treedepth_hits = sum(s[, "treedepth__"] >= max_treedepth),
      E_BFMI = bfmi
    )
  }))
  div <- sum(chains$divergences)
  tree <- sum(chains$treedepth_hits)
  min_bfmi <- if (all(is.finite(chains$E_BFMI))) min(chains$E_BFMI) else NA_real_
  ok <- isTRUE(complete && max_rhat <= 1.01 && min_neff >= 400 &&
                 div == 0 && tree == 0 && is.finite(min_bfmi) &&
                 min_bfmi >= 0.3)
  diagnostic_summary <- data.frame(
    max_Rhat = max_rhat, min_n_eff = min_neff,
    divergences = div, treedepth_hits = tree, min_E_BFMI = min_bfmi,
    diagnostics_ok = ok
  )
  list(summary = diagnostic_summary, chains = chains, parameters = parameters)
}

hs_plot <- function(curve, title = "Honshu: humidity tuning") {
  ggplot2::ggplot(curve, ggplot2::aes(x = step)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#90CAF9", alpha = 0.45
    ) +
    ggplot2::geom_line(ggplot2::aes(y = observed), colour = "black",
                       linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = pred_median), colour = "#0072B2",
                       linewidth = 0.55) +
    ggplot2::geom_vline(xintercept = c(52.5, 104.5, 156.5), linetype = "dotted") +
    ggplot2::facet_wrap(~id, ncol = 2) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 53, 105, 157, 208),
      labels = c("2016 W1", "2017 W1", "2018 W1", "2019 W1", "2019 W52")
    ) +
    ggplot2::labs(
      title = title, x = NULL, y = "Weekly reported cases",
      caption = paste("Black: observed; blue: predictive median;",
                      "shading: 95% predictive interval.")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8))
}


# Build the three humidity predictors. No new climate data are used.
hs_prepare <- function(inputs) {
  tr <- inputs$honshu_train
  te <- inputs$honshu_test
  feature_matrix <- function(d) {
    z <- as.numeric(d$q_anomaly_std)
    phase <- 2 * pi * d$week / 52
    cbind(humidity = z, humidity_sin = z * sin(phase),
          humidity_cos = z * cos(phase))
  }
  X_train <- feature_matrix(tr)
  X_test <- feature_matrix(te)
  base_data <- inputs$standata_humidity
  theta_sd <- setNames(c(0.3, 0.15, 0.15), colnames(X_train))
  dat <- list(
    N = nrow(tr), week = as.integer(tr$week), cases = as.integer(tr$cases),
    mu = base_data$mu, pop = base_data$pop, gamma = base_data$gamma,
    K = 3L, X = X_train, theta_sd = theta_sd,
    sigma_sd = 1, beta_sd = 2, phi_sd = 10,
    tau_sd = 200, omega_sd = 200, rho_sd = 0.02,
    S0_mean = 0.08, S0_sd = 0.02
  )
  list(standata = dat, X_test = X_test)
}

hs_code_md5 <- function(code) {
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  # Normalize CRLF and append exactly one final newline.
  code <- sub("\n+$", "", gsub("\r", "", code, fixed = TRUE))
  con <- file(f, open = "wb")
  writeBin(charToRaw(enc2utf8(paste0(code, "\n"))), con)
  close(con)
  unname(tools::md5sum(f))
}

# MAIN: the following code runs when you click Source.
Sys.setenv(PERTURBATION_ROOT = project_root)
source(file.path(project_root, "script", "setup.R"), encoding = "UTF-8")
perturbation_setup("analysis")
library(rstan)
library(ggplot2)
rstan::rstan_options(auto_write = TRUE)

library(dplyr)

source(
  file.path(project_root, "script", "script_data.R"),
  encoding = "UTF-8"
)

# Selected island: specific humidity in g/kg
island_weather <- japan_humidity_island %>%
  ungroup() %>%
  filter(island == island_id) %>%
  transmute(year, week, q_gkg = humidity * 1000)

# Reported cases and humidity, 2013-2019
island_data <- japan_rsv_island %>%
  ungroup() %>%
  filter(
    island == island_id,
    year >= 2013,
    year <= 2019
  ) %>%
  select(year, week, cases) %>%
  left_join(island_weather, by = c("year", "week")) %>%
  arrange(year, week)

# Reference: 2013-2015 same-week mean
q_reference_island <- island_data %>%
  filter(year <= 2015) %>%
  group_by(week) %>%
  summarise(q_mean = mean(q_gkg), .groups = "drop")

island_data <- island_data %>%
  left_join(q_reference_island, by = "week") %>%
  mutate(q_anomaly = q_gkg - q_mean)

q_scale_island <- sd(
  island_data$q_anomaly[island_data$year <= 2015]
)

stopifnot(is.finite(q_scale_island), q_scale_island > 0)

island_data <- island_data %>%
  mutate(q_anomaly_std = q_anomaly / q_scale_island)

island_train <- island_data %>% filter(year <= 2015)
island_test  <- island_data %>% filter(year >= 2016)

pop_island <- japan_pop_island %>%
  ungroup() %>%
  filter(island == island_id) %>%
  pull(pop)

stopifnot(
  nrow(island_train) == 156,
  nrow(island_test) == 208,
  anyDuplicated(island_data[c("year", "week")]) == 0,
  all(is.finite(island_data$cases)),
  all(is.finite(island_data$q_anomaly_std)),
  length(pop_island) == 1,
  is.finite(pop_island),
  pop_island > 0
)

# Keep these names for compatibility with the existing functions.
seasonal_inputs <- list(
  honshu_train = island_train,
  honshu_test = island_test,
  standata_humidity = list(
    mu = 1 / 80 / 52,
    pop = pop_island,
    gamma = 1
  ),
  q_reference = q_reference_island,
  q_scale = q_scale_island
)

message(
  "Island: ", island_id,
  " | Training weeks: ", nrow(island_train),
  " | Population: ", pop_island
)
seasonal_design <- hs_prepare(seasonal_inputs)
standata_humidity_seasonal <- seasonal_design$standata
honshu_train_seasonal <- seasonal_inputs$honshu_train
honshu_test_seasonal <- seasonal_inputs$honshu_test
X_test_seasonal <- seasonal_design$X_test
model_md5 <- hs_code_md5(stan_code)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
own_cache_file <- file.path(output_dir, "humidity_seasonal_fit_cache.rds")
batch_dir <- file.path(project_root, "analysis", "humidity_tuning_run1")
batch_cache_file <- file.path(batch_dir, "fit_cache", "humidity_seasonal.rds")
seasonal_cache <- NULL

if (!force_refit && file.exists(own_cache_file)) {
  candidate <- readRDS(own_cache_file)
  if (!identical(candidate$model_md5, model_md5) ||
      !isTRUE(all.equal(candidate$standata, standata_humidity_seasonal,
                        tolerance = 1e-12))) {
    stop("Saved standalone fit uses different inputs/model. Set force_refit <- TRUE to refit.")
  }
  seasonal_cache <- candidate
  message("Reusing the saved standalone posterior; no MCMC.")
}

if (FALSE) {
  candidate <- readRDS(batch_cache_file)
  manifest_file <- file.path(batch_dir, "run_manifest.rds")
  if (!file.exists(manifest_file)) stop("Original tuning manifest is missing: ", manifest_file)
  old_manifest <- readRDS(manifest_file)
  # These fingerprints identify the original model used by the saved tuning fit.
  original_raw_stan_md5 <- "7cea62caa0ecd59fc1fc2833fd38718c"
  original_normalized_stan_md5 <- "7cea62caa0ecd59fc1fc2833fd38718c"
  compatible <- identical(old_manifest$source_hashes[1], original_raw_stan_md5) &&
    identical(model_md5, original_normalized_stan_md5) &&
    identical(as.character(candidate$config$id), "humidity_seasonal") &&
    isTRUE(all.equal(candidate$standata, standata_humidity_seasonal, tolerance = 1e-12))
  if (!isTRUE(compatible)) {
    stop("The saved tuning fit does not match this seasonal model/input. ",
         "Set force_refit <- TRUE to fit this script's model.")
  }
  seasonal_cache <- candidate
  seasonal_cache$model_md5 <- model_md5
  seasonal_cache$sampling_settings <- old_manifest$settings
  message("Reusing humidity_seasonal posterior from the tuning run; no MCMC.")
}

if (is.null(seasonal_cache)) {
  model_humidity_seasonal <- rstan::stan_model(
    model_code = stan_code, model_name = "humidity_seasonal"
  )
  s <- sampling_settings
  inits <- lapply(seq_len(s$chains), hs_initial_values, K = 3L)
  started <- proc.time()[["elapsed"]]
  fit_humidity_seasonal <- rstan::sampling(
    model_humidity_seasonal,
    data = standata_humidity_seasonal,
    seed = s$seed, chains = s$chains, cores = s$cores,
    iter = s$iter, warmup = s$warmup, init = inits, refresh = 200,
    control = list(adapt_delta = s$adapt_delta,
                   max_treedepth = s$max_treedepth)
  )
  N <- standata_humidity_seasonal$N
  pars <- c("beta", "theta", "rho", "tau", "omega", "phi", "sigma", "S0", "I0",
            paste0(c("S", "I", "R"), "[", N, "]"), "lp__")
  seasonal_draw_array <- as.array(fit_humidity_seasonal, pars = pars)
  stopifnot(dim(seasonal_draw_array)[1] == s$iter - s$warmup,
            dim(seasonal_draw_array)[2] == s$chains)
  seasonal_cache <- list(
    model_md5 = model_md5, standata = standata_humidity_seasonal,
    sampling_settings = s, draw_array = seasonal_draw_array,
    draws = matrix(seasonal_draw_array, ncol = dim(seasonal_draw_array)[3],
                   dimnames = list(NULL, dimnames(seasonal_draw_array)[[3]])),
    diagnostics = hs_diagnostics(fit_humidity_seasonal, s$max_treedepth),
    fit_elapsed_seconds = proc.time()[["elapsed"]] - started
  )
  if (force_refit && file.exists(own_cache_file)) {
    backup <- paste0(own_cache_file, ".before_refit_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    if (!file.copy(own_cache_file, backup, overwrite = FALSE))
      stop("Could not back up the previous standalone fit.")
  }
}

# Save a compact local copy; the original batch cache remains unchanged.
if (force_refit || !file.exists(own_cache_file)) saveRDS(seasonal_cache, own_cache_file)
draws_humidity_seasonal <- seasonal_cache$draws
diagnostics_humidity_seasonal <- seasonal_cache$diagnostics
print(diagnostics_humidity_seasonal$summary)

if (!isTRUE(
  diagnostics_humidity_seasonal$summary$diagnostics_ok
)) {
  stop("MCMC diagnostics failed. Check chains before prediction.")
}

# Forecast continuously from the fitted end-of-2015 S/I/R states.
C_humidity_seasonal <- hs_predict(
  draws = draws_humidity_seasonal,
  week = honshu_test_seasonal$week,
  standata = standata_humidity_seasonal,
  X_test = X_test_seasonal
)

# Same NB2 predictive summaries and seed as the original tuning comparison.
result_humidity_seasonal <- hs_evaluate(
  C_humidity_seasonal, draws_humidity_seasonal,
  honshu_test_seasonal, id = "humidity_seasonal", seed = 202L
)
g_humidity_seasonal <- hs_plot(
  result_humidity_seasonal$curve,
  title = paste0(island_id, ": humidity_seasonal, 2016-2019")
)
print(diagnostics_humidity_seasonal$summary, row.names = FALSE, digits = 4)
print(result_humidity_seasonal$metrics, row.names = FALSE, digits = 4)
print(result_humidity_seasonal$peaks, row.names = FALSE, digits = 4)
print(g_humidity_seasonal)

utils::write.csv(result_humidity_seasonal$metrics,
                 file.path(output_dir, "metrics.csv"), row.names = FALSE)
utils::write.csv(result_humidity_seasonal$peaks,
                 file.path(output_dir, "peaks_2016_2019.csv"), row.names = FALSE)
utils::write.csv(result_humidity_seasonal$curve,
                 file.path(output_dir, "prediction_2016_2019.csv"), row.names = FALSE)
ggplot2::ggsave(file.path(output_dir, "humidity_seasonal_2016_2019.png"),
               g_humidity_seasonal, width = 11, height = 5, dpi = 150, bg = "white")
saveRDS(list(standata = standata_humidity_seasonal, X_test = X_test_seasonal,
             train = honshu_train_seasonal, test = honshu_test_seasonal,
             q_reference = seasonal_inputs$q_reference, q_scale = seasonal_inputs$q_scale,
             C = C_humidity_seasonal, result = result_humidity_seasonal,
             diagnostics = diagnostics_humidity_seasonal, plot = g_humidity_seasonal),
        file.path(output_dir, "humidity_seasonal_result.rds"))
message("Saved seasonal-only outputs in: ", output_dir)

