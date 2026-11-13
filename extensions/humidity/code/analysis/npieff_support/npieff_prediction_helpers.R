# Prediction and reporting helpers. Sourcing this file only defines functions.
# Full npieff trajectories come from the joint Stan posterior, including holdout.

np_predict <- function(draws, week, standata, X_test) {
  stopifnot(is.matrix(draws), is.numeric(draws), nrow(draws) > 0L,
            !is.null(colnames(draws)), !anyDuplicated(colnames(draws)),
            is.matrix(X_test), is.numeric(X_test),
            length(standata$K) == 1L, standata$K %in% c(0L, 3L),
            length(standata$N_train) == 1L,
            is.finite(standata$N_train), standata$N_train >= 2L,
            standata$N_train == as.integer(standata$N_train),
            length(standata$N_total) == 1L,
            is.finite(standata$N_total),
            standata$N_total == as.integer(standata$N_total),
            standata$N_total > standata$N_train,
            length(standata$Nnonpi) == 1L,
            is.finite(standata$Nnonpi), standata$Nnonpi >= 0L,
            standata$Nnonpi == as.integer(standata$Nnonpi),
            standata$Nnonpi <= standata$N_total,
            length(week) == standata$N_total - standata$N_train,
            length(week) == nrow(X_test), length(week) > 0L,
            ncol(X_test) == standata$K,
            all(is.finite(X_test)), all(is.finite(week)),
            all(week %in% 1:52), all(week == as.integer(week)),
            length(standata$pop) == 1L, is.finite(standata$pop),
            standata$pop > 0,
            length(standata$mu) == 1L, is.finite(standata$mu),
            standata$mu > 0,
            length(standata$gamma) == 1L, is.finite(standata$gamma),
            standata$gamma > 0)

  beta_names <- paste0("beta[", 1:52, "]")
  # paste0 with seq_len(0) does not reliably produce character(0).
  theta_names <- if (standata$K > 0L)
    paste0("theta[", seq_len(standata$K), "]") else character()
  state_names <- paste0(c("S", "I", "R"), "[", standata$N_train, "]")
  npieff_names <- paste0("npieff[", seq_len(standata$N_total), "]")
  required <- c(beta_names, theta_names, state_names,
                "rho", "tau", "omega", "phi", npieff_names)
  missing <- setdiff(required, colnames(draws))
  if (length(missing)) {
    stop("Missing prediction columns: ", paste(missing, collapse = ", "))
  }
  stopifnot(all(is.finite(draws[, required, drop = FALSE])))

  beta <- draws[, beta_names, drop = FALSE]
  theta <- draws[, theta_names, drop = FALSE]
  npieff <- draws[, npieff_names, drop = FALSE]
  S <- draws[, state_names[1]]
  I <- draws[, state_names[2]]
  R <- draws[, state_names[3]]
  rho <- draws[, "rho"]
  tau <- draws[, "tau"]
  omega <- draws[, "omega"]
  stopifnot(all(beta > 0), all(S >= 0), all(I >= 0), all(R >= 0),
            all(rho >= 0 & rho <= 1), all(tau > 0), all(omega >= 0),
            all(draws[, "phi"] > 0), all(npieff > 0))
  if (standata$Nnonpi > 0L &&
      any(npieff[, seq_len(standata$Nnonpi), drop = FALSE] != 1)) {
    stop("The pre-onset npieff path must equal one through Nnonpi.")
  }

  delta <- 1 / tau
  mu <- standata$mu
  gamma <- standata$gamma
  pop <- standata$pop
  C <- matrix(NA_real_, nrow(draws), length(week))

  for (t in seq_along(week)) {
    log_multiplier <- if (standata$K > 0L)
      as.vector(theta %*% X_test[t, ]) else rep(0, nrow(draws))
    beta_t <- beta[, week[t]] * exp(log_multiplier) *
      npieff[, standata$N_train + t]
    foi <- beta_t * (I + omega) / pop
    if (any(!is.finite(beta_t)) || any(beta_t < 0) ||
        any(!is.finite(foi)) || any(foi < 0)) {
      stop("Non-finite/negative transmission or force of infection at holdout step ",
           t, ". No invalid draws were removed or clamped.")
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


np_evaluate <- function(C, draws, test_data, id, seed = 202L) {
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

np_diagnostics <- function(fit, Nnonpi, max_treedepth = 15L) {
  stopifnot(length(Nnonpi) == 1L, is.finite(Nnonpi), Nnonpi >= 0L,
            Nnonpi == as.integer(Nnonpi))
  if (!inherits(fit, "stanfit")) stop("Expected a stanfit object.")
  ss <- rstan::summary(fit)$summary
  stopifnot(all(c("sd", "n_eff", "Rhat") %in% colnames(ss)))
  parameters <- data.frame(
    parameter = rownames(ss), sd = ss[, "sd"],
    n_eff = ss[, "n_eff"], Rhat = ss[, "Rhat"],
    row.names = NULL
  )
  # Only mathematically constant outputs are omitted. Stuck parameters fail.
  fixed_parameters <- c("C[1]", if (Nnonpi > 0L)
    paste0("npieff[", seq_len(Nnonpi), "]") else character())
  checked <- parameters[!parameters$parameter %in% fixed_parameters, , drop = FALSE]
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

np_plot_comparison <- function(curves) {
  required <- c("island", "model", "step", "observed", "lower",
                "pred_median", "upper")
  island_levels <- c("Hokkaido", "Shikoku", "Ryukyu")
  model_levels <- c("Baseline", "Humidity anomaly (seasonal + npieff)")
  stopifnot(is.data.frame(curves), nrow(curves) > 0L,
            all(required %in% names(curves)),
            all(as.character(curves$island) %in% island_levels),
            all(as.character(curves$model) %in% model_levels),
            all(is.finite(as.matrix(curves[, c("step", "observed", "lower",
                                               "pred_median", "upper")]))),
            all(curves$lower <= curves$pred_median),
            all(curves$pred_median <= curves$upper))
  curves$island <- factor(curves$island, levels = island_levels)
  curves$model <- factor(curves$model, levels = model_levels)
  curves <- curves[order(curves$island, curves$model, curves$step), ]

  # Export at width = 14, height = 10 inches in the calling workflow.
  # facet_grid free_y shares one vertical scale across both models per island.
  ggplot2::ggplot(curves, ggplot2::aes(x = step)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#90CAF9", alpha = 0.45
    ) +
    ggplot2::geom_line(ggplot2::aes(y = observed), colour = "black",
                       linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = pred_median), colour = "#0072B2",
                       linewidth = 0.55) +
    ggplot2::geom_vline(xintercept = 52.5, linetype = "dotted") +
    ggplot2::facet_grid(island ~ model, scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 27, 53, 79, 104),
      labels = c("2018 W1", "W27", "2019 W1", "W27", "W52")
    ) +
    ggplot2::labs(
      title = "2018-2019 forecasts: npieff in both models",
      subtitle = "Training: 2013-2017; continuous 104-week forecast",
      x = NULL, y = "Weekly reported cases",
      caption = paste("Black: observed; blue: predictive median;",
                      "shading: 95% predictive interval.")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8))
}
