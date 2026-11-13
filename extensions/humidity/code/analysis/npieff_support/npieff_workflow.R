# Shared workflow. Source npieff_prediction_helpers.R before calling these functions.
# Data: Park et al. CSVs only. NPI prior horizon is explicitly 2013-2019.
np_same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12))

np_read_data <- function(project_root) {
  p <- read.csv(file.path(project_root, "data", "japan_pop.csv"), encoding = "UTF-8")
  p$prefecture <- sub("-.*", "", p$prefecture)
  p$pop <- as.numeric(gsub(",", "", p$pop, fixed = TRUE))
  stopifnot(!anyDuplicated(p$prefecture), all(is.finite(p$pop)), all(p$pop > 0))
  r <- read.csv(file.path(project_root, "data", "japan_rsv.csv"), encoding = "UTF-8")
  r <- r[r$year %in% 2013:2019, ]
  r$week[r$week == 53] <- 52L
  stopifnot(all(is.finite(r$cases)), all(r$cases >= 0), all(r$prefecture %in% p$prefecture))
  r <- aggregate(cases ~ prefecture + year + week, r, function(x) round(mean(x)))
  r$island <- p$island[match(r$prefecture, p$prefecture)]
  r <- aggregate(cases ~ island + year + week, r, sum)
  q <- read.csv(file.path(project_root, "data_processed", "data_processed_humidity.csv"),
                 encoding = "UTF-8")
  q$key[q$key == "Hy\u014dgo"] <- "Hyogo"
  q$key[q$key == "Naoasaki"] <- "Nagasaki"
  q <- q[q$year %in% 2013:2019, ]
  list(population = p, cases = r, humidity = q)
}

np_make_inputs <- function(raw, island, model) {
  onset <- c(Hokkaido = 208L, Shikoku = 364L, Ryukyu = 312L)
  n_pref <- c(Hokkaido = 1L, Shikoku = 4L, Ryukyu = 1L)
  stopifnot(island %in% names(onset), model %in% c("baseline", "humidity_seasonal"))
  p <- raw$population[raw$population$island == island, ]
  stopifnot(nrow(p) == n_pref[[island]])
  q <- raw$humidity[raw$humidity$key %in% p$prefecture, ]
  stopifnot(nrow(q) == 364L * nrow(p),
            setequal(unique(q$key), p$prefecture),
            !anyNA(q[c("key", "year", "week", "humidity")]),
            !anyDuplicated(q[c("key", "year", "week")]), all(is.finite(q$humidity)))
  counts <- aggregate(humidity ~ year + week, q, length)
  stopifnot(nrow(counts) == 364L, all(counts$humidity == nrow(p)))
  q <- aggregate(humidity ~ year + week, q, mean)
  q$q_gkg <- q$humidity * 1000
  cases <- raw$cases[raw$cases$island == island, c("year", "week", "cases")]
  d <- merge(cases, q[c("year", "week", "q_gkg")], by = c("year", "week"), all = TRUE)
  d <- d[order(d$year, d$week), ]
  rownames(d) <- NULL
  stopifnot(nrow(d) == 364L, all(is.finite(as.matrix(d))),
            identical(as.integer(d$year), rep(2013:2019, each = 52L)),
            identical(as.integer(d$week), rep(1:52, 7L)),
            all(d$cases >= 0), all(d$cases == floor(d$cases)))
  tr <- d[d$year <= 2017, ]
  te <- d[d$year >= 2018, ]
  q_reference <- tapply(tr$q_gkg, tr$week, mean)
  q_scale <- sd(tr$q_gkg - q_reference[as.character(tr$week)])
  stopifnot(is.finite(q_scale), q_scale > 0)
  tr$q_anomaly_std <- as.numeric((tr$q_gkg - q_reference[as.character(tr$week)]) / q_scale)
  te$q_anomaly_std <- as.numeric((te$q_gkg - q_reference[as.character(te$week)]) / q_scale)
  K <- if (model == "baseline") 0L else 3L
  feature <- function(x) {
    if (K == 0L) return(matrix(numeric(), nrow(x), 0L))
    phase <- 2 * pi * x$week / 52
    z <- x$q_anomaly_std
    cbind(humidity = z, humidity_sin = z * sin(phase), humidity_cos = z * cos(phase))
  }
  # The Stan case array and covariate matrix contain TRAINING ROWS ONLY.
  stan <- list(N_train = 260L, N_total = 364L, Nnonpi = onset[[island]],
               npi_anchor = if (island == "Shikoku") 0L else 1L, K = K,
               week = as.integer(tr$week), cases = as.integer(tr$cases),
               X_train = feature(tr), theta_sd = if (K) c(0.3, 0.15, 0.15) else numeric(),
               mu = 1 / 80 / 52, pop = sum(p$pop), gamma = 1)
  list(island = island, model = model, standata = stan, train = tr, test = te,
       X_test = feature(te), q_reference = q_reference, q_scale = q_scale,
       prefectures = as.character(p$prefecture), train_years = 2013:2017,
       test_years = 2018:2019, npi_prior_end_year = 2019L)
}

np_initial_values <- function(chain, data) {
  M <- data$N_total - data$Nnonpi
  list(S0 = 0.07 + 0.005 * chain, I0 = 0.00008 + 0.00002 * chain,
       rho = 0.008 + 0.001 * chain, beta = rep(8, 52), sigma = 0.4,
       phi = 8 + chain, tau = 100 + 10 * chain, omega = 80 + 20 * chain,
       sigma_npi = 0.04 + 0.005 * chain,
       npi_first = if (M > 0L) array(0.95 + 0.025 * chain, dim = 1L) else array(numeric(), dim = 0L),
       z_npi = rep(0, max(M - 1L, 0L)), theta = rep(0, data$K))
}

np_save <- function(object, path) {
  if (file.exists(path)) {
    backup <- paste0(path, ".before_write_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    k <- 0L
    while (file.exists(backup)) {
      k <- k + 1L
      backup <- paste0(path, ".before_write_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", k)
    }
    if (!file.copy(path, backup, overwrite = FALSE)) stop("Cannot back up ", path)
  }
  saveRDS(object, path)
  invisible(path)
}

np_fit_one <- function(inputs, stan_file, model_md5, settings, output_root, force_refit, compiled) {
  folder <- file.path(output_root, tolower(inputs$island), inputs$model)
  cache_file <- file.path(folder, "fit_cache.rds")
  result_file <- file.path(folder, "result.rds")
  attempt_file <- file.path(folder, "latest_attempt.rds")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  attempt <- list(island = inputs$island, model = inputs$model,
                  model_md5 = model_md5, sampling_settings = settings,
                  started_at = Sys.time(), status = "running")
  saveRDS(attempt, attempt_file)
  finished <- FALSE
  on.exit({
    if (!finished) {
      attempt$status <- "interrupted_or_error"
      saveRDS(attempt, attempt_file)
    }
  }, add = TRUE)
  cache <- NULL
  if (!force_refit && file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    if (!identical(cache$schema, "npieff_paired_v1") ||
        !identical(cache$island, inputs$island) || !identical(cache$model, inputs$model) ||
        !identical(cache$model_md5, model_md5) || !np_same(cache$standata, inputs$standata) ||
        !np_same(cache$sampling_settings, settings)) {
      stop("Saved fit differs from this model/input/settings: ", inputs$island, "/", inputs$model,
           ". Review before setting force_refit <- TRUE.")
    }
    message("Reusing ", inputs$island, " / ", inputs$model, "; no MCMC.")
  }
  if (is.null(cache)) {
    if (!exists("model", envir = compiled, inherits = FALSE)) {
      rstan::rstan_options(auto_write = FALSE)
      compiled$model <- rstan::stan_model(file = stan_file, model_name = "sirs_npieff_holdout")
    }
    inits <- lapply(seq_len(settings$chains), np_initial_values, data = inputs$standata)
    message("Fitting ", inputs$island, " / ", inputs$model, " with npieff; training cases end in 2017.")
    started <- proc.time()[["elapsed"]]
    fit <- rstan::sampling(compiled$model, data = inputs$standata, init = inits,
                           seed = settings$seed, chains = settings$chains, cores = settings$cores,
                           iter = settings$iter, warmup = settings$warmup, refresh = 200,
                           control = list(adapt_delta = settings$adapt_delta,
                                          max_treedepth = settings$max_treedepth))
    pars <- c("beta", if (inputs$standata$K > 0L) "theta", "rho", "tau", "omega", "phi", "sigma",
              "S0", "I0", "sigma_npi", "npieff",
              paste0(c("S", "I", "R"), "[", inputs$standata$N_train, "]"), "lp__")
    a <- as.array(fit, pars = pars)
    stopifnot(dim(a)[1] == settings$iter - settings$warmup, dim(a)[2] == settings$chains)
    cache <- list(schema = "npieff_paired_v1", island = inputs$island, model = inputs$model,
                  model_md5 = model_md5, standata = inputs$standata, sampling_settings = settings,
                  draw_array = a, draws = matrix(a, ncol = dim(a)[3],
                                                 dimnames = list(NULL, dimnames(a)[[3]])),
                  diagnostics = np_diagnostics(fit, inputs$standata$Nnonpi, settings$max_treedepth),
                  fitted_at = Sys.time(), elapsed_seconds = proc.time()[["elapsed"]] - started,
                  npi_prior_end_year = 2019L)
    np_save(cache, cache_file)
  }
  print(cache$diagnostics$summary, row.names = FALSE, digits = 5)
  res <- c(inputs, list(schema = "npieff_paired_v1", model_md5 = model_md5,
                        fit_cache_md5 = unname(tools::md5sum(cache_file)),
                        diagnostics = cache$diagnostics, C = NULL, result = NULL,
                        status = "diagnostics_failed"))
  if (!isTRUE(cache$diagnostics$summary$diagnostics_ok)) {
    np_save(res, result_file)
    attempt$status <- "diagnostics_failed"
    saveRDS(attempt, attempt_file)
    finished <- TRUE
    message(inputs$island, " / ", inputs$model,
            ": saved failed diagnostics; no prediction generated. No automatic retry.")
    return(data.frame(island = inputs$island, model = inputs$model, status = "diagnostics_failed"))
  }
  res$C <- np_predict(cache$draws, inputs$test$week, inputs$standata, inputs$X_test)
  res$result <- np_evaluate(res$C, cache$draws, inputs$test, id = inputs$model, seed = 202L)
  res$status <- "ok"
  np_save(res, result_file)
  attempt$status <- "ok"
  attempt$fit_cache_md5 <- unname(tools::md5sum(cache_file))
  attempt$result_md5 <- unname(tools::md5sum(result_file))
  saveRDS(attempt, attempt_file)
  finished <- TRUE
  message("Saved predictions: ", result_file)
  data.frame(island = inputs$island, model = inputs$model, status = "ok")
}

np_run <- function(project_root, stan_file, output_root, islands, models,
                   settings, validation_only = FALSE, force_refit = FALSE) {
  stopifnot(all(islands %in% c("Hokkaido", "Shikoku", "Ryukyu")),
            all(models %in% c("baseline", "humidity_seasonal")),
            settings$chains == 4L, settings$cores >= 1L,
            settings$iter > settings$warmup, settings$warmup > 0L,
            settings$adapt_delta > 0, settings$adapt_delta < 1,
            settings$max_treedepth >= 1L, file.exists(stan_file))
  raw <- np_read_data(project_root)
  jobs <- unlist(lapply(islands, function(island) lapply(models, function(model)
    np_make_inputs(raw, island, model))), recursive = FALSE)
  input_summary <- do.call(rbind, lapply(jobs, function(x) data.frame(
    island = x$island, model = x$model, training_weeks = x$standata$N_train,
    forecast_weeks = x$standata$N_total - x$standata$N_train,
    active_npi_training = max(x$standata$N_train - x$standata$Nnonpi, 0L),
    active_npi_forecast = x$standata$N_total - max(x$standata$Nnonpi, x$standata$N_train))))
  print(input_summary, row.names = FALSE)
  if (validation_only) {
    message("INPUT VALIDATION ONLY: no compilation, MCMC, predictions, or output writes.")
    return(invisible(list(inputs = jobs, input_summary = input_summary)))
  }
  Sys.setenv(PERTURBATION_ROOT = project_root)
  source(file.path(project_root, "script", "setup.R"), encoding = "UTF-8")
  perturbation_setup("analysis")
  if (!requireNamespace("rstan", quietly = TRUE)) stop("Project rstan is unavailable.")
  model_md5 <- unname(tools::md5sum(stan_file))
  compiled <- new.env(parent = emptyenv())
  rows <- lapply(jobs, function(x) {
    s <- settings
    s$seed <- as.integer(settings$seed + 100L * (match(x$island, c("Hokkaido", "Shikoku", "Ryukyu")) - 1L) +
                           (match(x$model, c("baseline", "humidity_seasonal")) - 1L))
    tryCatch(np_fit_one(x, stan_file, model_md5, s, output_root, force_refit, compiled),
             error = function(e) {
               message(x$island, " / ", x$model, ": ", conditionMessage(e))
               data.frame(island = x$island, model = x$model, status = "error")
             })
  })
  status <- do.call(rbind, rows)
  print(status, row.names = FALSE)
  if (any(status$status != "ok")) warning("Some fits are not ready. Inspect the printed diagnostics/errors; no automatic refits.", call. = FALSE)
  invisible(list(status = status, input_summary = input_summary))
}

np_collect <- function(project_root, stan_file, output_root, islands, models) {
  raw <- np_read_data(project_root)
  current_md5 <- unname(tools::md5sum(stan_file))
  curves <- metrics <- checks <- list()
  for (island in islands) for (model in models) {
    folder <- file.path(output_root, tolower(island), model)
    cf <- file.path(folder, "fit_cache.rds")
    rf <- file.path(folder, "result.rds")
    af <- file.path(folder, "latest_attempt.rds")
    if (!all(file.exists(c(cf, rf, af)))) stop("Missing npieff fit/predictions: ", island, "/", model,
      ". First run fit_baseline_and_humidity_npieff_2013_2017.R. Old no-npieff results are not used.", call. = FALSE)
    attempt <- readRDS(af)
    if (!identical(attempt$status, "ok"))
      stop("Latest attempt is not ready: ", island, "/", model, " (", attempt$status,
           "). Earlier successful results will not be plotted.", call. = FALSE)
    fit <- readRDS(cf); res <- readRDS(rf)
    input <- np_make_inputs(raw, island, model)
    stopifnot(identical(fit$schema, "npieff_paired_v1"), identical(res$schema, "npieff_paired_v1"),
              identical(fit$island, island), identical(res$island, island),
              identical(fit$model, model), identical(res$model, model),
              identical(fit$model_md5, current_md5), identical(res$model_md5, current_md5),
              identical(attempt$model_md5, current_md5),
              identical(attempt$fit_cache_md5, unname(tools::md5sum(cf))),
              identical(attempt$result_md5, unname(tools::md5sum(rf))),
              np_same(attempt$sampling_settings, fit$sampling_settings),
              identical(res$fit_cache_md5, unname(tools::md5sum(cf))),
              np_same(fit$standata, res$standata), np_same(res$standata, input$standata),
              np_same(fit$diagnostics, res$diagnostics), np_same(res$X_test, input$X_test),
              np_same(res$test, input$test), np_same(res$train, input$train))
    if (!identical(res$status, "ok") || !isTRUE(fit$diagnostics$summary$diagnostics_ok))
      stop("MCMC diagnostics failed: ", island, "/", model, ". No comparison plotted.", call. = FALSE)
    curve <- as.data.frame(res$result$curve)
    stopifnot(identical(as.integer(curve$year), rep(2018:2019, each = 52)),
              identical(as.integer(curve$week), rep(1:52, 2)),
              identical(as.integer(curve$step), 1:104),
              np_same(as.numeric(curve$observed), as.numeric(input$test$cases)),
              all(is.finite(as.matrix(curve[c("observed", "pred_median", "lower", "upper")]))),
              all(curve$lower <= curve$pred_median & curve$pred_median <= curve$upper))
    label <- if (model == "baseline") "Baseline" else "Humidity anomaly (seasonal + npieff)"
    key <- paste(island, model)
    curve$island <- island; curve$model <- label
    curves[[key]] <- curve
    metrics[[key]] <- cbind(island = island, model = label,
                            res$result$metrics[c("mean_log_score", "MAE", "coverage_95")])
    checks[[key]] <- cbind(island = island, model = label, fit$diagnostics$summary)
  }
  list(curves = do.call(rbind, curves), metrics = do.call(rbind, metrics),
       diagnostics = do.call(rbind, checks))
}
