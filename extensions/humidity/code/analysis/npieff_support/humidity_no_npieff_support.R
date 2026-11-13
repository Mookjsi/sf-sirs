# Read and validate the restored seasonal-humidity model WITHOUT npieff.
# Baseline collection remains in npieff_workflow.R; baseline files are untouched.
hn_model_functions <- function(script_file) {
  expr <- parse(file = script_file, encoding = "UTF-8")
  env <- new.env(parent = globalenv())
  wanted <- c("hs_prepare", "hs_predict", "hs_evaluate", "hs_code_md5")
  for (x in expr) {
    if (!is.call(x) || !identical(x[[1]], as.name("<-")) || !is.symbol(x[[2]])) next
    nm <- as.character(x[[2]])
    if (nm == "stan_code") {
      stopifnot(is.character(x[[3]]), length(x[[3]]) == 1L)
      env$stan_code <- x[[3]]
    } else if (nm %in% wanted) {
      stopifnot(is.call(x[[3]]), identical(x[[3]][[1]], as.name("function")))
      eval(x, envir = env)
    }
  }
  stopifnot(all(vapply(c(wanted, "stan_code"), exists, logical(1),
                      envir = env, inherits = FALSE)),
            !grepl("\\b(npieff|sigma_npi|npi)\\b", env$stan_code))
  env
}

hn_collect <- function(project_root, island, script_dir) {
  filenames <- c(Hokkaido = "humidity_seasonal_hokkaido.R",
                  Shikoku = "humidity_seasonal_sikoku.R", Ryukyu = "humidity_seasonal_ryukyu.R")
  stopifnot(island %in% names(filenames))
  e <- hn_model_functions(file.path(script_dir, filenames[[island]]))
  md5 <- e$hs_code_md5(e$stan_code)
  input <- np_make_inputs(np_read_data(project_root), island, "humidity_seasonal")
  prepared <- e$hs_prepare(list(
    honshu_train = input$train, honshu_test = input$test,
    standata_humidity = input$standata[c("mu", "pop", "gamma")]
  ))
  folder <- file.path(project_root, "analysis", paste0("humidity_seasonal_", tolower(island)))
  cf <- file.path(folder, "humidity_seasonal_fit_cache.rds")
  rf <- file.path(folder, "humidity_seasonal_result.rds")
  af <- file.path(folder, "latest_attempt_no_npieff.rds")
  if (!all(file.exists(c(cf, rf)))) stop("Missing no-npieff humidity fit/predictions: ", island, call. = FALSE)
  if (file.exists(af)) {
    attempt <- readRDS(af)
    if (!identical(attempt$status, "ok"))
      stop("Latest humidity attempt is not ready: ", island, " (", attempt$status, ").", call. = FALSE)
    stopifnot(identical(attempt$model_md5, md5),
              identical(attempt$fit_cache_md5, unname(tools::md5sum(cf))),
              identical(attempt$result_md5, unname(tools::md5sum(rf))))
  }
  fit <- readRDS(cf)
  res <- readRDS(rf)
  stopifnot(identical(fit$model_md5, md5),
            np_same(fit$standata, prepared$standata),
            np_same(res$standata, prepared$standata),
            np_same(fit$diagnostics, res$diagnostics),
            np_same(res$X_test, prepared$X_test))
  if (!is.null(res$model_md5)) stopifnot(identical(res$model_md5, md5))
  if (!is.null(res$fit_cache_md5)) stopifnot(identical(res$fit_cache_md5, unname(tools::md5sum(cf))))
  for (period in c("train", "test")) for (nm in names(input[[period]])) {
    stopifnot(nm %in% names(res[[period]]),
              np_same(as.numeric(res[[period]][[nm]]), as.numeric(input[[period]][[nm]])))
  }
  ds <- fit$diagnostics$summary
  if (!isTRUE(ds$diagnostics_ok) || !isTRUE(ds$max_Rhat <= 1.01 && ds$min_n_eff >= 400 &&
      ds$divergences == 0 && ds$treedepth_hits == 0 && ds$min_E_BFMI >= 0.3)) {
    stop("MCMC diagnostics failed: ", island, "/humidity_seasonal (no npieff).", call. = FALSE)
  }
  # Legacy results have no cache hash: verify them by replaying the SAVED draws.
  # No compilation or sampling occurs. Newer result hashes are checked above too.
  C_check <- e$hs_predict(fit$draws, input$test$week, prepared$standata, prepared$X_test)
  if (!np_same(res$C, C_check)) stop("Saved humidity predictions do not match current posterior: ", island)
  predicted <- e$hs_evaluate(C_check, fit$draws, res$test, id = "humidity_seasonal", seed = 202L)
  if (!np_same(res$result$curve, predicted$curve))
    stop("Saved humidity curve differs from verified predictions: ", island)
  curve <- as.data.frame(predicted$curve)
  stopifnot(nrow(curve) == 104L,
            identical(as.integer(curve$year), rep(2018:2019, each = 52)),
            identical(as.integer(curve$week), rep(1:52, 2)),
            identical(as.integer(curve$step), 1:104),
            np_same(as.numeric(curve$observed), as.numeric(input$test$cases)),
            all(is.finite(as.matrix(curve[c("observed", "pred_median", "lower", "upper")]))),
            all(curve$lower <= curve$pred_median & curve$pred_median <= curve$upper))
  label <- "Humidity seasonal (no npieff)"
  curve$island <- island
  curve$model <- label
  list(curves = curve,
       metrics = cbind(island = island, model = label,
                        predicted$metrics[c("mean_log_score", "MAE", "coverage_95")]),
       diagnostics = cbind(island = island, model = label, ds))
}

hn_plot_comparison <- function(curves) {
  islands <- c("Hokkaido", "Shikoku", "Ryukyu")
  models <- c("Baseline + npieff", "Humidity seasonal (no npieff)")
  stopifnot(nrow(curves) > 0L, all(curves$island %in% islands), all(curves$model %in% models))
  curves$island <- factor(curves$island, levels = islands)
  curves$model <- factor(curves$model, levels = models)
  ggplot2::ggplot(curves, ggplot2::aes(x = step)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#90CAF9", alpha = .45) +
    ggplot2::geom_line(ggplot2::aes(y = observed), colour = "black", linewidth = .35) +
    ggplot2::geom_line(ggplot2::aes(y = pred_median), colour = "#0072B2", linewidth = .55) +
    ggplot2::geom_vline(xintercept = 52.5, linetype = "dotted") +
    ggplot2::facet_grid(island ~ model, scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(1,27,53,79,104),
                              labels = c("2018 W1","W27","2019 W1","W27","W52")) +
    ggplot2::labs(title = "Baseline + npieff vs seasonal humidity without npieff",
      subtitle = "Training: 2013-2017; conditional hindcast: 2018-2019",
      x = NULL, y = "Weekly reported cases",
      caption = "Black: observed; blue: predictive median; shading: 95% predictive interval.") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8))
}

