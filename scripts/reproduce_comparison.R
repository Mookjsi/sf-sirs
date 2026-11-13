args <- commandArgs(trailingOnly=TRUE)
if (!length(args) %in% c(6L, 8L) || !identical(args[c(1L,3L,5L)],
    c("--evaluation", "--park-weekly", "--out")) ||
    (length(args) == 8L && args[7L] != "--history"))
  stop("Usage: Rscript scripts/reproduce_comparison.R --evaluation INPUT/m2-evaluation.rds --park-weekly INPUT/park-weekly-fit.csv --out NEW_DIR [--history INPUT/phase28-history]")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."))
paths <- vapply(args[c(2L,4L)], normalizePath, character(1), mustWork=TRUE)
if (!identical(unname(tools::sha256sum(paths)), c(
    "45c5e155d900da17637dbb26a8d562af2315c2d863a78df05679b10728c65e18",
    "b4f864089114c2febeac795dec56fd69fc516b3013dbc8ba9c46a5408ee74d79")))
  stop("Unexpected M2 or Park summary input")
out <- args[6L]
if (file.exists(out)) stop("Output path already exists: ", out)
dir.create(out, recursive=TRUE)
out <- normalizePath(out)
setwd(root)
sys.source("R/load_phase6.R", globalenv())
load_spline_sirs_phase6(root, globalenv())
for (f in c("R/phase5_reporting.R", "R/phase10_peak_identity.R",
            "R/result_storage.R", "R/accepted_extension.R", "R/phase21.R",
            "R/phase22.R", "R/phase22_analysis.R", "R/phase24.R",
            "R/phase25.R", "R/phase26.R", "R/phase16_targeted.R",
            "R/phase27.R", "R/phase28.R", "R/frozen_reproduction.R"))
  sys.source(f, globalenv())
model <- frozen_m2_load(paths[1L])
errors <- frozen_m2_check(model)
park <- read.csv(paths[2L], stringsAsFactors=FALSE)
w <- model$weekly
if (nrow(park) != nrow(w) || nrow(park) != 573L ||
    !all(park[c("row", "year", "week", "seasonal_week", "observed")] ==
         w[c("row", "year", "week", "seasonal_week", "observed")]) ||
    any(!is.finite(park$expected_mean)) ||
    any(!is.finite(park$expected_median)))
  stop("Park and M2 common panel alignment failed")
clock <- frozen_m2_clock(model)
curves <- list(M2=w$expected, PARK_median=park$expected_median,
               PARK_mean=park$expected_mean)
metrics <- do.call(rbind, lapply(names(curves), function(k)
  p28_curve_metrics(w$observed, curves[[k]], model$score$mask,
    clock, if (k == "M2") "M2" else "PARK", k)))
timing <- do.call(rbind, lapply(names(curves), function(k)
  p28_timing(w$observed, curves[[k]], clock,
    if (k == "M2") "M2" else "PARK", k)))
movement <- do.call(rbind, lapply(split(timing,
  interaction(timing$model, timing$summary, drop=TRUE)), function(x) {
    a <- x[x$season %in% c(2015L, 2017L), ]
    a <- a[order(a$season), ]
    if (nrow(a) != 2L || !all(a$complete)) stop("Missing transition seasons")
    data.frame(model=a$model[1L], summary=a$summary[1L],
      observed_peak_movement=diff(a$observed_peak),
      fitted_peak_movement=diff(a$fitted_peak),
      peak_movement_error=diff(a$fitted_peak)-diff(a$observed_peak),
      observed_cog_movement=diff(a$observed_cog),
      fitted_cog_movement=diff(a$fitted_cog),
      cog_movement_error=diff(a$fitted_cog)-diff(a$observed_cog))
  }))
write.csv(metrics, file.path(out, "fit_metrics.csv"), row.names=FALSE)
write.csv(timing, file.path(out, "peak_metrics.csv"), row.names=FALSE)
write.csv(movement, file.path(out, "movement_metrics.csv"), row.names=FALSE)
write.csv(data.frame(component=names(errors), max_absolute_error=unname(errors)),
          file.path(out, "model_integrity.csv"), row.names=FALSE)
grDevices::pdf(file.path(out, "observed_fitted.pdf"), width=10, height=7)
graphics::par(mfrow=c(2,1))
for (years in list(2013:2019, 2020:2023)) {
  ix <- which(w$year %in% years)
  graphics::plot(ix, w$observed[ix], type="l", col="black",
                 xlab="Common panel row", ylab="Reported weekly cases")
  graphics::lines(ix, curves$M2[ix], col="steelblue")
  graphics::lines(ix, curves$PARK_median[ix], col="firebrick")
}
graphics::legend("topright", legend=c("Observed", "M2", "Park posterior median"),
                 col=c("black", "steelblue", "firebrick"), lty=1)
grDevices::dev.off()
if (length(args) == 8L) {
  history <- normalizePath(args[8L], mustWork=TRUE)
  historical_paths <- file.path(history,
    c("runtime_metrics.csv", "parameter_inventory.csv"))
  if (!identical(unname(tools::sha256sum(historical_paths)), c(
      "ddd3aa6ede80d1396a7de5204adcff71f33a6446817ff863cac73ec85cd686cd",
      "c9a5111d1ea23befd0ea70768d28337416177d793b3c08f4e764cb99a8d1f6bf")))
    stop("Historical runtime or coordinate inventory hash mismatch")
  runtime <- read.csv(historical_paths[1L])
  inventory <- read.csv(historical_paths[2L])
  file.copy(historical_paths, out)
  grDevices::pdf(file.path(out, "accuracy_runtime_coordinates.pdf"),
                 width=10, height=8)
  graphics::par(mfrow=c(2,2), mar=c(5,4,3,1))
  periods <- c("2013-2019", "2015-2017", "2020-2023", "full")
  selected <- function(model, summary) {
    q <- metrics[metrics$model == model & metrics$summary == summary, ]
    q[match(periods, q$period), ]
  }
  m2 <- selected("M2", "M2")
  median <- selected("PARK", "PARK_median")
  mean <- selected("PARK", "PARK_mean")
  graphics::barplot(rbind(m2$mae, median$mae), beside=TRUE,
    names.arg=periods, las=2, col=c("steelblue", "firebrick"),
    main="MAE: Park posterior median", ylab="Cases")
  graphics::barplot(rbind(m2$rmse, mean$rmse), beside=TRUE,
    names.arg=periods, las=2, col=c("steelblue", "firebrick"),
    main="RMSE: Park posterior mean", ylab="Cases")
  graphics::barplot(runtime$total_seconds, names.arg=runtime$model,
                    col=c("steelblue", "firebrick"), main="Original wall seconds")
  graphics::barplot(inventory$free_coordinates, names.arg=inventory$model,
                    col=c("steelblue", "firebrick"), main="Free coordinates")
  grDevices::dev.off()
}
writeLines(c(paste("m2_sha256", unname(tools::sha256sum(paths[1L]))),
             paste("park_summary_sha256", unname(tools::sha256sum(paths[2L]))),
             paste("R", R.version.string)), file.path(out, "manifest.txt"))
writeLines("COMPLETED", file.path(out, "STATUS"))
