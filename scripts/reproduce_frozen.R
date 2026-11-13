args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6L || args[1L] != "--evaluation" ||
    args[3L] != "--mode" || args[5L] != "--out" ||
    !args[4L] %in% c("check", "h2", "h3", "mechanisms"))
  stop("Usage: Rscript scripts/reproduce_frozen.R --evaluation INPUT/result.rds --mode check|h2|h3|mechanisms --out NEW_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."))
evaluation <- normalizePath(args[2L], mustWork=TRUE)
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
            "R/phase27.R", "R/frozen_reproduction.R")) sys.source(f, globalenv())
model <- frozen_m2_load(evaluation,
  "45c5e155d900da17637dbb26a8d562af2315c2d863a78df05679b10728c65e18")
errors <- frozen_m2_check(model)
write.csv(data.frame(component=names(errors), max_absolute_error=unname(errors),
                     tolerance=1e-6, pass=errors<=1e-6),
          file.path(out, "model_integrity.csv"), row.names=FALSE)
mode <- args[4L]
for (policy in c("raw", "blend")) {
  active <- if (policy == "raw") model else frozen_m2_blend(model)
  if (mode %in% c("h2", "mechanisms")) {
    h2 <- frozen_m2_h2(active)
    write.csv(h2$steps, file.path(out, paste0("h2_steps_", policy, ".csv")),
              row.names=FALSE)
    write.csv(h2$seasons, file.path(out, paste0("h2_seasons_", policy, ".csv")),
              row.names=FALSE)
  }
  if (mode %in% c("h3", "mechanisms")) {
    h3 <- frozen_m2_h3(active)
    write.csv(h3$steps, file.path(out, paste0("h3_steps_", policy, ".csv")),
              row.names=FALSE)
    write.csv(h3$seasons, file.path(out, paste0("h3_seasons_", policy, ".csv")),
              row.names=FALSE)
    write.csv(h3$finite, file.path(out, paste0("h3_finite_", policy, ".csv")),
              row.names=FALSE)
  }
}
writeLines(c(paste("evaluation_sha256", model$source_sha256),
             paste("mode", mode), paste("R", R.version.string),
             paste("run_time", format(Sys.time(), tz="UTC", usetz=TRUE))),
           file.path(out, "manifest.txt"))
writeLines("COMPLETED", file.path(out, "STATUS"))
