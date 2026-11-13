args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8L || !identical(args[c(1L,3L,5L,7L)],
    c("--evaluation", "--grid", "--parents", "--out")))
  stop("Usage: Rscript scripts/reproduce_h1.R --evaluation INPUT/result.rds --grid INPUT/peak_path_grid.csv --parents INPUT/h1_event_summary.csv --out NEW_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."))
paths <- vapply(args[c(2L,4L,6L)], normalizePath, character(1), mustWork=TRUE)
expected_hashes <- c("45c5e155d900da17637dbb26a8d562af2315c2d863a78df05679b10728c65e18",
                     "7a8913d47d1a98be409981709d23e5530f7130f83079707280c37be8e566c6dd",
                     "5fc580cce9e481694eb1f9d2fe30dc30081df2dd004b687e824e48ef1f68880e")
if (!identical(unname(tools::sha256sum(paths)), expected_hashes))
  stop("H1 source input hash mismatch")
out <- args[8L]
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
            "R/phase27.R", "R/m2_hypotheses.R", "R/frozen_reproduction.R",
            "config/phase25.R")) sys.source(f, globalenv())
model <- frozen_m2_load(paths[1L],
  "45c5e155d900da17637dbb26a8d562af2315c2d863a78df05679b10728c65e18")
frozen_m2_check(model)
grid <- read.csv(paths[2L], stringsAsFactors=FALSE)
parents <- read.csv(paths[3L], stringsAsFactors=FALSE)
result <- frozen_m2_h1(model, grid, parents)
keys <- c("policy", "season", "parent")
comparison <- merge(result$parents, parents, by=keys,
                    suffixes=c(".new", ".source"), sort=TRUE)
fields <- c("timing_event", "displacement", "exchange", "unresolved",
            "mixed", "absolute_movement")
if (nrow(comparison) != 1720L || any(vapply(fields, function(field) {
  a <- comparison[[paste0(field, ".new")]]
  b <- comparison[[paste0(field, ".source")]]
  if (is.numeric(a)) any(abs(a-b) > 0) else any(a != b)
}, logical(1)))) stop("Independent H1 parent classifications differ from final source")
write.csv(result$edges, file.path(out, "h1_edges.csv"), row.names=FALSE)
write.csv(result$parents, file.path(out, "h1_parents.csv"), row.names=FALSE)
write.csv(data.frame(parent_rows=nrow(comparison),
                     parent_classifications_match=TRUE),
          file.path(out, "h1_source_comparison.csv"), row.names=FALSE)
writeLines(c(paste("evaluation_sha256", model$source_sha256),
             paste("grid_sha256", unname(tools::sha256sum(paths[2L]))),
             paste("parents_sha256", unname(tools::sha256sum(paths[3L]))),
             paste("R", R.version.string)), file.path(out, "manifest.txt"))
writeLines("COMPLETED", file.path(out, "STATUS"))
