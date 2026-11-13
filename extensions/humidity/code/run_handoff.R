# Run one explicitly selected research entry point from the handoff root.
# This launcher does not select an experiment or change model settings.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop('Usage: Rscript --vanilla code/run_handoff.R "code/path/to/script.R"')
}
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "ARCHIVE_MANIFEST.csv"))) {
  stop("Start this command in the extracted sf-sirs-handoff directory.")
}
entry <- normalizePath(file.path(root, args[[1]]), winslash = "/", mustWork = TRUE)
if (!startsWith(entry, paste0(root, "/code/")) || !grepl("\\.R$", entry)) {
  stop("The entry point must be an R script inside this handoff's code directory.")
}
project <- file.path(root, "code", "perturbation")
analysis <- file.path(root, "code", "analysis")
Sys.setenv(HANDOFF_ROOT = root, PERTURBATION_ROOT = project,
           HANDOFF_ANALYSIS_ROOT = analysis)
lib <- file.path(project, ".R-library")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
message("Running selected entry point: ", entry)
source(entry, encoding = "UTF-8", chdir = FALSE)
