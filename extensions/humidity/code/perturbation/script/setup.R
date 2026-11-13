# Runtime setup only; the original scientific models remain in their original files.
perturbation_setup <- function(folder = NULL) {
  root <- Sys.getenv("PERTURBATION_ROOT", unset = Sys.getenv("PERTURBATION_ROOT"))
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  lib <- file.path(root, ".R-library")
  if (!dir.exists(lib)) stop("Project packages are missing: ", lib)
  .libPaths(c(lib, .Library))
  Sys.setenv(PERTURBATION_ROOT = root, R_LIBS = lib, R_LIBS_USER = lib)
  options(perturbation.root = root)
  expected <- c(rstan = "2.32.7", StanHeaders = "2.32.10")
  for (pkg in names(expected)) {
    disk_version <- as.character(utils::packageVersion(pkg, lib.loc = lib))
    if (disk_version != expected[[pkg]])
      stop(pkg, " must be ", expected[[pkg]], " in the project library; found ", disk_version)
    if (pkg %in% loadedNamespaces() &&
        as.character(getNamespaceVersion(pkg)) != expected[[pkg]])
      stop("An older session has already loaded ", pkg,
           ". Restart R in perturbation.Rproj, then run this file from its first line.")
  }
  if (!is.null(folder)) setwd(file.path(root, folder))
  invisible(root)
}

perturbation_validate <- function(x) {
  stopifnot(length(x$N) == 1L, x$N > 1L,
            length(x$Nnonpi) == 1L, x$Nnonpi >= 1L, x$Nnonpi < x$N,
            length(x$cases) == x$N, length(x$week) == x$N,
            all(is.finite(x$cases)), all(x$cases >= 0),
            all(x$cases == floor(x$cases)),
            all(x$week %in% 1:52),
            length(x$pop) == 1L, is.finite(x$pop), x$pop > 0,
            length(x$mu) == 1L, is.finite(x$mu), x$mu >= 0,
            length(x$gamma) == 1L, is.finite(x$gamma), x$gamma > 0)
  if (!is.null(x$type)) stopifnot(length(x$type) == x$N, all(x$type %in% 1:2))
  invisible(x)
}
