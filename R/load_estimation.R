# Explicit user loader. It only installs function definitions into the chosen environment.
load_spline_sirs_estimation <- function(project_root=".",envir=parent.frame()) {
  sys.source(file.path(project_root,"R/load_modules.R"),envir=envir)
  get("load_spline_sirs",envir=envir)(project_root,envir=envir)
  files <- c("estimation_config.R","parameter_transforms.R","penalized_objective.R","fit_single.R",
             "multistart.R","synthetic_validation.R","result_storage.R")
  for(f in files) sys.source(file.path(project_root,"R",f),envir=envir)
  invisible(files)
}
