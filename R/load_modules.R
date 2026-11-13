# Explicit user-invoked loader. It defines functions only and does not fit or simulate.
load_spline_sirs <- function(project_root=".",envir=parent.frame()) {
  files <- c("validation.R","data_adapter.R","sirs_core.R","nb2.R","cyclic_spline.R","fourier.R","model_interface.R")
  for(f in files) sys.source(file.path(project_root,"R",f),envir=envir)
  invisible(files)
}
