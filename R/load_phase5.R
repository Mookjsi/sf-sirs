# Definition-only loader. No data or fits are loaded by source().
load_spline_sirs_phase5 <- function(project_root=".",envir=parent.frame()) {
  sys.source(file.path(project_root,"R/load_phase3.R"),envir=envir)
  get("load_spline_sirs_phase3",envir=envir)(project_root,envir)
  for(f in c("phase5_workflow.R","phase5_reporting.R")) sys.source(file.path(project_root,"R",f),envir=envir)
  invisible(NULL)
}
