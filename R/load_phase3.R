# Definitions only; no data, fitting, tests or graphics on source.
load_spline_sirs_phase3 <- function(project_root=".",envir=parent.frame()) {
  sys.source(file.path(project_root,"R/load_estimation.R"),envir=envir)
  get("load_spline_sirs_estimation",envir=envir)(project_root,envir)
  for(f in c("phase3_data.R","phase3_contracts.R","phase3_metrics.R","phase3_forecast.R","phase3_selection.R",
             "phase3_reference.R","phase3_workflow.R","phase3_reporting.R","phase3_report_outputs.R"))
    sys.source(file.path(project_root,"R",f),envir=envir)
  invisible(NULL)
}
