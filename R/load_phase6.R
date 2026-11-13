# Definition-only loader; does not load fitted values or observations.
load_spline_sirs_phase6 <- function(project_root=".",envir=parent.frame()) {
  sys.source(file.path(project_root,"R/load_phase5.R"),envir=envir)
  get("load_spline_sirs_phase5",envir=envir)(project_root,envir)
  for(f in c("phase6_workflow.R","phase6_diagnostics.R","phase6_bayes.R","phase6_reporting.R","phase6_rescue.R","phase7b_gradient.R","phase7c_tangent.R","phase4_peaks.R","phase8_shift.R","phase8b_fine_grid.R","phase9_peak_sensitivity.R","phase10_peak_identity.R","phase10c_robustness.R","phase11_susceptibility_transmission.R","phase12_diagnostics.R"))
    sys.source(file.path(project_root,"R",f),envir=envir)
  invisible(NULL)
}
