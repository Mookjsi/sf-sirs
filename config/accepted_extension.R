# Extension-only choices; accepted historical parameters come from the saved fit.
accepted_extension_config <- function() list(version="accepted-extension-1",
  islands="Honshu",future_years=2020:2023,
  expected_weeks=stats::setNames(c(53L,52L,52L,52L),2020:2023),
  week53_policy="repeat_week52",week53_sensitivity="map_to_week1",
  future_external_multiplier=1,
  q_bounds=log(c(.05,20)),starts=list(seed=20261110L,
    log_q=list(c(0,0,0),c(log(.5),0,0),c(0,log(.5),0),c(0,0,log(.5)))),
  optimizer=list(maxit=2500L,factr=1e7,pgtol=0,ndeps=rep(1e-4,3L)),
  replay_tolerance=1e-7,stationarity_step=1e-4,
  near_best_objective_gap=2)
