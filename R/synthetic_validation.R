# Explicit user-invoked synthetic generation and descriptive assessment, not pass/fail recovery.
synthetic_contract_data <- function(spec) {
  integer_scalar(spec$years,"years");integer_scalar(spec$start_year,"start_year")
  scalar_number(spec$population,"population",0,strict_lower=TRUE)
  n <- 52L*spec$years
  d <- data.frame(island=spec$region,year=rep(seq.int(spec$start_year,length.out=spec$years),each=52L),
                  week=rep(1:52,spec$years),cases=rep(0,n))
  data <- adapt_preprocessed_island(d,data.frame(island=spec$region,pop=spec$population),spec$region,spec$post_npi)
  data$provenance$synthetic_contract_only <- TRUE
  data
}
generate_synthetic_dataset <- function(config) {
  if(!identical(config$data$kind,"synthetic_spec")) stop("Synthetic specification required",call.=FALSE)
  data <- synthetic_contract_data(config$data); truth <- config$truth
  scalar_number(truth$phi,"truth phi",0,strict_lower=TRUE)
  basis <- new_cyclic_basis(config$model$K)
  beta <- spline_beta_by_week(truth$free_coordinates,truth$beta0,basis)
  states <- simulate_weekly_sirs(beta,make_initial_state(truth$S0,truth$I0,data$fixed_parameters$pop),
    c(data$fixed_parameters,list(omega=truth$omega,delta=truth$delta,rho=truth$rho)),
    data$seasonal_week,data$npieff)
  observations <- with_phase2_seed(config$seed,function()
    stats::rnbinom(length(data$cases),mu=states$expected_cases,size=truth$phi))
  if(any(!is.finite(observations))) stop("Nonfinite synthetic NB2 observations",call.=FALSE)
  data$cases <- observations;data$time_markers$cases <- observations
  data$provenance$synthetic_contract_only <- FALSE
  list(type="synthetic_dataset",data=data,truth=list(parameters=truth,beta_by_week=beta,
    states=states,basis=basis,seed=config$seed,
    rng=c("Mersenne-Twister","Inversion","Rejection"),
    initial_row="Mean zero sentinel, hence simulated y[1]=0; likelihood still excludes row 1"),
    scenario=config$scenario,config=config)
}
assess_synthetic_fit <- function(dataset,result) {
  if(!identical(dataset$type,"synthetic_dataset")) stop("Synthetic truth artifact required",call.=FALSE)
  runs <- if(identical(result$type,"multistart_fit")) result$runs else list(result)
  data <- dataset$data
  if(!all(vapply(runs,function(x) identical(x$comparison_context$cases,data$cases) &&
    identical(x$comparison_context$weeks,data$seasonal_week) &&
    identical(x$comparison_context$npieff,data$npieff) && identical(x$comparison_context$fixed,data$fixed_parameters) &&
    identical(x$comparison_context$mask,data$likelihood_include),logical(1)))) stop("Dataset/result mismatch",call.=FALSE)
  metrics <- lapply(seq_along(runs),function(i) {
    r <- runs[[i]]; e <- r$evaluation
    row <- data.frame(start=i,status=r$status,valid=isTRUE(r$valid),objective=NA_real_,
      loglik=NA_real_,penalty=NA_real_,beta0=NA_real_,phi=NA_real_,
      beta_rmse=NA_real_,expected_cases_rmse=NA_real_,S0=NA_real_,rho=NA_real_)
    if(!is.null(e)) {row$objective <- e$objective;row$loglik <- e$loglik;row$penalty <- e$penalty}
    if(isTRUE(r$valid) && isTRUE(e$valid)) {
      row$beta0 <- e$parameters$beta0;row$phi <- e$parameters$phi
      row$S0 <- e$parameters$S0;row$rho <- e$parameters$rho
      row$beta_rmse <- sqrt(mean((e$beta_by_week-dataset$truth$beta_by_week)^2))
      ids <- data$likelihood_include
      row$expected_cases_rmse <- sqrt(mean((e$states$expected_cases[ids]-dataset$truth$states$expected_cases[ids])^2))
    }
    row
  })
  metrics <- do.call(rbind,metrics)
  between <- data.frame(start_a=integer(),start_b=integer(),status_a=character(),status_b=character(),
    objective_difference=numeric(),beta0_difference=numeric(),phi_difference=numeric(),
    S0_difference=numeric(),rho_difference=numeric(),beta_curve_rmse=numeric(),expected_curve_rmse=numeric())
  ids <- which(metrics$valid)
  if(length(ids)>=2L) {
    pairs <- utils::combn(ids,2L)
    between <- do.call(rbind,lapply(seq_len(ncol(pairs)),function(i) {
      a <- pairs[1,i];b <- pairs[2,i];x <- runs[[a]]$evaluation;y <- runs[[b]]$evaluation
      data.frame(start_a=a,start_b=b,status_a=runs[[a]]$status,status_b=runs[[b]]$status,
        objective_difference=x$objective-y$objective,beta0_difference=x$parameters$beta0-y$parameters$beta0,
        phi_difference=x$parameters$phi-y$parameters$phi,
        S0_difference=x$parameters$S0-y$parameters$S0,rho_difference=x$parameters$rho-y$parameters$rho,
        beta_curve_rmse=sqrt(mean((x$beta_by_week-y$beta_by_week)^2)),
        expected_curve_rmse=sqrt(mean((x$states$expected_cases[data$likelihood_include]-
                                      y$states$expected_cases[data$likelihood_include])^2)))
    }))
  }
  list(type="synthetic_assessment",per_start=metrics,status_counts=table(metrics$status),
       invalid_fraction=mean(!metrics$valid),between_starts=between,truth=dataset$truth$parameters,
       interpretation="Descriptive recovery on one realization; no strict statistical pass threshold or identifiability claim")
}
