# Penalized maximum likelihood using Phase 1 functions, without priors or Jacobians.
evaluate_candidate <- function(q,problem) {
  answer <- list(valid=FALSE,parameters_valid=FALSE,loglik=NA_real_,penalty=NA_real_,
    objective=Inf,n_observations=problem$n_observations,failure_stage=NULL,
    failure_reason=NULL,state_transition_failure=NULL,parameters=NULL,states=NULL,beta_by_week=NULL)
  stage <- "parameters"
  tryCatch({
    p <- decode_parameters(q,problem$map); answer$parameters <- p; answer$parameters_valid <- TRUE
    stage <- if(isTRUE(problem$map$fourier)) "fourier" else "spline"
    profile <- if(isTRUE(problem$map$constant)) list(beta_by_week=rep(p$beta0,52),
      coefficients=rep(0,problem$basis$K),log_normalizer=0) else if(isTRUE(problem$map$fourier))
      fourier_profile(p$free_coordinates,p$beta0,problem$basis) else spline_profile(p$free_coordinates,p$beta0,problem$basis)
    answer$beta_by_week <- profile$beta_by_week
    if(isTRUE(problem$map$fourier)) {
      answer$fourier_coefficients <- setNames(profile$coefficients,problem$basis$coefficient_names)
      answer$a0 <- profile$a0
    } else answer$spline_coefficients <- profile$coefficients
    answer$log_normalizer <- profile$log_normalizer
    stage <- "penalty"
    answer$penalty <- if(isTRUE(problem$map$fourier)) 0 else cyclic_penalty(profile$coefficients,problem$model$lambda)
    stage <- "initial_state"
    initial <- make_initial_state(p$S0,p$I0,problem$data$fixed_parameters$pop)
    stage <- "state_transition"
    # Reuse the Phase 1 recurrence; no duplicated S/I/R equations.
    dynamic <- if(isTRUE(problem$map$joint_mode)) list(omega=p$omega,delta=p$delta,rho=p$rho,mu=p$mu,gamma=p$gamma) else
      list(omega=problem$model$omega,delta=problem$model$delta,rho=p$rho)
    fixed <- utils::modifyList(problem$data$fixed_parameters,dynamic)
    answer$states <- simulate_weekly_sirs(profile$beta_by_week,initial,fixed,
      problem$data$seasonal_week,problem$data$npieff)
    stage <- "likelihood"
    ll <- nb2_loglik(problem$data$cases,answer$states$expected_cases,p$phi,problem$data$likelihood_include)
    answer$loglik <- ll$total; answer$per_week_loglik <- ll$per_week
    answer$included_per_week <- ll$included_per_week
    answer$objective <- -answer$loglik+answer$penalty
    stage <- "objective_range"
    if(!is.finite(answer$objective) || answer$objective>=problem$estimation$failure_value)
      stop("Nonfinite objective or objective reaches the explicit optimizer failure sentinel",call.=FALSE)
    answer$valid <- TRUE
    answer
  },error=function(e) {
    answer$failure_stage <- stage; answer$failure_reason <- conditionMessage(e)
    if(stage=="state_transition") answer$state_transition_failure <- conditionMessage(e)
    answer
  })
}
make_objective <- function(problem) {
  tracker <- new.env(parent=emptyenv())
  tracker$evaluations <- 0L; tracker$invalid <- 0L; tracker$reasons <- integer(); tracker$examples <- list()
  fn <- function(q) {
    tracker$evaluations <- tracker$evaluations+1L
    d <- evaluate_candidate(q,problem)
    if(d$valid) return(d$objective)
    tracker$invalid <- tracker$invalid+1L
    key <- paste(d$failure_stage,d$failure_reason,sep=": ")
    if(!key %in% names(tracker$reasons)) tracker$reasons[key] <- 0L
    tracker$reasons[key] <- tracker$reasons[key]+1L
    if(length(tracker$examples)<problem$estimation$failure_examples)
      tracker$examples[[length(tracker$examples)+1L]] <- list(coordinates=q,stage=d$failure_stage,reason=d$failure_reason)
    problem$estimation$failure_value
  }
  list(fn=fn,diagnostics=function() list(evaluations=tracker$evaluations,invalid=tracker$invalid,
                                        reasons=tracker$reasons,examples=tracker$examples))
}
