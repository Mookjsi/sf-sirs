# Open-loop deterministic propagation. Future cases are absent from this API.
p3_forecast <- function(train,fit,future) {
  validate_phase2_data(train)
  if(!isTRUE(fit$valid) || !identical(fit$status,"CONVERGED_VALID")) stop("Forecast requires a valid converged endpoint")
  if("cases" %in% names(future) || "observed" %in% names(future)) stop("Future cases must not enter propagation")
  n <- length(train$cases);m <- length(future$weeks)
  if(!m || !identical(as.integer(future$indices),seq.int(n+1L,n+m))) stop("Future indices must immediately follow training")
  seasonal_weeks(future$weeks);finite_numeric(future$npieff,"future NPI",m)
  ordinal <- function(x) as.integer(x$year)*52L+as.integer(x$week)-1L
  marks <- rbind(p3_markers(train),future$markers)
  if(nrow(marks)!=n+m || any(diff(ordinal(marks))!=1L) ||
     !identical(as.integer(marks$week),as.integer(c(train$seasonal_week,future$weeks))))
    stop("Future calendar markers are discontinuous or misaligned",call.=FALSE)
  # Match the fit to precisely its training-only input, never to held-out outcomes.
  ctx <- fit$comparison_context
  if(!identical(ctx$cases,train$cases) || !identical(ctx$weeks,train$seasonal_week) ||
     !identical(ctx$npieff,train$npieff) || !identical(ctx$fixed,train$fixed_parameters))
    stop("Fit/training contract mismatch",call.=FALSE)
  p <- fit$evaluation$parameters;model <- ctx$model
  initial <- make_initial_state(p$S0,p$I0,train$fixed_parameters$pop)
  # Replay the historical prefix in the unchanged Phase 1 core, then continue.
  # The terminal mass need not equal pop; resetting it would change the recurrence.
  all_states <- simulate_weekly_sirs(fit$evaluation$beta_by_week,initial,
    c(train$fixed_parameters,list(omega=model$omega,delta=model$delta,rho=p$rho)),
    c(train$seasonal_week,future$weeks),c(train$npieff,future$npieff))
  if(!isTRUE(all.equal(all_states[seq_len(n),,drop=FALSE],fit$evaluation$states,tolerance=1e-10)))
    stop("Historical replay differs from fitted training states",call.=FALSE)
  ii <- seq.int(n+1L,n+m)
  list(type="open_loop_forecast",origin=n,boundary_state=all_states[n,c("S","I","R"),drop=FALSE],
    states=all_states[ii,,drop=FALSE],indices=future$indices,markers=future$markers,
    beta_by_week=fit$evaluation$beta_by_week,phi=p$phi,likelihood_include=rep(TRUE,m),
    future_input=future$source,method="Unchanged Phase 1 prefix replay and continuation; no case assimilation")
}
p3_score_forecast <- function(forecast,outcomes,policy) {
  if(!identical(as.integer(outcomes$indices),as.integer(forecast$indices)) ||
     !identical(outcomes$markers,forecast$markers)) stop("Scoring outcomes do not match forecast rows")
  scores <- p3_metrics(outcomes$cases,forecast$states$expected_cases,forecast$phi)
  c(forecast,list(scores=scores,peaks=p3_peaks(forecast$markers,outcomes$cases,
      forecast$states$expected_cases,policy)))
}
p3_outcomes <- function(data,ids) {
  list(indices=ids,cases=data$cases[ids],markers=p3_markers(data,ids))
}
