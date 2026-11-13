# Single-start fitting. Source() only defines functions; the user explicitly invokes fitting.
endpoint_status <- function(valid,convergence,reported_sentinel=FALSE) {
  if(!isTRUE(valid) || isTRUE(reported_sentinel)) return("INVALID")
  if(!is.numeric(convergence) || length(convergence)!=1L || !is.finite(convergence)) return("OPTIMIZER_ERROR")
  if(convergence==0L) "CONVERGED_VALID" else "NONCONVERGED_VALID"
}
fit_prepared_problem <- function(problem,initial_parameters,optimizer=stats::optim) {
  if(isTRUE(problem$data$provenance$synthetic_contract_only))
    stop("Contract placeholders are not observations; generate an explicit dataset first",call.=FALSE)
  result <- list(type="single_fit",status="INVALID_START",valid=FALSE,initial_parameters=initial_parameters,
    initial_coordinates=NULL,coordinates=NULL,convergence=NA_integer_,message=NULL,
    counts=NULL,warnings=character(),evaluation=NULL,reported_sentinel=FALSE,
    comparison_context=problem$comparison_context,optimizer_settings=problem$estimation,
    map=problem$map,basis=problem$basis,lambda=problem$model$lambda,
    seed=NULL,seed_usage="No randomness in the supported single-start optim methods")
  q <- tryCatch(encode_parameters(initial_parameters,problem$map),error=function(e) e)
  if(inherits(q,"error")) {result$message <- conditionMessage(q);return(result)}
  result$initial_coordinates <- q
  result$initial_parameter_record <- parameter_record(problem,decode_parameters(q,problem$map))
  start <- evaluate_candidate(q,problem); result$initial_evaluation <- start
  if(!start$valid) {result$message <- start$failure_reason;return(result)}
  objective <- make_objective(problem)
  control <- c(problem$estimation$control,list(ndeps=problem$ndeps,parscale=problem$parscale))
  args <- list(par=q,fn=objective$fn,method=problem$estimation$method,control=control,hessian=FALSE)
  if(problem$estimation$method=="L-BFGS-B") {args$lower <- problem$map$lower;args$upper <- problem$map$upper}
  warnings <- character()
  optimized <- tryCatch(withCallingHandlers(do.call(optimizer,args),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e) e)
  result$warnings <- warnings; result$failure_diagnostics <- objective$diagnostics()
  result$effective_optim_control <- control
  if(inherits(optimized,"error")) {result$status <- "OPTIMIZER_ERROR";result$message <- conditionMessage(optimized);return(result)}
  if(!is.list(optimized) || is.null(optimized$par)) {
    result$status <- "OPTIMIZER_ERROR";result$message <- "Optimizer returned no parameter vector";return(result)
  }
  result$coordinates <- optimized$par; result$convergence <- optimized$convergence
  result$message <- optimized$message; result$counts <- optimized$counts; result$optimizer_value <- optimized$value
  final <- evaluate_candidate(optimized$par,problem)
  sentinel <- length(optimized$value)!=1L || !is.finite(optimized$value) || optimized$value>=problem$estimation$failure_value
  result$reported_sentinel <- sentinel
  result$evaluation <- final; result$valid <- final$valid && !sentinel
  result$status <- endpoint_status(final$valid,optimized$convergence,sentinel)
  if(!final$valid) result$message <- final$failure_reason
  boundary_evaluated <- is.numeric(optimized$par) && length(optimized$par)==problem$map$dimension && all(is.finite(optimized$par))
  hits <- if(boundary_evaluated) problem$map$order[
    abs(optimized$par-problem$map$lower)<=problem$estimation$boundary_tolerance |
    abs(optimized$par-problem$map$upper)<=problem$estimation$boundary_tolerance] else character()
  result$boundary <- list(coordinates=hits,evaluated=boundary_evaluated,
    lower=problem$map$lower,upper=problem$map$upper,tolerance=problem$estimation$boundary_tolerance)
  if(final$parameters_valid) result$parameter_record <- parameter_record(problem,final$parameters)
  result$interpretation <- "Local penalized likelihood endpoint only; no global optimum, identifiability or uncertainty claim"
  result
}
fit_spline_sirs <- function(data,model_config,estimation_config,initial_parameters,optimizer=stats::optim) {
  problem <- prepare_estimation_problem(data,model_config,estimation_config)
  fit_prepared_problem(problem,initial_parameters,optimizer)
}
