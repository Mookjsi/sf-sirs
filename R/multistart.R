# Sequential starts, explicit RNG seed, and separate converged/nonconverged selections.
with_phase2_seed <- function(seed,fun) {
  integer_scalar(seed,"seed",0L)
  old_kind <- RNGkind(); existed <- exists(".Random.seed",envir=globalenv(),inherits=FALSE)
  if(existed) old_seed <- get(".Random.seed",envir=globalenv())
  on.exit({do.call(RNGkind,as.list(old_kind));
    if(existed) assign(".Random.seed",old_seed,envir=globalenv()) else
      if(exists(".Random.seed",envir=globalenv(),inherits=FALSE)) rm(".Random.seed",envir=globalenv())},add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)
  fun()
}
validate_starts_config <- function(starts_config) {
  if(!is.list(starts_config)) stop("Set starts_config explicitly for multistart",call.=FALSE)
  scalar_number(starts_config$perturb_sd,"perturb_sd",0)
  integer_scalar(starts_config$n_perturbed,"n_perturbed",0L)
  integer_scalar(starts_config$seed,"starts seed",0L)
  for(n in c("beta0_candidates","phi_candidates")) {
    finite_numeric(starts_config[[n]],n)
    if(any(starts_config[[n]]<=0)) stop(n," must be positive",call.=FALSE)
  }
  if(!is.null(starts_config$user_starts) && !is.list(starts_config$user_starts))
    stop("user_starts must be a list of physical parameter lists",call.=FALSE)
  invisible(starts_config)
}
generate_starts <- function(base_parameters,problem,starts_config) {
  validate_starts_config(starts_config)
  with_phase2_seed(starts_config$seed,function() {
    starts <- list(); grid <- expand.grid(beta0=unique(starts_config$beta0_candidates),phi=unique(starts_config$phi_candidates))
    perturbations <- if(isTRUE(problem$map$constant)) 0L else starts_config$n_perturbed
    for(i in seq_len(nrow(grid))) for(j in 0:perturbations) {
      p <- base_parameters;p$beta0 <- grid$beta0[i];p$phi <- grid$phi[i]
      p$free_coordinates <- if(j==0L) rep(0,problem$basis$K-1L) else
        stats::rnorm(problem$basis$K-1L,sd=starts_config$perturb_sd)
      starts[[length(starts)+1L]] <- list(label=paste0("grid",i,"-",if(j==0) "constant" else paste0("perturbed",j)),parameters=p)
    }
    if(!is.null(starts_config$user_starts)) for(i in seq_along(starts_config$user_starts))
      starts[[length(starts)+1L]] <- list(label=paste0("user",i),parameters=starts_config$user_starts[[i]])
    starts
  })
}
select_multistart <- function(runs) {
  if(!length(runs)) stop("No starts provided",call.=FALSE)
  context <- runs[[1]]$comparison_context
  if(!all(vapply(runs,function(x) identical(x$comparison_context,context),logical(1))))
    stop("Do not compare different data/K/lambda/fixed-parameter contexts",call.=FALSE)
  eligible <- vapply(runs,function(x) isTRUE(x$valid) && !isTRUE(x$reported_sentinel) &&
    !is.null(x$evaluation) && isTRUE(x$evaluation$valid) && is.finite(x$evaluation$objective),logical(1))
  converged <- which(eligible & vapply(runs,function(x) identical(x$status,"CONVERGED_VALID"),logical(1)))
  nonconverged <- which(eligible & vapply(runs,function(x) identical(x$status,"NONCONVERGED_VALID"),logical(1)))
  best <- function(ids) if(length(ids)) ids[which.min(vapply(runs[ids],function(x) x$evaluation$objective,numeric(1)))] else NA_integer_
  ci <- best(converged); ni <- best(nonconverged)
  list(status=if(!is.na(ci)) "CONVERGED_CANDIDATE_SELECTED" else if(!is.na(ni)) "NO_CONVERGED_CANDIDATE" else "ALL_FAILED",
       best_converged_index=ci,best_nonconverged_index=ni,
       representative=if(!is.na(ci)) runs[[ci]] else NULL,
       best_nonconverged=if(!is.na(ni)) runs[[ni]] else NULL)
}
fit_multistart <- function(data,model_config,estimation_config,starts_config,initial_parameters,optimizer=stats::optim,explicit_starts=NULL) {
  problem <- prepare_estimation_problem(data,model_config,estimation_config)
  if(is.null(explicit_starts)) starts <- generate_starts(initial_parameters,problem,starts_config) else {
    validate_starts_config(starts_config)
    if(!is.list(explicit_starts) || !length(explicit_starts)) stop("Explicit starts must be a nonempty list")
    starts <- explicit_starts
    labels <- vapply(starts,function(x) {
      if(!is.list(x$parameters) || !is.character(x$label) || length(x$label)!=1L || is.na(x$label) || !nzchar(x$label))
        stop("Every explicit start needs a label and physical parameter list")
      x$label
    },character(1))
    if(anyDuplicated(labels)) stop("Explicit start labels must be unique")
  }
  runs <- vector("list",length(starts))
  for(i in seq_along(starts)) {
    runs[[i]] <- fit_prepared_problem(problem,starts[[i]]$parameters,optimizer)
    runs[[i]]$start_label <- starts[[i]]$label
  }
  c(list(type="multistart_fit",runs=runs,starts=starts,seed=starts_config$seed,
         comparison_context=problem$comparison_context,settings=estimation_config),select_multistart(runs))
}
