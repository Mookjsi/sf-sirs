# Hyperparameter selection uses B only. There is no C argument in this module.
p3_candidates <- function(candidates) {
  if(!is.data.frame(candidates) || !all(c("id","K","lambda") %in% names(candidates)) ||
     !nrow(candidates) || nrow(candidates)>9L || anyNA(candidates)) stop("Supply 1..9 explicit candidates")
  if(!is.character(candidates$id) || any(!nzchar(candidates$id)) || anyDuplicated(candidates$id)) stop("Unique candidate IDs required")
  for(i in seq_len(nrow(candidates))) {
    integer_scalar(candidates$K[i],"candidate K",4L);scalar_number(candidates$lambda[i],"candidate lambda",0)
  }
  if(anyDuplicated(candidates[c("K","lambda")])) stop("Duplicate K/lambda candidates")
  invisible(candidates)
}
p3_model <- function(cfg,candidate=NULL,constant=FALSE) {
  m <- cfg$model;m$transmission <- if(constant) "constant" else "spline"
  m$K <- if(constant) 4L else candidate$K[1]
  m$lambda <- if(constant) 0 else candidate$lambda[1]
  # K=4 is an internal zero-profile carrier for M0, never four fitted shape coefficients.
  m$fixed_sources <- lapply(cfg$fixed_provenance,function(x) paste(x$category,x$citation,sep=": "))
  m
}
p3_fit <- function(train,cfg,candidate=NULL,constant=FALSE,optimizer=stats::optim) {
  model <- p3_model(cfg,candidate,constant)
  init <- cfg$initial_parameters;init$free_coordinates <- rep(0,model$K-1L)
  starts <- cfg$starts
  # Constant starts vary beta0/phi only. User spline starts are not silently mapped to M0.
  if(constant) {
    starts$n_perturbed <- 0L
    if(length(starts$user_starts)) for(i in seq_along(starts$user_starts)) {
      if(length(starts$user_starts[[i]]$free_coordinates) && any(starts$user_starts[[i]]$free_coordinates!=0))
        stop("For matched Phase 3 use beta0/phi grids, not spline-specific user starts")
      starts$user_starts[[i]]$free_coordinates <- rep(0,model$K-1L)
    }
  }
  fit_multistart(train,model,cfg$estimation,starts,init,optimizer)
}
p3_candidate_campaign <- function(train,future,outcomes,cfg,optimizer=stats::optim) {
  # Callers pass A and B separately; C cannot enter fits or selection.
  if(!identical(outcomes$indices,future$indices)) stop("B outcome indices mismatch")
  p3_candidates(cfg$candidates)
  records <- vector("list",nrow(cfg$candidates))
  for(i in seq_len(nrow(cfg$candidates))) {
    candidate <- cfg$candidates[i,,drop=FALSE]
    record <- list(candidate=candidate,status="FAILED",fit=NULL,forecast=NULL,reason=NULL)
    records[[i]] <- tryCatch({
      record$fit <- p3_fit(train,cfg,candidate,optimizer=optimizer)
      record$status <- record$fit$status
      if(is.null(record$fit$representative)) {
        record$reason <- "No converged valid start";record
      } else {
        record$forecast <- p3_score_forecast(p3_forecast(train,record$fit$representative,future),outcomes,cfg$peaks)
        if(!is.finite(record$forecast$scores$summary$mean_nb2_nll)) stop("Nonfinite B predictive score")
        record$status <- "AVAILABLE";record
      }
    },error=function(e) {record$status <- "FAILED";record$reason <- conditionMessage(e);record})
  }
  records
}
p3_baseline_record <- function(train,future,outcomes,cfg,optimizer=stats::optim) {
  record <- list(model="M0",status="FAILED",fit=NULL,forecast=NULL,reason=NULL)
  tryCatch({
    record$fit <- p3_fit(train,cfg,constant=TRUE,optimizer=optimizer)
    record$status <- record$fit$status
    if(is.null(record$fit$representative)) {
      record$reason <- "No converged valid M0 start";return(record)
    }
    record$forecast <- p3_score_forecast(p3_forecast(train,record$fit$representative,future),outcomes,cfg$peaks)
    if(!is.finite(record$forecast$scores$summary$mean_nb2_nll)) stop("Nonfinite M0 B score")
    record$status <- "AVAILABLE";record
  },error=function(e) {record$status <- "FAILED";record$reason <- conditionMessage(e);record})
}
p3_candidate_table <- function(records) {
  do.call(rbind,lapply(records,function(x) {
    fit <- x$fit
    data.frame(id=x$candidate$id,K=x$candidate$K,lambda=x$candidate$lambda,status=x$status,
      starts=if(is.null(fit)) 0L else length(fit$runs),
      converged=if(is.null(fit)) 0L else sum(vapply(fit$runs,function(r) identical(r$status,"CONVERGED_VALID"),logical(1))),
      nonconverged=if(is.null(fit)) 0L else sum(vapply(fit$runs,function(r) identical(r$status,"NONCONVERGED_VALID"),logical(1))),
      failed=if(is.null(fit)) NA_integer_ else sum(vapply(fit$runs,function(r) !r$status %in% c("CONVERGED_VALID","NONCONVERGED_VALID"),logical(1))),
      training_objective=if(is.null(fit$representative)) NA_real_ else fit$representative$evaluation$objective,
      mean_nb2_nll=if(is.null(x$forecast)) NA_real_ else x$forecast$scores$summary$mean_nb2_nll,
      reason=if(is.null(x$reason)) "" else x$reason)
  }))
}
p3_choose <- function(table,rule) {
  scalar_number(rule$tie_tolerance,"tie tolerance",0)
  if(!identical(rule$rule,"smaller_K_then_config_order")) stop("Unknown selection rule")
  eligible <- which(table$status=="AVAILABLE" & is.finite(table$mean_nb2_nll))
  if(!length(eligible)) return(list(status="SELECTION_FAILED_ALL_UNAVAILABLE",id=NULL,rule=rule))
  best <- min(table$mean_nb2_nll[eligible])
  tied <- eligible[table$mean_nb2_nll[eligible]<=best+rule$tie_tolerance]
  selected <- tied[order(table$K[tied],tied)][1]
  list(status="SELECTED_ON_B",id=table$id[selected],candidate=table[selected,c("id","K","lambda"),drop=FALSE],
    best_B_score=best,tied_ids=table$id[tied],rule=rule,
    interpretation="Absolute NLL tolerance, smaller K then original config order; lambda is not ordered across K")
}
