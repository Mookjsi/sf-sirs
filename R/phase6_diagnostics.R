# Explicit diagnostic action: profiles refit nuisance parameters, reports never do.
p6_physical <- function(q,map) {
  p <- decode_parameters(q,map)
  c(beta0=p$beta0,phi=p$phi,S0=p$S0,I0=p$I0,rho=p$rho,omega=p$omega,delta=p$delta,
    mu=p$mu,gamma=p$gamma,setNames(p$free_coordinates,map$shape_order))
}
p6_hessian <- function(fit,problem,step,condition_limit,component="likelihood") {
  if(!component %in% c("likelihood","penalized")) stop("Unknown Hessian component")
  q <- fit$coordinates;d <- length(q);h <- rep(step,d)
  tryCatch({
    if(length(fit$boundary$coordinates) || any(q-2*h<=problem$map$lower | q+2*h>=problem$map$upper)) stop("Boundary or stencil outside guards")
    # Probe precisely the stencil used by the central mixed/diagonal differences.
    fn <- function(v) {e <- evaluate_candidate(v,problem);if(!e$valid) stop("Invalid Hessian stencil");if(component=="likelihood") -e$loglik else e$objective}
    center <- fn(q);H <- matrix(0,d,d,dimnames=list(names(q),names(q)))
    eye <- diag(h,d)
    for(i in seq_len(d)) {
      H[i,i] <- (fn(q+eye[,i])-2*center+fn(q-eye[,i]))/h[i]^2
      if(i>1L) for(j in seq_len(i-1L)) H[i,j] <- H[j,i] <-
        (fn(q+eye[,i]+eye[,j])-fn(q+eye[,i]-eye[,j])-fn(q-eye[,i]+eye[,j])+fn(q-eye[,i]-eye[,j]))/(4*h[i]*h[j])
    }
    eigenvalues <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
    condition <- if(min(eigenvalues)>0) max(eigenvalues)/min(eigenvalues) else Inf
    usable <- is.finite(condition) && condition<condition_limit
    covariance <- correlation <- uncertainty <- NULL
    if(usable) {
      V <- solve(H);base <- p6_physical(q,problem$map)
      J <- vapply(seq_len(d),function(i) (p6_physical(q+eye[,i],problem$map)-p6_physical(q-eye[,i],problem$map))/(2*h[i]),numeric(length(base)))
      rownames(J) <- names(base);covariance <- J%*%V%*%t(J)
      varying <- diag(covariance)>0
      for(n in c("mu","gamma")) if(!problem$map$joint[[n]]) varying[names(base)==n] <- FALSE
      correlation <- matrix(NA_real_,length(base),length(base),dimnames=list(names(base),names(base)))
      correlation[varying,varying] <- stats::cov2cor(covariance[varying,varying,drop=FALSE])
      se <- sqrt(pmax(0,diag(covariance)))
      uncertainty <- data.frame(parameter=names(base),estimate=base,local_SE=ifelse(varying,se,NA_real_),
        definition="Local delta-method approximation; fixed values have NA; no coverage certification")
    }
    # Numerical gradient is recorded separately: positive curvature need not be stationary.
    gradient <- vapply(seq_len(d),function(i) (fn(q+eye[,i])-fn(q-eye[,i]))/(2*h[i]),numeric(1))
    geometry_covariance <- NULL
    if(!isTRUE(problem$map$fourier)) {
      geometry_covariance <- covariance;covariance <- NULL;uncertainty <- NULL
    }
    answer <- list(status=if(usable) "LOCAL_CURVATURE_AVAILABLE_NOT_IDENTIFIABILITY_PROOF" else "NONPOSITIVE_OR_ILL_CONDITIONED_NO_INVERSE",
      step=step,hessian=H,eigenvalues=eigenvalues,condition=condition,gradient=gradient,
      covariance=covariance,correlation=correlation,uncertainty=uncertainty)
    if(!isTRUE(problem$map$fourier)) {
      answer$component <- component;answer$geometry_covariance <- geometry_covariance
      answer$interpretation <- "Inverse curvature geometry at a penalized endpoint, NOT ordinary frequentist covariance; likelihood gradient need not vanish"
    }
    answer
  },error=function(e) list(status="UNAVAILABLE",step=step,reason=conditionMessage(e)))
}
p6_profile <- function(cfg,data,name,fit,target,optimizer=stats::optim,criterion=if(name=="SPL") "penalized" else "likelihood") {
  if(!criterion %in% c("likelihood","penalized")) stop("Unknown profile criterion")
  rows <- list();records <- list();values <- sort(unique(c(cfg$diagnostics$profile[[target]],fit$evaluation$parameters[[target]])))
  for(i in seq_along(values)) {
    v <- values[i]
    records[[i]] <- tryCatch({
      model <- p6_model(cfg,name);model[[target]] <- v;model$fixed_sources[[target]] <- "Explicit profile likelihood constraint"
      if(name=="SPL" && criterion=="likelihood") model$lambda <- 0
      estimation <- cfg$estimation;estimation$free[[target]] <- FALSE
      problem <- prepare_estimation_problem(data,model,estimation)
      start <- fit$evaluation$parameters
      if(target=="S0") start$I0 <- (1-v)*start$I0/(1-start$S0)
      if(target=="I0") start$S0 <- (1-v)*start$S0/(1-start$I0)
      start[[target]] <- v;start$profile_target <- target
      local <- cfg;local$starts$user <- list(A=NULL,AS=NULL,SPL=NULL)
      starts <- p6_starts(local,problem,name,start,cfg$diagnostics$profile_starts)
      f <- fit_multistart(data,model,estimation,p6_start_config(cfg),start,optimizer,explicit_starts=starts)
      r <- f$representative
      list(fit=f,status=if(is.null(r)) "UNAVAILABLE" else "NUISANCE_REOPTIMIZED",loglik=if(is.null(r)) NA_real_ else r$evaluation$loglik,
        boundary=if(is.null(r)) "" else paste(r$boundary$coordinates,collapse=";"),
        objective=if(is.null(r)) NA_real_ else r$evaluation$objective,penalty=if(is.null(r)) NA_real_ else r$evaluation$penalty)
    },error=function(e) list(status="FAILED",loglik=NA_real_,reason=conditionMessage(e),boundary=""))
    r <- records[[i]]
    rows[[i]] <- data.frame(model=name,parameter=target,value=v,status=r$status,logLik=r$loglik,
      deviance=2*(fit$evaluation$loglik-r$loglik),better_than_saved_MLE=is.finite(r$loglik) && r$loglik>fit$evaluation$loglik+1e-6,
      boundary=r$boundary,reason=if(is.null(r$reason)) "" else r$reason)
    if(name=="SPL") {
      rows[[i]]$criterion <- criterion
      rows[[i]]$penalty <- if(is.null(r$penalty)) NA_real_ else r$penalty
      rows[[i]]$objective <- if(is.null(r$objective)) NA_real_ else r$objective
      baseline <- if(criterion=="likelihood") -fit$evaluation$loglik else fit$evaluation$objective
      rows[[i]]$objective_contrast <- 2*(rows[[i]]$objective-baseline)
      rows[[i]]$reference <- "Saved penalized endpoint; contrasts are NOT calibrated LR or confidence intervals"
      rows[[i]]$better_than_saved_endpoint <- rows[[i]]$better_than_saved_MLE
    }
  }
  list(table=do.call(rbind,rows),records=records,
    interpretation="Finite-grid nuisance reoptimization, not a certified confidence interval; failed profiles retained")
}
p6_diagnostic_targets <- function(cfg) {
  targets <- if(is.null(cfg$diagnostics$targets)) c("S0","I0","rho") else cfg$diagnostics$targets
  if(!is.character(targets) || !length(targets) || anyDuplicated(targets) || any(!targets %in% c("S0","I0","rho"))) stop("Select profile targets S0/I0/rho")
  targets
}
p6_diagnostic_choices <- function(cfg,mode=NULL) {
  if(is.null(mode)) mode <- if(is.null(cfg$diagnostics$mode)) "full" else cfg$diagnostics$mode
  if(length(mode)!=1L || is.na(mode) || !mode %in% c("full","hessian-only","profile-only")) stop("Unknown diagnostic mode")
  if(!is.null(cfg$diagnostics$spline_profile) && !cfg$diagnostics$spline_profile %in% c("both","likelihood","penalized")) stop("Unknown spline profile criterion")
  mode
}
p6_diagnose <- function(result,optimizer=stats::optim,mode=NULL) {
  if(!identical(result$type,"phase6_fit")) stop("Joint MLE fit required")
  cfg <- result$config;records <- list();mode <- p6_diagnostic_choices(cfg,mode)
  for(name in names(result$models)) {
    fit <- result$models[[name]]$fit$representative
    if(is.null(fit)) {records[[name]] <- list(status="FIT_UNAVAILABLE");next}
    problem <- p6_problem(cfg,result$data,name)
    hessians <- profiles <- list()
    if(mode!="profile-only") {
      hessians <- lapply(cfg$diagnostics$hessian_steps,function(step) p6_hessian(fit,problem,step,cfg$diagnostics$condition_limit))
      if(name=="SPL") {
        penalized <- lapply(cfg$diagnostics$hessian_steps,function(step) p6_hessian(fit,problem,step,cfg$diagnostics$condition_limit,"penalized"))
        for(i in seq_along(hessians)) hessians[[i]]$component <- "likelihood"
        for(i in seq_along(penalized)) penalized[[i]]$component <- "penalized"
        hessians <- c(hessians,penalized)
      }
    }
    if(mode!="hessian-only") for(n in p6_diagnostic_targets(cfg)) {
      if(name!="SPL") profiles[[n]] <- p6_profile(cfg,result$data,name,fit,n,optimizer) else {
        choice <- cfg$diagnostics$spline_profile
        criteria <- if(is.null(choice) || choice=="both") c("likelihood","penalized") else choice
        for(criterion in criteria) profiles[[paste(n,criterion,sep="_")]] <- p6_profile(cfg,result$data,name,fit,n,optimizer,criterion)
      }
    }
    runs <- result$models[[name]]$fit$runs
    good <- Filter(function(r) identical(r$status,"CONVERGED_VALID") && r$evaluation$objective<=fit$evaluation$objective+cfg$diagnostics$near_optimal_nll,runs)
    endpoints <- if(length(good)) do.call(rbind,lapply(good,function(r) p6_physical(r$coordinates,r$map))) else NULL
    endpoint_cor <- if(length(good)>=3L) suppressWarnings(stats::cor(endpoints)) else NULL
    records[[name]] <- list(status="DIAGNOSTICS_STORED_NOT_CERTIFIED",hessians=hessians,profiles=profiles,
      near_optimal_endpoints=endpoints,endpoint_correlation=endpoint_cor,
      warning="Between-start correlations are optimization geometry, not sampling uncertainty. Inspect step stability and profile branches.")
  }
  list(type="phase6_diagnostics",status="DIAGNOSTICS_STORED",records=records,mode=mode,
    config=result$config,identity=result$identity,source_fit_hash=result$artifact_hash)
}
