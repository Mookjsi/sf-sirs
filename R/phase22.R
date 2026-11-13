# Phase 22 uses the production parameter chart, recurrence, NB2 likelihood and penalties.
p22_progress <- function(message,...) {
  callback <- getOption("phase22.progress")
  if(is.function(callback)) callback(sprintf(message,...))
  invisible(NULL)
}
p22_inputs <- function(cfg) {
  paths <- unlist(cfg[c("reference","replay","panel","plan","fit","evaluation")])
  if(any(!file.exists(paths))) stop("Missing frozen input: ",paste(paths[!file.exists(paths)],collapse=", "))
  saved <- lapply(paths,read_phase2_result)
  names(saved) <- names(paths)
  ref <- ae_reference(saved$reference,cfg$reference)
  replay <- saved$replay$result
  if(!identical(replay$status,"HISTORICAL_REPLAY_PASSED") ||
     !identical(replay$source_hash,ref$artifact_hash) ||
     !all(replay$historical_replay_checks$passed)) stop("Historical replay gate failed")
  live <- ae_replay(ref,cfg$replay_tolerance)
  if(!identical(live$status,"HISTORICAL_REPLAY_PASSED")) stop("Current historical replay failed")
  plan <- saved$plan$result
  bounds <- plan$settings$q_log_bounds
  if(length(bounds)!=2L || any(!is.finite(bounds)) ||
     !identical(as.numeric(bounds),as.numeric(log(c(.05,20)))))
    stop("Saved Phase 21 q bounds differ from accepted extension")
  registry <- plan$registry
  row <- which(registry$id==cfg$selected)
  if(length(row)!=1L || registry$family[row]!="E2-SPLINE" ||
     registry$spacing[row]!=8L || registry$lambda[row]!=.1) stop("Selected candidate configuration mismatch")
  basis <- plan$bases[[cfg$selected]]
  if(is.null(basis) || basis$dimension!=29L || nrow(basis$X)!=209L)
    stop("Saved K8 basis does not mean 8-week knot spacing with 29 coordinates")
  frozen <- saved$fit$result$fits[[cfg$selected]]
  if(is.null(frozen$theta) || !identical(frozen$status,"NUMERICALLY_VALID"))
    stop("Selected frozen correction unavailable")
  old_env <- new.env(parent=baseenv())
  sys.source("config/accepted_extension.R",old_env)
  calendar <- ae_calendar(saved$panel$result$island_cases,old_env$accepted_extension_config())
  if(!identical(calendar,plan$calendar_bridge)) stop("Phase 21 calendar mismatch")
  list(ref=ref,replay=live,calendar=calendar,basis=basis,
       candidate=p21_candidate(registry,row),frozen=frozen,
       q_log_bounds=bounds,
       hashes=tools::md5sum(paths),saved=saved)
}

p22_problem <- function(inputs,family="SPL",free_biology=FALSE,
                        fixed_nuisance=FALSE) {
  ref <- inputs$ref; cfg <- ref$source$config
  if(free_biology) {
    cfg$estimation$joint$omega <- TRUE
    cfg$estimation$joint$delta <- TRUE
  } else {
    cfg$estimation$joint$omega <- FALSE
    cfg$estimation$joint$delta <- FALSE
  }
  cfg$estimation$joint$mu <- FALSE
  cfg$estimation$joint$gamma <- FALSE
  cfg$fixed$omega <- ref$parameters$omega
  cfg$fixed$delta <- ref$parameters$delta
  if(fixed_nuisance) cfg$estimation$free <- list(S0=FALSE,I0=FALSE,rho=FALSE)
  cfg$spline$lambda <- 1
  data <- inputs$replay$data
  model <- p6_model(cfg,family)
  model$omega <- ref$parameters$omega
  model$delta <- ref$parameters$delta
  model$fixed_sources$omega <- "Accepted Phase 7C lambda=1 endpoint"
  model$fixed_sources$delta <- "Accepted Phase 7C lambda=1 endpoint"
  if(fixed_nuisance) {
    for(k in c("S0","I0","rho")) {
      model[[k]] <- ref$parameters[[k]]
      model$fixed_sources[[k]] <- "Accepted Phase 7C lambda=1 endpoint"
    }
  }
  problem <- prepare_estimation_problem(data,model,cfg$estimation)
  list(problem=problem,cfg=cfg)
}

p22_score <- function(coords,post,inputs,setup,mask=NULL) {
  problem <- setup$problem; basis <- inputs$basis; calendar <- inputs$calendar
  n <- length(problem$data$cases); nf <- nrow(calendar)
  if(is.null(mask)) mask <- c(FALSE,rep(TRUE,n-1L+nf))
  if(length(mask)!=n+nf || anyNA(mask) || !is.logical(mask) || mask[1L]) stop("Invalid common likelihood mask")
  if(!p21_feasible(inputs$candidate,basis,post,
    list(q_log_bounds=inputs$q_log_bounds))) stop("Post correction outside saved continuous bounds")
  p <- decode_parameters(coords,problem$map)
  profile <- if(isTRUE(problem$map$fourier))
    fourier_profile(p$free_coordinates,p$beta0,problem$basis) else
    spline_profile(p$free_coordinates,p$beta0,problem$basis)
  g <- as.vector(basis$X%*%post)
  beta <- profile$beta_by_week
  pop <- problem$data$fixed_parameters$pop
  fixed <- list(pop=pop,mu=ref_or(p$mu,problem$data$fixed_parameters$mu),
    gamma=ref_or(p$gamma,problem$data$fixed_parameters$gamma),
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  weeks <- c(problem$data$seasonal_week,calendar$seasonal_week)
  transition <- c(rep(1,n),exp(g))
  states <- simulate_weekly_sirs(beta,make_initial_state(p$S0,p$I0,pop),fixed,
    weeks,transition)
  observed <- c(problem$data$cases,calendar$observed)
  likelihood <- nb2_loglik(observed,states$expected_cases,p$phi,mask)
  seasonal_penalty <- if(isTRUE(problem$map$fourier)) 0 else
    cyclic_penalty(profile$coefficients,problem$model$lambda)
  post_penalty <- p21_penalty(inputs$candidate,basis,post)
  list(objective=-likelihood$total+seasonal_penalty+post_penalty,
    loglik=likelihood$total,historical_loglik=sum(likelihood$included_per_week[seq_len(n)]),
    extension_loglik=sum(likelihood$included_per_week[n+seq_len(nf)]),
    seasonal_penalty=seasonal_penalty,post_penalty=post_penalty,
    states=states,parameters=p,beta=beta,q=exp(g),
    per_week=likelihood$per_week,mask=mask,observed=observed)
}

ref_or <- function(x,y) if(is.null(x)) y else x

p22_fit <- function(inputs,setup,cfg,mask=NULL,seed=cfg$seed,
                    allow_full_start=TRUE,extra_starts=list()) {
  problem <- setup$problem; basis <- inputs$basis
  if(isTRUE(problem$map$fourier)) {
    w <- 0:51;h <- problem$basis$harmonics
    X <- do.call(cbind,lapply(seq_len(h),function(j)
      cbind(cos(2*pi*j*w/52),sin(2*pi*j*w/52))))
    fitted <- inputs$ref$evaluation$parameters
    fitted$free_coordinates <- as.vector(qr.solve(X,log(inputs$ref$beta)))
    historical <- encode_parameters(fitted,problem$map)
  } else historical <- inputs$ref$record$retained_coordinates[problem$map$order]
  if(anyNA(historical)) stop("Historical coordinate chart mismatch")
  baseline <- c(historical,if(allow_full_start) inputs$frozen$theta else
    rep(0,inputs$basis$dimension))
  d <- length(historical)
  controls <- cbind(matrix(0,nrow(basis$control),d),basis$control)
  # Keep the full continuous spline feasible using saved Bernstein controls.
  ui <- rbind(cbind(diag(d),matrix(0,d,basis$dimension)),
              cbind(-diag(d),matrix(0,d,basis$dimension)),controls,-controls)
  ci <- c(problem$map$lower,-problem$map$upper,
          rep(inputs$q_log_bounds[1],nrow(controls)),
          rep(-inputs$q_log_bounds[2],nrow(controls)))
  evaluations <- 0L;best_objective <- Inf;started <- proc.time()[[3L]]
  fn <- function(x) {
    evaluations <<- evaluations+1L
    z <- tryCatch(p22_score(x[seq_len(d)],x[-seq_len(d)],inputs,setup,mask),error=function(e) NULL)
    value <- if(is.null(z) || !is.finite(z$objective)) 1e20 else z$objective
    if(is.finite(value)) best_objective <<- min(best_objective,value)
    if(evaluations%%500L==0L)
      p22_progress("fit objective evaluations=%d; best=%0.5f; elapsed=%0.0fs",
        evaluations,best_objective,proc.time()[[3L]]-started)
    value
  }
  set.seed(seed)
  starts <- vector("list",cfg$starts)
  starts[[1L]] <- baseline
  design <- setup$cfg$initial
  design$omega <- inputs$ref$parameters$omega
  design$delta <- inputs$ref$parameters$delta
  if(!problem$map$free$S0) design$S0 <- inputs$ref$parameters$S0
  if(!problem$map$free$I0) design$I0 <- inputs$ref$parameters$I0
  if(!problem$map$free$rho) design$rho <- inputs$ref$parameters$rho
  design$free_coordinates <- rep(0,problem$basis$K-1L)
  independent <- encode_parameters(design,problem$map)
  if(cfg$starts>=2L) for(i in 2:cfg$starts) {
    z <- baseline
    z[seq_len(d)] <- pmin(problem$map$upper-1e-5,pmax(problem$map$lower+1e-5,
      if(i==2L) independent else
          z[seq_len(d)]+rnorm(d,sd=cfg$perturb_sd)))
    z[-seq_len(d)] <- 0
    starts[[i]] <- z
  }
  if(length(extra_starts)) {
    if(any(!vapply(extra_starts,function(z) is.numeric(z) &&
      length(z)==ncol(ui) && all(is.finite(z)),logical(1))))
      stop("Invalid supplied optimizer start")
    starts <- c(starts,extra_starts)
  }
  attempts <- lapply(seq_along(starts),function(i) {
    x <- starts[[i]]
    p22_progress("fit start %d/%d: %d seasonal/biological and %d post coordinates",
      i,length(starts),d,basis$dimension)
    if(any(as.vector(ui%*%x-ci)<=0))
      return(list(start=i,status="INFEASIBLE_START",coordinates=x))
    grad <- function(y) p21_numeric_gradient(y,fn,
      function(z) all(as.vector(ui%*%z-ci)>0),cfg$gradient_step)
    opt <- tryCatch(stats::constrOptim(x,fn,grad=grad,ui=ui,ci=ci,
      control=list(maxit=cfg$maxit,reltol=1e-8),outer.iterations=20L),error=function(e)e)
    if(inherits(opt,"error")) {
      p22_progress("fit start %d/%d optimizer error: %s",i,length(starts),conditionMessage(opt))
      return(list(start=i,status="OPTIMIZER_ERROR",message=conditionMessage(opt)))
    }
    score <- tryCatch(p22_score(opt$par[seq_len(d)],opt$par[-seq_len(d)],inputs,setup,mask),error=function(e)e)
    valid <- !inherits(score,"error") && is.finite(score$objective)
    p22_progress("fit start %d/%d ended: convergence=%s objective=%s",
      i,length(starts),as.character(opt$convergence),
      if(valid) format(score$objective,digits=10) else "invalid")
    list(start=i,status=if(valid && opt$convergence==0L) "NUMERICALLY_VALID" else
      if(valid) "NONCONVERGED_ENDPOINT" else "INVALID_ENDPOINT",
      initial=x,coordinates=opt$par,
      convergence=opt$convergence,outer_iterations=opt$outer.iterations,
      barrier_value=opt$barrier.value,
      message=opt$message,reported=opt$value,
      replay_difference=if(valid) score$objective-opt$value else NA_real_,
      score=if(valid) score else NULL)
  })
  ok <- which(vapply(attempts,function(z) identical(z$status,"NUMERICALLY_VALID"),logical(1)))
  if(!length(ok)) return(list(status="FAILED",attempts=attempts))
  best <- ok[which.min(vapply(attempts[ok],function(z) z$score$objective,numeric(1)))]
  list(status="ENDPOINTS_RECORDED",attempts=attempts,best=best,
       competitive=ok[vapply(attempts[ok],function(z)
         z$score$objective-attempts[[best]]$score$objective<=2,logical(1))],
       setup=setup,seed=seed)
}
