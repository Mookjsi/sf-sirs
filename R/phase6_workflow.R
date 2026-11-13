# Independent joint MLE using the existing core, chart, NB2 and optimizer.
p6_models <- function(cfg) {
  models <- if(is.null(cfg$models)) c("A","AS") else cfg$models
  if(!is.character(models) || !length(models) || anyDuplicated(models) ||
     any(!models %in% c("A","AS","SPL"))) stop("Select unique A/AS/SPL models")
  if("AS" %in% models && "A" %in% models && match("AS",models)<match("A",models)) stop("A must precede AS for nested warm start")
  models
}
p6_model <- function(cfg,name) {
  if(name=="SPL") {
    K <- if(is.null(cfg$spline$K)) 8L else cfg$spline$K
    lambda <- if(is.null(cfg$spline$lambda)) 1 else cfg$spline$lambda
    integer_scalar(K,"SPL K",4L);scalar_number(lambda,"SPL lambda",0)
    return(list(K=K,lambda=lambda,transmission="spline",fixed_sources=cfg$fixed_sources))
  }
  if(!name %in% c("A","AS")) stop("Phase 6 model must be A, AS or SPL")
  list(K=8L,lambda=0,transmission="fourier",harmonics=if(name=="A") 1L else 2L,
    fixed_sources=cfg$fixed_sources)
}
p6_data <- function(cfg) {
  d <- p3_load_data(cfg)
  d$fixed_parameters$mu <- cfg$fixed$mu;d$fixed_parameters$gamma <- cfg$fixed$gamma
  d
}
p6_problem <- function(cfg,data,name) prepare_estimation_problem(data,p6_model(cfg,name),cfg$estimation)
p6_start_config <- function(cfg) list(seed=cfg$starts$seed,beta0_candidates=cfg$initial$beta0,
  phi_candidates=cfg$initial$phi,n_perturbed=0L,perturb_sd=cfg$starts$sd,user_starts=NULL)
p6_starts <- function(cfg,problem,name,base=NULL,n=cfg$starts$n) {
  if(is.null(base)) {base <- cfg$initial;base$free_coordinates <- rep(0,problem$basis$K-1L)}
  q <- encode_parameters(base,problem$map)
  with_phase2_seed(cfg$starts$seed,function() {
    starts <- list(list(label="independent-design",parameters=base))
    attempts <- 0L
    while(length(starts)<n) {
      attempts <- attempts+1L;if(attempts>10000L) stop("Start rejection limit: review guards/design")
      v <- q+stats::rnorm(length(q),sd=cfg$starts$sd)
      if(any(v<=problem$map$lower | v>=problem$map$upper)) next
      starts[[length(starts)+1L]] <- list(label=paste0("seeded-",length(starts)),parameters=decode_parameters(v,problem$map))
    }
    if(is.null(base$profile_target) && length(cfg$starts$user[[name]])) for(i in seq_along(cfg$starts$user[[name]]))
      starts[[length(starts)+1L]] <- list(label=paste0("user-",i),parameters=cfg$starts$user[[name]][[i]])
    for(s in starts) encode_parameters(s$parameters,problem$map)
    starts
  })
}
p6_check <- function(cfg,data) {
  p3_text(cfg$campaign_id,"campaign_id")
  if(!grepl("^[A-Za-z0-9_-]+$",cfg$campaign_id) || !identical(cfg$independent,TRUE)) stop("Declare independent design and safe campaign ID")
  if(!identical(cfg$estimation$free,list(S0=TRUE,I0=TRUE,rho=TRUE)) ||
    !isTRUE(cfg$estimation$joint$omega) || !isTRUE(cfg$estimation$joint$delta)) stop("Joint S0/I0/rho/omega/delta must be free")
  validate_phase2_data(data)
  if(!identical(cfg$data$include_npi,FALSE) || !is.null(cfg$data$post_npi) || any(data$npieff!=1)) stop("Phase 6 no-NPI contract violated")
  # This path guard supplements (and cannot replace) scientific provenance review.
  paths <- p3_inputs(cfg)
  if(any(grepl("posterior|/fits/|/extract/|/artifacts/",paths,ignore.case=TRUE))) stop("Fitted-value input paths forbidden in independent Phase 6")
  w <- cfg$expected_window
  if(length(data$cases)!=w$n || data$definition$region!=w$region || data$processed_week_marker[1]!=w$first ||
    tail(data$processed_week_marker,1)!=w$last) stop("Unexpected independent estimation data window")
  for(n in c("mu","gamma")) {
    scalar_number(cfg$fixed[[n]],n,0,strict_lower=TRUE);p3_text(cfg$fixed_sources[[n]],paste(n,"source"))
    if(data$fixed_parameters[[n]]!=cfg$fixed[[n]]) stop("Fixed-rate data/config mismatch")
  }
  integer_scalar(cfg$starts$n,"starts",2L);integer_scalar(cfg$starts$seed,"start seed",0L)
  scalar_number(cfg$starts$sd,"start sd",0,strict_lower=TRUE)
  p6_validate_priors(cfg)
  for(n in c("seed","chains","parallel_chains","iter_warmup","iter_sampling","max_treedepth","predictive_draws"))
    integer_scalar(cfg$bayes[[n]],paste("bayes",n),1L)
  if(cfg$bayes$parallel_chains>cfg$bayes$chains) stop("parallel_chains exceeds chains")
  scalar_number(cfg$bayes$adapt_delta,"adapt_delta",0,1,strict_lower=TRUE)
  if(cfg$bayes$adapt_delta>=1) stop("adapt_delta must be below one")
  integer_scalar(cfg$prior_predictive$n,"prior draws",1L)
  integer_scalar(cfg$prior_predictive$seed,"prior seed",0L)
  scalar_number(cfg$diagnostics$condition_limit,"condition limit",1,strict_lower=TRUE)
  scalar_number(cfg$diagnostics$near_optimal_nll,"near optimal NLL threshold",0)
  for(n in c("mesh","derivative_tolerance","root_tolerance")) scalar_number(cfg$peak_detection[[n]],n,0,strict_lower=TRUE)
  design <- lapply(p6_models(cfg),function(name) {
    p <- p6_problem(cfg,data,name);s <- p6_starts(cfg,p,name)
    data.frame(model=name,dimension=p$map$dimension,starts=length(s))
  })
  p6_diagnostic_choices(cfg)
  for(n in p6_diagnostic_targets(cfg)) {
    v <- cfg$diagnostics$profile[[n]];finite_numeric(v,paste("profile",n))
    if(any(v<=0 | v>=1) || any(diff(v)<=0)) stop("Profiles require increasing interior fractions")
  }
  finite_numeric(cfg$diagnostics$hessian_steps,"hessian steps")
  if(any(cfg$diagnostics$hessian_steps<=0)) stop("Positive Hessian steps required")
  integer_scalar(cfg$diagnostics$profile_starts,"profile starts",2L)
  list(type="phase6_check",status="CONFIG_CHECKED_NO_TRAJECTORY_OR_FIT",start_design=do.call(rbind,design),
    data_audit=p3_data_audit(data,cfg),warning=paste(cfg$limitations,collapse="; "))
}
p6_fit <- function(cfg,data,optimizer=stats::optim) {
  records <- list()
  for(name in p6_models(cfg)) records[[name]] <- tryCatch({
    problem <- p6_problem(cfg,data,name);starts <- p6_starts(cfg,problem,name)
    warm_status <- "NOT_APPLICABLE";difference <- NA_real_
    annual <- records$A$fit$representative
    if(name=="AS" && !is.null(annual)) {
      warm <- annual$evaluation$parameters;warm$free_coordinates <- c(warm$free_coordinates,0,0)
      e <- evaluate_candidate(encode_parameters(warm,problem$map),problem)
      difference <- e$loglik-annual$evaluation$loglik
      if(!e$valid || !is.finite(difference) || abs(difference)>1e-8*max(1,abs(e$loglik))) stop("A-in-AS likelihood identity failed")
      starts[[length(starts)+1L]] <- list(label="nested-annual-zero-semiannual",parameters=warm)
      warm_status <- "ANNUAL_ENDPOINT_EMBEDDED"
    }
    fit <- fit_multistart(data,problem$model,cfg$estimation,p6_start_config(cfg),cfg$initial,
      optimizer=optimizer,explicit_starts=starts)
    list(status=if(is.null(fit$representative)) fit$status else "AVAILABLE",fit=fit,
      warm_start=warm_status,warm_loglik_difference=difference,reason=if(is.null(fit$representative)) "No converged valid endpoint" else "")
  },error=function(e) list(status="FAILED",reason=conditionMessage(e),fit=NULL,warm_start="UNAVAILABLE"))
  table <- p5_comparison(records)
  table$AIC_BIC_definition <- "Joint unpenalized local maximum likelihood; all free nuisance coordinates counted; global optimum unverified"
  if("SPL" %in% table$model) {
    table$AIC_BIC_definition[table$model=="SPL"] <- "NOT_COMPARABLE: penalized spline effective degrees of freedom unavailable; in-sample likelihood descriptive only"
    table$penalized_objective <- vapply(records,function(r) if(is.null(r$fit$representative)) NA_real_ else r$fit$representative$evaluation$objective,numeric(1))
    table$roughness <- vapply(records,function(r) {
      f <- r$fit$representative
      if(is.null(f) || isTRUE(f$map$fourier)) return(if(is.null(f)) NA_real_ else 0)
      cyclic_penalty(f$evaluation$spline_coefficients,1)
    },numeric(1))
  }
  lr <- p5_lr(table);lr$nominal_p_value <- lr$p_value;lr$p_value <- NA_real_
  lr$definition <- "Exploratory chi-square(2) reference only: regularity, identifiability, interiority and global optimization unverified; no calibrated p-value"
  list(type="phase6_fit",status=if(all(table$status=="AVAILABLE")) "LOCAL_ENDPOINTS_AVAILABLE" else "FAILURES_RETAINED",
    config=cfg,models=records,comparison=table,lr=lr,data=data,markers=p3_markers(data),observed=data$cases,
    mask=data$likelihood_include,data_audit=p3_data_audit(data,cfg))
}
p6_simulate <- function(cfg,data,name,parameters,seed) {
  problem <- p6_problem(cfg,data,name)
  e <- p6_propagate(encode_parameters(parameters,problem$map),problem)
  y <- with_phase2_seed(seed,function() stats::rnbinom(length(data$cases),size=parameters$phi,mu=e$states$expected_cases))
  if(any(!is.finite(y))) stop("Synthetic NB2 draw failed")
  synthetic <- data;synthetic$cases <- y;synthetic$time_markers$cases <- y
  synthetic$provenance <- list(source="Explicit Phase 6 synthetic draw")
  list(data=synthetic,truth=parameters,states=e$states,seed=seed)
}
p6_propagate <- function(q,problem) {
  p <- decode_parameters(q,problem$map)
  profile <- if(isTRUE(problem$map$fourier)) fourier_profile(p$free_coordinates,p$beta0,problem$basis) else
    spline_profile(p$free_coordinates,p$beta0,problem$basis)
  fixed <- list(pop=problem$data$fixed_parameters$pop,mu=p$mu,gamma=p$gamma,
    omega=p$omega,delta=p$delta,rho=p$rho)
  states <- simulate_weekly_sirs(profile$beta_by_week,make_initial_state(p$S0,p$I0,fixed$pop),
    fixed,problem$data$seasonal_week,problem$data$npieff)
  list(parameters=p,states=states)
}
p6_recover <- function(cfg,data,optimizer=stats::optim) {
  if(is.null(cfg$recovery$truth)) stop("Set explicit independent recovery.truth for every selected model first")
  integer_scalar(cfg$recovery$replicates,"recovery replicates",1L)
  records <- list()
  for(name in p6_models(cfg)) for(i in seq_len(cfg$recovery$replicates)) {
    seed <- cfg$recovery$seed+i+if(name=="AS") 100000L else if(name=="SPL") 200000L else 0L
    records[[paste(name,i,sep="-")]] <- tryCatch({
      s <- p6_simulate(cfg,data,name,cfg$recovery$truth[[name]],seed)
      f <- p6_fit(cfg,s$data,optimizer)
      r <- f$models[[name]]$fit$representative
      truth <- cfg$recovery$truth[[name]]
      list(model=name,replicate=i,seed=seed,truth=truth,synthetic=s,fit=f,
        errors=if(is.null(r)) NULL else p6_physical(r$coordinates,r$map)-p6_physical(encode_parameters(truth,r$map),r$map),
        status=if(is.null(r)) "FIT_UNAVAILABLE" else "RECOVERY_ERRORS_STORED_NOT_CERTIFIED")
    },error=function(e) list(model=name,replicate=i,seed=seed,status="FAILED",reason=conditionMessage(e)))
  }
  list(type="phase6_recovery",status="SYNTHETIC_EXPERIMENT_RECORDS_STORED",records=records,
    interpretation="Inspect bias, dispersion, failures, boundary hits and all starts; fitting success is not successful recovery")
}
p6_identity <- function(root,cfg,path) {
  files <- unique(c(list.files(file.path(root,"R"),"[.]R$",full.names=TRUE),
    file.path(root,c("scripts/phase6.R","scripts/phase6.sh","stan/phase6_fourier_sirs.stan","config/phase2_defaults.R")),
    path,cfg$config_dependencies,p3_inputs(cfg)))
  tools::md5sum(normalizePath(files,mustWork=TRUE))
}
