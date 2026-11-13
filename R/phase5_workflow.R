# Phase 5 orchestration reuses Phase 2 fitting and Phase 3 data/scoring/forecast code.
p5_model <- function(cfg,name) {
  if(!name %in% c("A","AS","SPL")) stop("Unknown Phase 5 model")
  model <- cfg$model
  model$fixed_sources <- lapply(cfg$fixed_provenance,function(x) paste(x$category,x$citation,sep=": "))
  model$K <- 8L;model$lambda <- if(name=="SPL") 1 else 0
  model$transmission <- if(name=="SPL") "spline" else "fourier"
  if(name!="SPL") model$harmonics <- if(name=="A") 1L else 2L
  model
}
p5_start_config <- function(cfg,name) {
  starts <- cfg$starts
  if(name!="SPL") {
    starts$seed <- cfg$fourier$seed;starts$perturb_sd <- cfg$fourier$coefficient_scale
    size <- length(unique(starts$beta0_candidates))*length(unique(starts$phi_candidates))
    starts$n_perturbed <- as.integer(max(1,ceiling(cfg$fourier$n_starts/size)-1))
    starts$user_starts <- NULL
  }
  starts
}
p5_starts <- function(cfg,problem,name) {
  initial <- cfg$initial_parameters;initial$free_coordinates <- rep(0,problem$basis$K-1L)
  if(name!="SPL" && !is.null(cfg$fourier$user_starts[[name]])) {
    rows <- cfg$fourier$user_starts[[name]]
    if(!is.list(rows) || length(rows)<8L) stop("Fourier user starts must provide at least eight physical parameter lists")
    starts <- lapply(seq_along(rows),function(i) list(label=paste0("user-",i),parameters=rows[[i]]))
  } else {
    starts <- generate_starts(initial,problem,p5_start_config(cfg,name))
    if(name!="SPL") starts <- starts[seq_len(cfg$fourier$n_starts)]
  }
  for(start in starts) encode_parameters(start$parameters,problem$map)
  starts
}
p5_check <- function(cfg,data) {
  p3_text(cfg$campaign_id,"campaign_id")
  if(!grepl("^[A-Za-z0-9_-]+$",cfg$campaign_id)) stop("Filename-safe campaign_id required")
  p3_text(cfg$analysis_type,"analysis_type")
  if(!cfg$analysis_type %in% c("descriptive","conditional_future","strict_future")) stop("Invalid analysis_type")
  if(!isTRUE(cfg$model$K==8L) || !isTRUE(cfg$model$lambda==1)) stop("Phase 5 SPL reference is fixed at K=8, lambda=1")
  validate_phase2_data(data)
  if(!identical(cfg$data$include_npi,FALSE) || !is.null(cfg$data$post_npi) || any(data$npieff!=1) ||
     length(data$cases)>data$definition$Nnonpi) stop("Phase 5 prohibits NPI effects and post-NPI observations")
  if(!is.null(cfg$expected_window)) {
    w <- cfg$expected_window
    if(!identical(data$definition$region,w$region) || length(data$cases)!=w$analysis_end ||
       data$processed_week_marker[1]!=w$first || tail(data$processed_week_marker,1)!=w$last)
      stop("Data do not match the prespecified region/window")
  }
  validate_starts_config(cfg$starts)
  if(length(cfg$starts$user_starts)) stop("Use fourier.user_starts for named Fourier overrides; SPL uses matched grid starts")
  integer_scalar(cfg$fourier$n_starts,"Fourier n_starts",8L)
  integer_scalar(cfg$fourier$seed,"Fourier seed",0L)
  scalar_number(cfg$fourier$coefficient_scale,"Fourier coefficient scale",0)
  if(!is.list(cfg$fourier$user_starts) || !identical(names(cfg$fourier$user_starts),c("A","AS"))) stop("fourier.user_starts must have A and AS entries")
  for(n in c("mesh","derivative_tolerance","root_tolerance")) scalar_number(cfg$peak_detection[[n]],paste("peak_detection",n),0,strict_lower=TRUE)
  if(cfg$peak_detection$mesh>1) stop("Beta maxima detection mesh must be <= 1 week")
  provenance <- p3_provenance(cfg,length(data$cases),FALSE)
  checks <- list()
  for(name in c("A","AS","SPL")) {
    problem <- prepare_estimation_problem(data,p5_model(cfg,name),cfg$estimation)
    starts <- p5_starts(cfg,problem,name)
    keys <- vapply(starts,function(x) paste(encode_parameters(x$parameters,problem$map),collapse="|"),character(1))
    checks[[name]] <- data.frame(model=name,dimension=problem$map$dimension,starts=length(starts),unique_starts=length(unique(keys)))
  }
  if(!is.logical(cfg$rolling$enabled) || length(cfg$rolling$enabled)!=1L || is.na(cfg$rolling$enabled)) stop("Set rolling.enabled")
  if(cfg$rolling$enabled) {
    if(!cfg$analysis_type %in% c("strict_future","conditional_future")) stop("Rolling requires strict_future or conditional_future")
    origins <- cfg$rolling$origins;finite_numeric(origins,"rolling origins")
    integer_scalar(cfg$rolling$horizon,"rolling horizon",1L)
    if(any(origins!=floor(origins) | origins<52) || any(diff(origins)<=0) || any(origins+cfg$rolling$horizon>length(data$cases)))
      stop("Rolling origins must increase and leave the complete declared horizon")
    for(origin in origins) {
      p3_provenance(cfg,origin,TRUE);p3_information(cfg,origin)
    }
  }
  list(type="phase5_check",status="CONFIG_CHECKED_NO_FIT",start_design=do.call(rbind,checks),
    provenance=provenance,data_audit=p3_data_audit(data,cfg),
    warning=p5_warning(cfg,provenance))
}
p5_warning <- function(cfg,provenance) {
  conditional <- any(provenance$category %in% c("full_data","user_conditional"))
  paste(c(if(conditional) "CONDITIONAL ON FULL-DATA/USER-CONDITIONAL NUISANCE INPUTS; not strict held-out evidence",
    cfg$limitations,"Approximate conditional LR; regularity and global optimum unverified; harmonic structure is not causal evidence"),collapse="; ")
}
p5_fit_models <- function(data,cfg,optimizer=stats::optim) {
  records <- list()
  for(name in c("A","AS","SPL")) {
    records[[name]] <- tryCatch({
      model <- p5_model(cfg,name);problem <- prepare_estimation_problem(data,model,cfg$estimation)
      starts <- p5_starts(cfg,problem,name);warm_status <- "NOT_APPLICABLE";warm_loglik_difference <- NA_real_
      if(name=="AS") {
        annual <- records$A$fit$representative
        if(is.null(annual)) warm_status <- "ANNUAL_UNAVAILABLE" else {
          warm <- annual$evaluation$parameters
          warm$free_coordinates <- c(unname(warm$free_coordinates),0,0)
          embedded <- evaluate_candidate(encode_parameters(warm,problem$map),problem)
          if(!isTRUE(embedded$valid)) stop("Nested annual warm endpoint is not valid in AS")
          warm_loglik_difference <- embedded$loglik-annual$evaluation$loglik
          if(abs(warm_loglik_difference)>1e-8*max(1,abs(annual$evaluation$loglik)))
            stop("Nested annual warm start does not reproduce annual likelihood")
          starts[[length(starts)+1L]] <- list(label="nested-annual-a2-b2-zero",parameters=warm)
          warm_status <- "ANNUAL_ENDPOINT_EMBEDDED_WITH_ZERO_SEMIANNUAL"
        }
      }
      fit <- fit_multistart(data,model,cfg$estimation,p5_start_config(cfg,name),
        cfg$initial_parameters,optimizer=optimizer,explicit_starts=starts)
      list(model=name,status=if(is.null(fit$representative)) fit$status else "AVAILABLE",fit=fit,
        warm_start=warm_status,warm_loglik_difference=warm_loglik_difference,
        reason=if(is.null(fit$representative)) "No converged valid endpoint" else "")
    },error=function(e) list(model=name,status="FAILED",fit=NULL,warm_start="UNAVAILABLE",reason=conditionMessage(e)))
  }
  p3_matched(records$A,records$AS);p3_matched(records$A,records$SPL)
  records
}
p5_comparison <- function(models) {
  do.call(rbind,lapply(names(models),function(name) {
    r <- models[[name]];fit <- r$fit$representative;available <- identical(r$status,"AVAILABLE")
    k <- if(available) fit$map$dimension else NA_integer_
    n <- if(available) sum(fit$comparison_context$mask) else NA_integer_
    loglik <- if(available) fit$evaluation$loglik else NA_real_
    data.frame(model=name,status=r$status,reason=r$reason,n=n,estimated_parameters=k,
      logLik=loglik,unpenalized_NLL=-loglik,
      AIC=if(name!="SPL" && available) 2*k-2*loglik else NA_real_,
      BIC=if(name!="SPL" && available) log(n)*k-2*loglik else NA_real_,
      AIC_BIC_definition=if(name=="SPL") "NOT_COMPARABLE: penalized spline effective degrees of freedom unavailable" else "Unpenalized Fourier maximum likelihood, conditional on fixed nuisance values",
      penalty=if(available) fit$evaluation$penalty else NA_real_,
      starts=length(r$fit$runs),converged=if(is.null(r$fit)) 0L else sum(vapply(r$fit$runs,function(x) identical(x$status,"CONVERGED_VALID"),logical(1))),
      boundary=if(available) paste(fit$boundary$coordinates,collapse=";") else "",
      convergence=if(available) fit$convergence else NA_integer_,
      nested_warm_status=r$warm_start,
      nested_warm_loglik_difference=if(is.null(r$warm_loglik_difference)) NA_real_ else r$warm_loglik_difference)
  }))
}
p5_lr <- function(table,tolerance=1e-8) {
  a <- table[table$model=="A",];as <- table[table$model=="AS",]
  available <- nrow(a)==1L && nrow(as)==1L && a$status=="AVAILABLE" && as$status=="AVAILABLE"
  comparable <- available && a$n==as$n && as$estimated_parameters-a$estimated_parameters==2L
  lr <- if(comparable) 2*(as$logLik-a$logLik) else NA_real_
  status <- if(!comparable) "UNAVAILABLE_OR_NOT_COMPARABLE" else if(lr< -tolerance) "AS_WORSE_OPTIMIZATION_UNRESOLVED" else
    if(nzchar(a$boundary) || nzchar(as$boundary)) "APPROXIMATE_WITH_BOUNDARY_WARNING" else "APPROXIMATE_REGULARITY_UNVERIFIED"
  data.frame(LR=lr,df=2L,p_value=if(comparable && lr>= -tolerance) stats::pchisq(max(0,lr),df=2,lower.tail=FALSE) else NA_real_,
    status=status,definition="Approximate conditional likelihood-ratio test given fixed nuisance values; not causal proof; no bootstrap calibration")
}
p5_fit <- function(cfg,data,optimizer=stats::optim) {
  models <- p5_fit_models(data,cfg,optimizer);table <- p5_comparison(models)
  list(type="phase5_fit",status=if(all(table$status=="AVAILABLE")) "ALL_MODELS_AVAILABLE" else "MODEL_FAILURES_RETAINED",
    models=models,comparison=table,lr=p5_lr(table),data_audit=p3_data_audit(data,cfg),
    markers=p3_markers(data),observed=data$cases,mask=data$likelihood_include,
    provenance=p3_provenance(cfg,length(data$cases),FALSE),warning=p5_warning(cfg,p3_provenance(cfg,length(data$cases),FALSE)))
}
p5_rolling <- function(cfg,data,optimizer=stats::optim) {
  if(!isTRUE(cfg$rolling$enabled)) stop("Rolling validation requires rolling.enabled=TRUE")
  folds <- lapply(cfg$rolling$origins,function(origin) {
    train <- p3_prefix(data,as.integer(origin));ids <- seq.int(origin+1L,origin+cfg$rolling$horizon)
    local <- cfg;local$future <- list(fold=list(values=rep(1,length(ids)),kind="known",citation="Original pre-NPI unit multiplier",known_at_index=0L,last_used_index=0L))
    policy <- p3_evaluation_policy(local,data,"fold",origin,ids)
    models <- p5_fit_models(train,cfg,optimizer)
    future <- p3_future_input(data,ids,local,"fold",origin)
    for(name in names(models)) models[[name]] <- tryCatch({
      r <- models[[name]]
      if(r$status=="AVAILABLE") r$forecast <- p3_score_forecast(p3_forecast(train,r$fit$representative,future),
        p3_outcomes(data,ids),list(incomplete_years="exclude",ties="first_chronological"))
      r
    },error=function(e) {r <- models[[name]];r$status <- "FORECAST_FAILED";r$reason <- conditionMessage(e);r})
    list(origin=origin,indices=ids,policy=policy,models=models)
  })
  list(type="phase5_rolling",status="ROLLING_FOLDS_STORED",folds=folds,
    warning="Overlapping test horizons are not independent; no strict held-out claim with full-data nuisance conditioning")
}
p5_identity <- function(root,cfg,config_path) {
  files <- unique(c(list.files(file.path(root,"R"),pattern="[.]R$",full.names=TRUE),
    file.path(root,"scripts/phase5.R"),file.path(root,"scripts/phase5.sh"),
    file.path(root,"config/phase2_defaults.R"),config_path,cfg$config_dependencies,p3_inputs(cfg)))
  if(any(!file.exists(files))) stop("Phase 5 provenance input missing: ",paste(files[!file.exists(files)],collapse=", "))
  tools::md5sum(normalizePath(files,mustWork=TRUE))
}
p5_inherit_config <- function(path) {
  env <- new.env(parent=parent.frame());sys.source(path,envir=env)
  cfg <- env$phase5_config()
  cfg$config_dependencies <- unique(c(cfg$config_dependencies,normalizePath(path,mustWork=TRUE)))
  cfg
}
