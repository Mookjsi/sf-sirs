# Proper independent coordinate priors; optional cmdstanr backend is explicit only.
p6_prior_names <- function() c("log_beta0","log_phi","log_S_over_R","log_I_over_R","logit_rho","log_omega","log_delta","log_mu","log_gamma","coefficient")
p6_validate_priors <- function(cfg) {
  if(!identical(sort(names(cfg$priors)),sort(p6_prior_names()))) stop("Explicit proper priors required for every coordinate family")
  for(n in p6_prior_names()) {
    v <- cfg$priors[[n]];finite_numeric(v,paste("prior",n),2L)
    if(v[2]<=0) stop("Proper normal prior scale must be positive")
  }
  invisible(NULL)
}
p6_unbounded_problem <- function(cfg,data,name) {
  if(!name %in% c("A","AS")) stop("Bayesian Phase 6 interfaces support Fourier A/AS only; SPL Bayesian inference is not implemented")
  p <- p6_problem(cfg,data,name)
  p$map$lower[] <- -Inf;p$map$upper[] <- Inf
  p
}
p6_prior_predictive <- function(cfg,data,name) {
  p6_validate_priors(cfg);integer_scalar(cfg$prior_predictive$n,"prior draws",1L)
  problem <- p6_unbounded_problem(cfg,data,name)
  records <- with_phase2_seed(cfg$prior_predictive$seed,function() lapply(seq_len(cfg$prior_predictive$n),function(i) {
    q <- setNames(vapply(problem$map$order,function(n) {
      v <- cfg$priors[[if(n %in% problem$map$shape_order) "coefficient" else n]]
      stats::rnorm(1,v[1],v[2])
    },numeric(1)),problem$map$order)
    e <- tryCatch(p6_propagate(q,problem),error=function(e) e)
    if(inherits(e,"error")) return(list(draw=i,status="INVALID_PRIOR_DRAW",coordinates=q,reason=conditionMessage(e)))
    y <- stats::rnbinom(length(data$cases),size=e$parameters$phi,mu=e$states$expected_cases)
    list(draw=i,status=if(all(is.finite(y))) "GENERATED" else "NB2_RNG_FAILED",coordinates=q,
      parameters=e$parameters,expected=e$states$expected_cases,replicated=y)
  }))
  list(type="phase6_prior_predictive",status="PRIOR_DRAWS_STORED_WITH_FAILURES",model=name,records=records,
    warning="No rejection/resampling to hide implausible prior predictions. Review before sampling; these priors are design assumptions.")
}
p6_stan_data <- function(cfg,data,name) {
  p6_validate_priors(cfg)
  if(any(data$cases>.Machine$integer.max)) stop("Stan count exceeds signed integer support")
  list(N=length(data$cases),H=if(name=="A") 1L else 2L,y=as.integer(data$cases),week=as.integer(data$seasonal_week),
    pop=data$fixed_parameters$pop,mu_fixed=cfg$fixed$mu,gamma_fixed=cfg$fixed$gamma,
    free_mu=as.integer(cfg$estimation$joint$mu),free_gamma=as.integer(cfg$estimation$joint$gamma),
    prior_mean=vapply(cfg$priors[p6_prior_names()],`[`,numeric(1),1L),
    prior_sd=vapply(cfg$priors[p6_prior_names()],`[`,numeric(1),2L))
}
p6_cmdstan_backend <- function(stan_data,cfg,source,out) {
  if(!requireNamespace("cmdstanr",quietly=TRUE) || !requireNamespace("posterior",quietly=TRUE))
    stop("Bayesian sampling requires user-installed cmdstanr, CmdStan and posterior; no installation is automatic")
  target <- file.path(out,"phase6_fourier_sirs.stan")
  if(!file.copy(source,target)) stop("Could not preserve Stan source in new output")
  model <- cmdstanr::cmdstan_model(target)
  b <- cfg$bayes
  fit <- model$sample(data=stan_data,seed=b$seed,init=p6_stan_init(cfg,stan_data$H),chains=b$chains,parallel_chains=b$parallel_chains,
    iter_warmup=b$iter_warmup,iter_sampling=b$iter_sampling,adapt_delta=b$adapt_delta,
    max_treedepth=b$max_treedepth,output_dir=out,refresh=100,sig_figs=17)
  # Retain joint draws and chain ordering; no CmdStan external pointer required later.
  list(draws=as.array(fit$draws(format="draws_array")),sampler=as.array(fit$sampler_diagnostics(format="draws_array")),
    versions=list(cmdstan=as.character(cmdstanr::cmdstan_version()),cmdstanr=as.character(utils::packageVersion("cmdstanr"))),
    stan_source_hash=tools::md5sum(target),csv_files=fit$output_files())
}
p6_stan_init <- function(cfg,H) {
  # Independent prior-center jitter, not an MLE or original posterior warm start.
  with_phase2_seed(cfg$bayes$seed,function() lapply(seq_len(cfg$bayes$chains),function(i) {
    draw <- function(n) stats::rnorm(1,cfg$priors[[n]][1],.1)
    list(log_beta0=draw("log_beta0"),log_phi=draw("log_phi"),
      initial_logratios=c(draw("log_S_over_R"),draw("log_I_over_R")),
      logit_rho=draw("logit_rho"),log_omega=draw("log_omega"),log_delta=draw("log_delta"),
      coefficients=stats::rnorm(2L*H,cfg$priors$coefficient[1],.1),
      log_mu=if(cfg$estimation$joint$mu) draw("log_mu") else numeric(),
      log_gamma=if(cfg$estimation$joint$gamma) draw("log_gamma") else numeric())
  }))
}
p6_bayes <- function(cfg,data,name,source,out,train_end=NULL,backend=p6_cmdstan_backend) {
  if(!name %in% c("A","AS")) stop("Select --model A or AS")
  if(!is.null(train_end)) {
    integer_scalar(train_end,"training end",2L)
    if(train_end>length(data$cases)) stop("Training end beyond data")
    data <- p3_prefix(data,train_end)
  }
  s <- backend(p6_stan_data(cfg,data,name),cfg,source,out)
  if(length(dim(s$draws))!=3L || length(dim(s$sampler))!=3L) stop("Backend must return iteration x chain x variable arrays")
  list(type="phase6_bayes",status="DRAWS_STORED_NOT_DIAGNOSTICALLY_CERTIFIED",model=name,
    training_end=length(data$cases),training=data,posterior=s,
    config=cfg,prior_definition="Proper normals on unconstrained coordinates; no MLE guard truncation or Jacobian for these coordinate priors")
}
p6_draw_matrix <- function(x) {
  d <- dim(x);matrix(x,nrow=d[1]*d[2],ncol=d[3],dimnames=list(NULL,dimnames(x)[[3]]))
}
p6_draw_q <- function(row,problem) {
  get <- function(n) {if(!n %in% names(row)) stop("Missing posterior coordinate ",n);unname(row[n])}
  q <- setNames(numeric(problem$map$dimension),problem$map$order)
  for(n in names(q)) q[n] <- if(n %in% problem$map$shape_order) get(paste0("coefficients[",match(n,problem$map$shape_order),"]")) else
    if(n=="log_S_over_R") get("initial_logratios[1]") else if(n=="log_I_over_R") get("initial_logratios[2]") else
    if(n %in% c("log_mu","log_gamma")) get(paste0(n,"[1]")) else get(n)
  q
}
p6_posterior_predictive <- function(result,data=result$training,indices=which(data$likelihood_include)) {
  cfg <- result$config;problem <- p6_unbounded_problem(cfg,data,result$model)
  m <- p6_draw_matrix(result$posterior$draws)
  n <- min(nrow(m),cfg$bayes$predictive_draws)
  integer_scalar(n,"predictive draws",1L)
  selected <- with_phase2_seed(cfg$bayes$seed+1L,function() sample.int(nrow(m),n))
  records <- with_phase2_seed(cfg$bayes$seed+2L,function() lapply(selected,function(i) tryCatch({
    q <- p6_draw_q(m[i,],problem);e <- p6_propagate(q,problem)
    y <- stats::rnbinom(length(indices),size=e$parameters$phi,mu=e$states$expected_cases[indices])
    if(any(!is.finite(y))) stop("NB2 RNG failed")
    list(draw=i,status="GENERATED",expected=e$states$expected_cases[indices],replicated=y,
      log_lik=stats::dnbinom(data$cases[indices],mu=e$states$expected_cases[indices],size=e$parameters$phi,log=TRUE))
  },error=function(e) list(draw=i,status="FAILED",reason=conditionMessage(e)))))
  ok <- vapply(records,function(r) r$status=="GENERATED",logical(1))
  table <- scores <- NULL
  # Any failed selected draw makes aggregated predictive summaries unavailable.
  if(all(ok)) {
    y <- do.call(rbind,lapply(records,`[[`,"replicated"));expected <- do.call(rbind,lapply(records,`[[`,"expected"))
    table <- data.frame(index=indices,observed=data$cases[indices],mean=colMeans(expected),
      q025=apply(y,2,stats::quantile,.025),median=apply(y,2,stats::median),q975=apply(y,2,stats::quantile,.975))
    L <- do.call(rbind,lapply(records,`[[`,"log_lik"))
    lme <- function(v) {a <- max(v);if(a== -Inf) -Inf else a+log(mean(exp(v-a)))}
    scores <- data.frame(index=indices,integrated_log_predictive=apply(L,2,lme))
  }
  list(status=if(all(ok)) "JOINT_DRAW_PREDICTIVE_STORED" else "FAILED_DRAWS_NO_AGGREGATION",records=records,table=table,scores=scores,
    warning="Posterior predictive summaries are joint-draw simulations, not simulation at marginal medians; interval coverage is descriptive")
}
p6_sampler_diagnostics <- function(result) {
  s <- result$posterior$sampler;names <- dimnames(s)[[3]]
  required <- c("divergent__","treedepth__","energy__")
  if(!all(required %in% names)) stop("Sampler diagnostics incomplete")
  chains <- lapply(seq_len(dim(s)[2]),function(j) {
    energy <- s[,j,"energy__"];variance <- stats::var(energy)
    data.frame(chain=j,divergences=sum(s[,j,"divergent__"]),
      max_depth_hits=sum(s[,j,"treedepth__"]>=result$config$bayes$max_treedepth),
      BFMI=if(is.finite(variance) && variance>0) mean(diff(energy)^2)/variance else NA_real_)
  })
  do.call(rbind,chains)
}
p6_bayes_diagnose <- function(result) {
  if(!requireNamespace("posterior",quietly=TRUE)) stop("Install posterior explicitly to compute rank R-hat and ESS")
  a <- posterior::as_draws_array(result$posterior$draws)
  vars <- dimnames(result$posterior$draws)[[3]]
  # All free coordinates and derived physical parameters; omit deterministic trajectories.
  keep <- !grepl("^(states|beta\\[|g\\[|log_lik\\[)",vars)
  summary <- as.data.frame(posterior::summarise_draws(posterior::subset_draws(a,variable=vars[keep])))
  chains <- p6_sampler_diagnostics(result)
  p <- p6_unbounded_problem(result$config,result$training,result$model)
  m <- p6_draw_matrix(result$posterior$draws)
  physical <- t(vapply(seq_len(nrow(m)),function(i) p6_physical(p6_draw_q(m[i,],p),p$map),numeric(length(p6_physical(p6_draw_q(m[1,],p),p$map)))))
  colnames(physical) <- names(p6_physical(p6_draw_q(m[1,],p),p$map))
  variable <- apply(physical,2,stats::sd)>0
  correlations <- if(any(variable)) stats::cor(physical[,variable,drop=FALSE]) else NULL
  # A stuck free coordinate must FAIL, even when all chains give the same constant.
  required <- c("log_beta0","log_phi","initial_logratios[1]","initial_logratios[2]","logit_rho","log_omega","log_delta",
    paste0("coefficients[",seq_len(p$basis$K-1L),"]"),
    if(result$config$estimation$joint$mu) "log_mu[1]",if(result$config$estimation$joint$gamma) "log_gamma[1]")
  monitor <- summary$variable %in% required
  replay <- p6_stan_replay(result)
  pass <- dim(a)[2]>=4L && all(required %in% summary$variable) && isTRUE(all(is.finite(summary$rhat[monitor]) & summary$rhat[monitor]<=1.01 &
    summary$ess_bulk[monitor]>=400 & summary$ess_tail[monitor]>=400) &&
    all(chains$divergences==0 & chains$max_depth_hits==0 & is.finite(chains$BFMI) & chains$BFMI>=.3)) && all(replay$status=="SAVED_STAN_R_REPLAY_MATCH")
  list(type="phase6_bayes_diagnostics",status=if(pass) "NUMERICAL_THRESHOLDS_PASSED_NOT_MODEL_VALIDATION" else "REVIEW_REQUIRED",
    summary=summary,chains=chains,correlations=correlations,replay=replay,predictive=p6_posterior_predictive(result),
    config=result$config,identity=result$identity,source_fit_hash=result$artifact_hash)
}
p6_stan_replay <- function(result) {
  problem <- p6_unbounded_problem(result$config,result$training,result$model)
  m <- p6_draw_matrix(result$posterior$draws);N <- length(result$training$cases)
  variables <- as.vector(outer(seq_len(N),1:5,function(i,j) paste0("states[",i,",",j,"]")))
  ids <- unique(as.integer(round(seq(1,nrow(m),length.out=min(20L,nrow(m))))))
  do.call(rbind,lapply(ids,function(i) tryCatch({
    if(!all(variables %in% colnames(m))) stop("Saved Stan states missing; numerical port verification unavailable")
    e <- p6_propagate(p6_draw_q(m[i,],problem),problem)
    R <- as.numeric(as.matrix(e$states[,c("S","I","R","incidence","expected_cases")]))
    stan <- unname(m[i,variables]);difference <- max(abs(R-stan))
    pass <- all(is.finite(stan)) && all(abs(R-stan)<=1e-5+1e-8*pmax(abs(R),abs(stan)))
    data.frame(draw=i,status=if(pass) "SAVED_STAN_R_REPLAY_MATCH" else "STAN_R_MISMATCH",max_absolute_difference=difference,reason="")
  },error=function(e) data.frame(draw=i,status="REPLAY_UNAVAILABLE",max_absolute_difference=NA_real_,reason=conditionMessage(e)))))
}
p6_bayes_predict <- function(result,data) {
  cfg <- result$config;p <- cfg$prediction
  integer_scalar(p$origin,"prediction origin",2L);integer_scalar(p$horizon,"prediction horizon",1L)
  p3_text(p$protocol_id,"predeclared predictive protocol ID")
  if(result$training_end!=p$origin || p$origin+p$horizon>length(data$cases) ||
     !identical(result$training$cases,data$cases[seq_len(p$origin)]) ||
     !identical(result$training$seasonal_week,data$seasonal_week[seq_len(p$origin)]) ||
     !identical(result$training$fixed_parameters,data$fixed_parameters)) stop("Training-prefix/origin identity mismatch")
  ids <- seq.int(p$origin+1L,p$origin+p$horizon)
  prediction <- p6_posterior_predictive(result,p3_prefix(data,max(ids)),ids)
  list(type="phase6_bayes_predict",status=prediction$status,model=result$model,origin=p$origin,
    protocol_id=p$protocol_id,indices=ids,prediction=prediction,
    warning="Blocked posterior predictive continuation without state updates; compare identical declared folds/priors only after sampler review; overlapping folds dependent")
}
