# Saved-only exports. No SIRS propagation, differentiation, optimizer or sampling.
p6_report <- function(result,out,diagnostics=NULL) {
  if(!result$type %in% c("phase6_fit","phase6_bayes")) stop("Phase 6 fit or Bayesian artifact required")
  if(!is.null(diagnostics) && (!identical(diagnostics$config,result$config) ||
    !identical(diagnostics$identity,result$identity) || !identical(diagnostics$source_fit_hash,result$artifact_hash)))
    stop("Diagnostics do not belong to this exact saved fit")
  csv <- function(x,n) utils::write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE)
  mat <- function(x,n) if(!is.null(x)) csv(data.frame(parameter=rownames(x),x,check.names=FALSE),n)
  cfg <- result$config
  fixed <- data.frame(parameter=c("pop","mu","gamma","npieff","period"),
    status=c("fixed",if(cfg$estimation$joint$mu) "free" else "fixed",if(cfg$estimation$joint$gamma) "free" else "fixed","fixed","fixed"),
    source=c(cfg$fixed_sources$pop,cfg$fixed_sources$mu,cfg$fixed_sources$gamma,cfg$fixed_sources$npieff,"Prespecified 52-week prepared seasonal calendar"))
  csv(fixed,"fixed_assumptions")
  if(result$type=="phase6_bayes") {
    if(!is.null(diagnostics) && diagnostics$type!="phase6_bayes_diagnostics") stop("Bayesian diagnostics required")
    if(!is.null(diagnostics)) {
      csv(diagnostics$summary,"posterior_summary");csv(diagnostics$chains,"sampler_diagnostics")
      csv(diagnostics$replay,"stan_R_replay")
      mat(diagnostics$correlations,"posterior_correlations")
      csv(diagnostics$predictive$table,"posterior_predictive")
    }
    # Algebraic summaries from saved joint draws only (no trajectories).
    m <- p6_draw_matrix(result$posterior$draws)
    problem <- p6_unbounded_problem(cfg,result$training,result$model)
    physical <- t(vapply(seq_len(nrow(m)),function(i) p6_physical(p6_draw_q(m[i,],problem),problem$map),
      numeric(length(p6_physical(p6_draw_q(m[1,],problem),problem$map)))))
    colnames(physical) <- names(p6_physical(p6_draw_q(m[1,],problem),problem$map))
    csv(data.frame(parameter=colnames(physical),mean=colMeans(physical),sd=apply(physical,2,stats::sd),
      q025=apply(physical,2,stats::quantile,.025),median=apply(physical,2,stats::median),
      q975=apply(physical,2,stats::quantile,.975),diagnostic_status=if(is.null(diagnostics)) "UNREVIEWED" else diagnostics$status),
      "posterior_parameter_uncertainty")
    annual <- sqrt(m[,"coefficients[1]"]^2+m[,"coefficients[2]"]^2)
    semi <- if(result$model=="AS") sqrt(m[,"coefficients[3]"]^2+m[,"coefficients[4]"]^2) else rep(0,nrow(m))
    csv(data.frame(model=result$model,draw=seq_len(nrow(m)),annual_log_amplitude=annual,semiannual_log_amplitude=semi),"posterior_harmonic_amplitudes")
    beta_names <- paste0("beta[",1:52,"]")
    if(all(beta_names %in% colnames(m))) {
      beta <- m[,beta_names,drop=FALSE]
      csv(data.frame(week=1:52,mean=colMeans(beta),q025=apply(beta,2,stats::quantile,.025),
        median=apply(beta,2,stats::median),q975=apply(beta,2,stats::quantile,.975)),"posterior_beta_curves")
    }
    # Continuous beta maxima are algebraic evaluations of saved joint coefficients.
    ids <- unique(as.integer(round(seq(1,nrow(m),length.out=min(cfg$bayes$predictive_draws,nrow(m))))))
    beta_peaks <- lapply(ids,function(i) {
      parameters <- decode_parameters(p6_draw_q(m[i,],problem),problem$map)
      fit <- list(evaluation=list(parameters=parameters),map=problem$map,basis=problem$basis)
      detection <- p5_beta_peaks(fit,cfg$peak_detection)
      data.frame(draw=i,status=detection$status,count=detection$count,
        week_locations=paste(detection$table$week_coordinate,collapse=";"))
    })
    csv(do.call(rbind,beta_peaks),"posterior_transmission_peaks")
    expected_names <- paste0("states[",seq_along(result$training$cases),",5]")
    if(all(expected_names %in% colnames(m))) {
      mean_expected <- colMeans(m[,expected_names,drop=FALSE]);ids <- which(result$training$likelihood_include)
      error <- mean_expected[ids]-result$training$cases[ids]
      csv(data.frame(index=seq_along(mean_expected),observed=result$training$cases,posterior_mean_expected=mean_expected,
        scored=result$training$likelihood_include),"posterior_cases_expected")
      csv(data.frame(n=length(ids),MAE=mean(abs(error)),RMSE=sqrt(mean(error^2)),bias=mean(error),
        definition="Descriptive fit of posterior mean expected counts; not likelihood at marginal posterior medians"),"case_fit_metrics")
    }
    writeLines(c("# Phase 6 Bayesian saved-result report","",
      paste("Model:",result$model,"Training end:",result$training_end),
      if(is.null(diagnostics)) "Diagnostics unavailable: run bayes-diagnose explicitly." else diagnostics$status,
      "Posterior amplitudes are not causal effects. No AIC/BIC or LR is computed from posterior summaries.",
      "Use predeclared blocked train-prefix predictive evaluation; no pointwise PSIS/LOO claim for dependent time points.",
      "All summaries are conditional on the proper design priors and fixed demographic/recovery assumptions."),file.path(out,"report.md"))
  } else {
    if(!is.null(diagnostics) && diagnostics$type!="phase6_diagnostics") stop("Joint MLE diagnostics required")
    csv(result$comparison,"model_comparison");csv(result$lr,"likelihood_ratio")
    csv(p5_coefficients(result$models),"fourier_coefficients");csv(p5_diagnostics(result$models),"all_start_diagnostics")
    csv(result$data_audit,"data_audit")
    parameters <- curves <- peaks <- amplitudes <- series <- metrics <- list()
    for(name in names(result$models)) {
      f <- result$models[[name]]$fit$representative;if(is.null(f)) next
      t <- f$parameter_record;t$model <- name;parameters[[name]] <- t
      curves[[name]] <- data.frame(model=name,week=1:52,beta=f$evaluation$beta_by_week)
      b <- p5_beta_peaks(f,cfg$peak_detection)
      t <- b$table;t$model <- rep(name,nrow(t));t$status <- rep(b$status,nrow(t));peaks[[name]] <- t
      z <- f$evaluation$parameters$free_coordinates
      if(name!="SPL") amplitudes[[name]] <- data.frame(model=name,annual_log_amplitude=sqrt(sum(z[1:2]^2)),
        semiannual_log_amplitude=if(name=="AS") sqrt(sum(z[3:4]^2)) else 0)
      series[[name]] <- data.frame(model=name,result$markers,observed=result$observed,scored=result$mask,
        expected=f$evaluation$states$expected_cases,incidence=f$evaluation$states$incidence,I=f$evaluation$states$I)
      ids <- which(result$mask)
      t <- p3_metrics(result$observed[ids],f$evaluation$states$expected_cases[ids],f$evaluation$parameters$phi)$summary
      t$model <- name;metrics[[name]] <- t
    }
    spl <- result$models$SPL$fit$representative
    if(!is.null(spl)) csv(data.frame(coefficient=seq_along(spl$evaluation$spline_coefficients),
      value=spl$evaluation$spline_coefficients,K=spl$basis$K,lambda=spl$lambda,
      constraint="c=c(z,-sum(z))"),"spline_coefficients")
    csv(p3_bind_rows(parameters),"parameter_estimates");csv(p3_bind_rows(curves),"beta_curves")
    csv(p3_bind_rows(peaks),"transmission_peaks");csv(p3_bind_rows(amplitudes),"harmonic_amplitudes")
    csv(p3_bind_rows(series),"cases_expected");csv(p3_bind_rows(metrics),"case_fit_metrics")
    if(!is.null(diagnostics)) for(name in names(diagnostics$records)) {
      d <- diagnostics$records[[name]]
      for(i in seq_along(d$hessians)) {
        h <- d$hessians[[i]];prefix <- paste(name,"hessian",i,sep="_")
        if(!is.null(h$component)) prefix <- paste(name,h$component,"hessian",i,sep="_")
        csv(data.frame(status=h$status,step=h$step,condition=if(is.null(h$condition)) NA_real_ else h$condition,
          reason=if(is.null(h$reason)) "" else h$reason),paste0(prefix,"_status"))
        csv(data.frame(eigenvalue=h$eigenvalues),paste0(prefix,"_eigenvalues"))
        mat(h$hessian,paste0(prefix,"_matrix"));mat(h$correlation,paste0(prefix,"_correlation"))
        csv(h$uncertainty,paste0(prefix,"_local_uncertainty"))
        mat(h$geometry_covariance,paste0(prefix,"_inverse_curvature_geometry_NOT_frequentist_covariance"))
        if(!is.null(h$component)) csv(data.frame(component=h$component,
          interpretation="Curvature at saved penalized endpoint; NOT ordinary frequentist covariance"),paste0(prefix,"_definition"))
      }
      for(n in names(d$profiles)) csv(d$profiles[[n]]$table,paste(name,n,"profile",sep="_"))
      mat(d$endpoint_correlation,paste0(name,"_near_optimal_start_correlation"))
    }
    writeLines(c("# Phase 6 independent joint-MLE saved-result report","",
      paste("Campaign:",cfg$campaign_id,"| scored observations:",sum(result$mask)),
      "A: annual log-Fourier; AS: annual + semiannual log-Fourier. Same SIRS/NB2/data/mask.",
      if("SPL" %in% names(result$models)) "SPL: independent joint penalized cyclic spline. AIC/BIC are NA; no AS-vs-SPL LR. Unpenalized fit is descriptive, not complexity-adjusted superiority.",
      "All estimated nuisance coordinates are counted in Fourier AIC/BIC. These describe selected local endpoints.",
      "LR nominal_p_value is an exploratory chi-square(2) reference only; calibrated p_value is NA.",
      "Global optimality, regularity and identifiability remain unresolved; numerical convergence is insufficient.",
      if(is.null(diagnostics)) "Uncertainty/identifiability diagnostics unavailable: run diagnose explicitly." else
        "Profiles reoptimize nuisance coordinates; inspect better-than-saved MLE flags, failures, bounds and Hessian step stability. Local SEs are approximate.",
      "Amplitudes are log-transmission harmonics, not causal effects or automatic evidence of two case peaks.",
      cfg$limitations),file.path(out,"report.md"))
  }
  list(type="phase6_report",status="SAVED_ONLY_REPORT_WRITTEN",files=list.files(out))
}
