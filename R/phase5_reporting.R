# Unpenalized comparison and descriptive beta diagnostics from saved fit endpoints.
p5_beta_evaluator <- function(fit) {
  p <- fit$evaluation$parameters
  if(isTRUE(fit$map$fourier)) {
    function(x,derivative=0L) fourier_profile(p$free_coordinates,p$beta0,fit$basis,x,derivative)$beta_by_week
  } else {
    coefficients <- constrained_coefficients(p$free_coordinates,fit$basis)
    function(x,derivative=0L) evaluate_spline_beta(x,coefficients,p$beta0,fit$basis,derivative)
  }
}
p5_beta_peaks <- function(fit,settings) {
  evaluate <- p5_beta_evaluator(fit)
  mesh <- sort(unique(c(seq(0,52,by=settings$mesh),52)))
  derivative <- evaluate(mesh,1L);tol <- settings$derivative_tolerance
  if(max(abs(derivative))<=tol) return(list(status="CONSTANT_OR_NUMERICALLY_FLAT",count=0L,
    table=data.frame(x=numeric(),week_coordinate=numeric(),beta=numeric())))
  roots <- numeric()
  for(i in seq_len(length(mesh)-1L)) {
    if(abs(derivative[i])<=tol) roots <- c(roots,mesh[i])
    if(derivative[i]*derivative[i+1L]<0) roots <- c(roots,stats::uniroot(function(x) evaluate(x,1L),
      mesh[i:(i+1L)],tol=settings$root_tolerance)$root)
  }
  roots <- sort(roots %% 52)
  if(length(roots)) roots <- roots[c(TRUE,diff(roots)>10*settings$root_tolerance)]
  if(length(roots)>1 && 52-tail(roots,1)+roots[1]<=10*settings$root_tolerance) roots <- head(roots,-1L)
  probe <- settings$mesh/10
  keep <- vapply(roots,function(x) evaluate(x-probe,1L)>tol && evaluate(x+probe,1L)< -tol,logical(1))
  roots <- roots[keep]
  list(status="NUMERICAL_LOCAL_MAXIMA",count=length(roots),
    table=data.frame(x=roots,week_coordinate=roots+1,beta=if(length(roots)) evaluate(roots) else numeric()),
    limitation="Finite mesh/sign detection with analytic derivatives; flat/tangential stationary regions may be unresolved")
}
p5_coefficients <- function(models) {
  rows <- list()
  for(name in intersect(c("A","AS"),names(models))) {
    fit <- models[[name]]$fit$representative
    if(is.null(fit)) next
    e <- fit$evaluation;coeff <- setNames(e$parameters$free_coordinates,fit$basis$coefficient_names)
    rows[[length(rows)+1L]] <- data.frame(model=name,parameter=c("beta0","a0",names(coeff)),
      value=c(e$parameters$beta0,e$a0,coeff),status=c("estimated","derived_centered_intercept",rep("estimated",length(coeff))),
      definition="x=w-1; log beta=log(beta0)+g(x)-log(mean(exp(g(0:51))))")
  }
  p3_bind_rows(rows)
}
p5_diagnostics <- function(models) {
  rows <- list()
  scalar <- function(x) if(is.null(x) || !length(x)) NA_real_ else x
  for(name in names(models)) {
    record <- models[[name]];runs <- record$fit$runs
    if(!length(runs)) runs <- list(list(status=record$status,message=record$reason))
    for(i in seq_along(runs)) {
      r <- runs[[i]]
      rows[[length(rows)+1L]] <- data.frame(model=name,start=i,start_id=if(is.null(r$start_label)) "" else r$start_label,
        status=r$status,convergence=scalar(r$convergence),objective=scalar(r$evaluation$objective),
        unpenalized_NLL= -scalar(r$evaluation$loglik),penalty=scalar(r$evaluation$penalty),
        message=paste(c(r$message,record$reason),collapse="; "),warnings=paste(r$warnings,collapse="; "),
        boundary=paste(r$boundary$coordinates,collapse="; "),evaluations=scalar(r$failure_diagnostics$evaluations),
        invalid_evaluations=scalar(r$failure_diagnostics$invalid),
        initial_coordinates=paste(r$initial_coordinates,collapse=";"),final_coordinates=paste(r$coordinates,collapse=";"),
        lower=paste(r$map$lower,collapse=";"),upper=paste(r$map$upper,collapse=";"))
    }
  }
  p3_bind_rows(rows)
}
p5_report <- function(result,out,rolling=NULL) {
  if(!identical(result$type,"phase5_fit")) stop("Report requires a phase5_fit artifact")
  if(!is.null(rolling) && (!identical(rolling$type,"phase5_rolling") || !identical(rolling$config,result$config) ||
                         !identical(rolling$identity,result$identity))) stop("Rolling report identity/config mismatch")
  cfg <- result$config;table <- result$comparison;models <- result$models
  p3_matched(models$A,models$AS);p3_matched(models$A,models$SPL)
  curves <- list();peaks <- list();amplitudes <- list();parameters <- list();case_series <- list()
  case_metrics <- list();case_peaks <- list()
  for(i in seq_len(nrow(table))) {
    name <- table$model[i];fit <- models[[name]]$fit$representative
    table$beta_peak_count[i] <- NA_integer_;table$beta_peak_locations[i] <- ""
    table$beta_peak_status[i] <- "UNAVAILABLE"
    if(is.null(fit)) next
    detection <- p5_beta_peaks(fit,cfg$peak_detection)
    table$beta_peak_count[i] <- detection$count
    table$beta_peak_locations[i] <- paste(detection$table$week_coordinate,collapse=";")
    table$beta_peak_status[i] <- detection$status
    curves[[name]] <- data.frame(model=name,week=1:52,x=0:51,beta=fit$evaluation$beta_by_week,
      shape=fit$evaluation$beta_by_week/fit$evaluation$parameters$beta0)
    if(nrow(detection$table)) {t <- detection$table;t$model <- name;peaks[[name]] <- t}
    t <- fit$parameter_record;t$model <- name;parameters[[name]] <- t
    case_series[[name]] <- data.frame(model=name,result$markers,observed=result$observed,
      expected=fit$evaluation$states$expected_cases,incidence=fit$evaluation$states$incidence,
      prevalence_I=fit$evaluation$states$I,scored=result$mask,check.names=FALSE)
    ids <- which(result$mask)
    t <- p3_metrics(result$observed[ids],fit$evaluation$states$expected_cases[ids],fit$evaluation$parameters$phi)$summary
    t$model <- name;case_metrics[[name]] <- t
    t <- p3_peaks(result$markers[ids,,drop=FALSE],result$observed[ids],fit$evaluation$states$expected_cases[ids])
    t$model <- name;case_peaks[[name]] <- t
    if(name!="SPL") {
      coefficients <- fit$evaluation$parameters$free_coordinates
      amplitudes[[name]] <- data.frame(model=name,annual_log_amplitude=sqrt(coefficients[1]^2+coefficients[2]^2),
        semiannual_log_amplitude=if(name=="AS") sqrt(coefficients[3]^2+coefficients[4]^2) else 0,
        definition="Amplitudes of log-transmission harmonics; not case amplitudes or a causal effect")
    }
  }
  csv <- function(x,file) utils::write.csv(x,file.path(out,file),row.names=FALSE)
  csv(table,"model_comparison.csv");csv(result$lr,"likelihood_ratio.csv")
  csv(p5_coefficients(models),"fourier_coefficients.csv")
  beta <- p3_bind_rows(curves);csv(beta,"beta_curves.csv")
  csv(p3_bind_rows(peaks),"beta_local_maxima.csv");csv(p3_bind_rows(amplitudes),"harmonic_amplitudes.csv")
  csv(p5_diagnostics(models),"fit_diagnostics.csv");csv(p3_bind_rows(parameters),"parameter_estimates.csv")
  csv(p3_bind_rows(case_series),"cases_expected.csv");csv(result$data_audit,"data_audit.csv")
  csv(p3_bind_rows(case_metrics),"case_fit_metrics.csv");csv(p3_bind_rows(case_peaks),"yearly_case_peaks.csv")
  csv(result$provenance,"fixed_provenance.csv")
  grDevices::png(file.path(out,"beta_curves.png"),width=1200,height=750,res=120)
  tryCatch({
    if(!"beta" %in% names(beta)) {graphics::plot.new();graphics::title("All fitted beta profiles unavailable")}
    else {
      colors <- c(A="blue",AS="red",SPL="darkgreen")
      graphics::plot(range(beta$week),range(beta$beta),type="n",xlab="Prepared seasonal week w (x=w-1)",
        ylab="Transmission rate (week^-1)",main=paste(cfg$data$region,cfg$campaign_id),
        sub="Conditional descriptive profiles; harmonic structure is not a transmission mechanism")
      for(name in names(curves)) graphics::lines(curves[[name]]$week,curves[[name]]$beta,col=colors[[name]],lwd=2)
      graphics::legend("topright",names(curves),col=colors[names(curves)],lty=1,lwd=2,bty="n")
    }
  },finally=grDevices::dev.off())
  rolling_rows <- list()
  if(!is.null(rolling)) for(fold in rolling$folds) for(name in names(fold$models)) {
    r <- fold$models[[name]]
    t <- if(identical(r$status,"AVAILABLE") && !is.null(r$forecast)) r$forecast$scores$summary else data.frame(reason=r$reason)
    t$model <- name;t$origin <- fold$origin;t$first_test <- min(fold$indices);t$last_test <- max(fold$indices)
    t$status <- r$status;t$evaluation_class <- if(r$status=="AVAILABLE") fold$policy$class else "UNAVAILABLE"
    rolling_rows[[length(rolling_rows)+1L]] <- t
  }
  csv(p3_bind_rows(rolling_rows),"rolling_metrics.csv")
  lr <- result$lr
  lines <- c(paste0("# Phase 5: ",cfg$campaign_id),"",
    paste("Region:",cfg$data$region,"| prepared coverage:",result$data_audit$first_marker,"through",result$data_audit$last_marker),
    paste("Observations scored:",sum(result$mask),"(initial row excluded)"),"",
    "Models: A annual log-Fourier; AS annual + semiannual log-Fourier; SPL K=8, lambda=1 cyclic spline.",
    "Normalization: beta0 is the arithmetic mean over 52 prepared weeks; x=w-1; a0 is derived.","",
    paste("LR:",lr$LR,"| df:",lr$df,"| approximate conditional p-value:",lr$p_value,"|",lr$status),"",
    result$warning,"",
    "## Methodological limitations","",
    "- All fits share cases, masks, dynamics, fixed nuisance policy, NB2 and optimizer controls. AS adds an annual nested warm start.",
    "- Fourier models have no roughness penalty; AIC/BIC include beta0, harmonic coefficients, phi and any free nuisance coordinates.",
    "- SPL unpenalized NLL is descriptive; naive spline AIC/BIC are NA because effective degrees of freedom are unavailable.",
    "- Global optimization and chi-squared regularity are unverified. Boundary hits and singular/weak identification can invalidate nominal LR calibration.",
    "- A significantly nonzero second harmonic in log(beta) does not prove two beta peaks, case peaks, a causal mechanism or predictive improvement.",
    "- Even annual-only log forcing creates higher spectral components in beta through exponentiation, and nonlinear SIRS can create complex case dynamics.",
    "- Local beta maxima use an analytic-derivative sign mesh with seam deduplication; flat or tangential regions can remain unresolved.",
    "- Full-data nuisance medians condition the comparison on observations beyond the fit/forecast origin. Source BFMI warning remains unresolved.",
    "- Rolling folds use training-only likelihood/parameter fitting and no teacher forcing. Full-data nuisance conditioning prevents strict held-out interpretation; overlapping horizons are dependent.",
    "- Existing Phase 3 artifacts are not modified. Report reads saved endpoints; no fit or SIRS propagation occurs here.","",
    "See model_comparison.csv, likelihood_ratio.csv, harmonic_amplitudes.csv, beta_local_maxima.csv and fit_diagnostics.csv.")
  writeLines(lines,file.path(out,"report.md"))
  list(type="phase5_report",status="REPORT_WRITTEN_NO_FIT",comparison=table,lr=lr,files=list.files(out))
}

# Optional periodic prominence layer. Existing derivative-root finder is unchanged.
p5_periodic_prominence <- function(peaks,evaluate,mesh=.1,period=52,height_tolerance=1e-10) {
  scalar_number(mesh,"prominence mesh",0,strict_lower=TRUE)
  scalar_number(period,"period",0,strict_lower=TRUE)
  scalar_number(height_tolerance,"height tolerance",0)
  peaks$prominence <- peaks$relative_prominence <- numeric(nrow(peaks))
  for(i in seq_len(nrow(peaks))) {
    x <- peaks$x[i];height <- peaks$beta[i]
    higher <- which(peaks$beta>height+height_tolerance)
    # Equal highest peaks use a full turn in both directions; no seam endpoint.
    left <- if(length(higher)) min((x-peaks$x[higher])%%period) else period
    right <- if(length(higher)) min((peaks$x[higher]-x)%%period) else period
    floor <- function(distance,direction) {
      steps <- sort(unique(c(seq(0,distance,by=mesh),distance)))
      min(evaluate((x+direction*steps)%%period))
    }
    base <- max(floor(left,-1),floor(right,1))
    peaks$prominence[i] <- max(0,height-base)
    peaks$relative_prominence[i] <- peaks$prominence[i]/height
  }
  peaks
}
p5_beta_periodic_peaks <- function(fit,settings,min_prominence=0,min_relative_prominence=0) {
  scalar_number(min_prominence,"minimum prominence",0)
  scalar_number(min_relative_prominence,"minimum relative prominence",0,1)
  result <- p5_beta_peaks(fit,settings)
  result$table <- p5_periodic_prominence(result$table,p5_beta_evaluator(fit),settings$mesh)
  result$filtered <- result$table[result$table$prominence>=min_prominence &
    result$table$relative_prominence>=min_relative_prominence,,drop=FALSE]
  result$prominence_definition <- "Sampled periodic topographic prominence; traverse to nearest strictly higher detected maximum in each direction, or a full turn for global/equal-highest maxima"
  result
}
p5_circular_peak_match <- function(reference,target,max_shift=4,period=52,edge_allowed=NULL,edge_cost=NULL) {
  scalar_number(period,"peak period",0,strict_lower=TRUE)
  scalar_number(max_shift,"peak match shift",0,period/2)
  if(length(reference)) finite_numeric(reference,"reference peaks")
  if(length(target)) finite_numeric(target,"target peaks")
  a <- order(reference%%period);b <- order(target%%period)
  reference <- (reference%%period)[a];target <- (target%%period)[b]
  if(length(target)>16L) stop("Exact peak matcher supports at most 16 target peaks")
  distance <- abs(outer(reference,target,`-`));distance <- pmin(distance,period-distance)
  allowed <- matrix(TRUE,length(reference),length(target));cost <- distance
  if(!is.null(edge_allowed)) {
    if(!is.logical(edge_allowed) || anyNA(edge_allowed) || !identical(dim(edge_allowed),dim(distance))) stop("Invalid peak edge mask")
    allowed <- edge_allowed[a,b,drop=FALSE]
  }
  if(!is.null(edge_cost)) {
    if(!is.numeric(edge_cost) || any(!is.finite(edge_cost)) || any(edge_cost<0) || !identical(dim(edge_cost),dim(distance))) stop("Invalid peak edge costs")
    cost <- edge_cost[a,b,drop=FALSE]
  }
  memo <- new.env(parent=emptyenv())
  solve <- function(i,mask) {
    if(i>length(reference)) return(list(n=0L,cost=0,pairs=matrix(integer(),0,2)))
    key <- paste(i,mask,sep=":");if(exists(key,memo,inherits=FALSE)) return(get(key,memo))
    best <- solve(i+1L,mask)
    for(j in seq_along(target)) {
      bit <- bitwShiftL(1L,j-1L)
      if(bitwAnd(mask,bit)!=0L || distance[i,j]>max_shift || !allowed[i,j]) next
      next_match <- solve(i+1L,bitwOr(mask,bit))
      candidate <- list(n=next_match$n+1L,cost=next_match$cost+cost[i,j],pairs=rbind(c(i,j),next_match$pairs))
      if(candidate$n>best$n || (candidate$n==best$n && candidate$cost<best$cost-1e-12)) best <- candidate
    }
    assign(key,best,memo);best
  }
  pairs <- solve(1L,0L)$pairs
  rows <- list()
  if(nrow(pairs)) for(i in seq_len(nrow(pairs))) {
    a <- reference[pairs[i,1]];b <- target[pairs[i,2]]
    rows[[length(rows)+1L]] <- data.frame(reference_x=a,target_x=b,circular_distance=min(abs(a-b),period-abs(a-b)),
      signed_shift=(b-a+period/2)%%period-period/2,status="MATCHED")
  }
  for(i in setdiff(seq_along(reference),pairs[,1])) rows[[length(rows)+1L]] <- data.frame(reference_x=reference[i],target_x=NA_real_,circular_distance=NA_real_,signed_shift=NA_real_,status="UNMATCHED_REFERENCE")
  for(j in setdiff(seq_along(target),pairs[,2])) rows[[length(rows)+1L]] <- data.frame(reference_x=NA_real_,target_x=target[j],circular_distance=NA_real_,signed_shift=NA_real_,status="UNMATCHED_TARGET")
  if(!length(rows)) return(data.frame(reference_x=numeric(),target_x=numeric(),circular_distance=numeric(),signed_shift=numeric(),status=character()))
  do.call(rbind,rows)
}
