# Independent numerical rescue. All likelihood/state/penalty/chart math is shared.
p6r_close <- function(a,b,tol) is.numeric(a) && is.numeric(b) && length(a)==length(b) &&
  all(is.finite(a)) && all(is.finite(b)) && all(abs(a-b)<=tol*(1+pmax(abs(a),abs(b))))
p6r_settings <- function(cfg) {
  r <- cfg$rescue
  for(n in c("seed","n_local")) integer_scalar(r[[n]],paste("rescue",n),0L)
  integer_scalar(r$grid_n,"projection grid",52L)
  if(r$grid_n%%52L!=0) stop("Projection grid_n must be a multiple of 52")
  integer_scalar(r$max_cycles,"block cycles",1L)
  if(!is.logical(r$blockwise) || length(r$blockwise)!=1L || is.na(r$blockwise)) stop("Declare rescue blockwise TRUE/FALSE")
  for(n in c("local_sd","absolute_stop","relative_stop","boundary_distance","gradient_agreement",
    "gradient_tolerance","agreement_objective","agreement_coordinates")) scalar_number(r[[n]],paste("rescue",n),0)
  scalar_number(r$audit_tolerance,"audit tolerance",0,strict_lower=TRUE)
  finite_numeric(r$gradient_steps,"gradient steps")
  if(length(r$gradient_steps)<2 || any(r$gradient_steps<=0) || anyDuplicated(r$gradient_steps)) stop("At least two distinct positive gradient steps required")
  r
}
p6r_validate <- function(cfg,source,data) {
  if(!identical(source$type,"phase6_fit") || !identical(source$config$independent,TRUE)) stop("Independent phase6_fit source required")
  p6r_settings(cfg)
  # Only numerical controls, starts, diagnostics, campaign name may differ.
  for(n in c("data","expected_window","fixed","fixed_sources","independent"))
    if(!identical(cfg[[n]],source$config[[n]])) stop("Rescue scientific assumption mismatch: ",n)
  for(n in c("free","joint","bounds","coordinate_bounds"))
    if(!identical(cfg$estimation[[n]],source$config$estimation[[n]])) stop("Rescue support/free-policy mismatch: ",n)
  for(n in c("cases","seasonal_week","npieff","likelihood_include","fixed_parameters","time_markers"))
    if(!identical(data[[n]],source$data[[n]])) stop("Rescue observation/window/population mismatch: ",n)
  for(n in c("input_hashes","preprocessing_id","case_unit","manifest"))
    if(!identical(data$provenance[[n]],source$data$provenance[[n]])) stop("Rescue input/preprocessing provenance mismatch: ",n)
  if(!identical(source$observed,data$cases) || !identical(source$mask,data$likelihood_include)) stop("Saved observation/mask mismatch")
  for(name in c("AS","SPL")) {
    f <- source$models[[name]]$fit$representative
    if(is.null(f) || !isTRUE(f$valid) || !isTRUE(f$evaluation$valid)) stop("Valid saved ",name," endpoint required")
    p <- p6_problem(cfg,data,name)
    if(!identical(p6_model(cfg,name),p6_model(source$config,name)) || !identical(f$basis,p$basis) ||
      !identical(f$map$order,p$map$order) || !identical(f$map$lower,p$map$lower) || !identical(f$map$upper,p$map$upper) ||
      !identical(f$map,p$map) || !identical(f$lambda,p$model$lambda))
      stop("Saved model/basis/map/K/lambda mismatch: ",name)
    context <- f$comparison_context
    if(!identical(context$cases,data$cases) || !identical(context$mask,data$likelihood_include) ||
      !identical(context$weeks,data$seasonal_week) || !identical(context$npieff,data$npieff) ||
      !identical(context$fixed,data$fixed_parameters)) stop("Saved endpoint data contract mismatch: ",name)
    check_optimizer_coordinates(f$coordinates,p$map)
  }
  invisible(TRUE)
}
p6r_audit <- function(cfg,source,data) {
  p6r_validate(cfg,source,data);tol <- cfg$rescue$audit_tolerance
  rows <- lapply(c("AS","SPL"),function(name) {
    f <- source$models[[name]]$fit$representative;p <- p6_problem(cfg,data,name)
    e <- evaluate_candidate(f$coordinates,p);old <- f$evaluation
    if(!e$valid || !p6r_close(c(e$loglik,e$penalty,e$objective),c(old$loglik,old$penalty,old$objective),tol) ||
      !p6r_close(e$beta_by_week,old$beta_by_week,tol) ||
      !p6r_close(as.numeric(as.matrix(e$states)),as.numeric(as.matrix(old$states)),tol) ||
      !p6r_close(old$objective,-old$loglik+old$penalty,tol)) stop("Saved objective/state reconstruction mismatch: ",name)
    cached <- p6_physical(f$coordinates,p$map)
    physical <- old$parameters
    # Compare physical values directly; avoid inverse-logit roundoff at a guard endpoint.
    values <- c(beta0=physical$beta0,phi=physical$phi,S0=physical$S0,I0=physical$I0,rho=physical$rho,
      omega=physical$omega,delta=physical$delta,mu=physical$mu,gamma=physical$gamma,
      setNames(physical$free_coordinates,p$map$shape_order))
    if(!p6r_close(cached,values,tol)) stop("Saved decoded parameter mismatch: ",name)
    data.frame(model=name,saved_NLL=-old$loglik,reconstructed_NLL=-e$loglik,penalty=e$penalty,
      objective=e$objective,residual=e$objective-(-e$loglik+e$penalty),status="SAVED_REPLAY_MATCH")
  })
  do.call(rbind,rows)
}
p6r_project <- function(annual_semiannual,problem,grid_n=1040L) {
  if(!isTRUE(annual_semiannual$map$fourier) || annual_semiannual$basis$harmonics!=2L ||
    isTRUE(problem$map$fourier) || !inherits(problem$basis,"cyclic_cubic_basis")) stop("Projection requires saved AS and target SPL")
  integer_scalar(grid_n,"projection grid",52L)
  if(grid_n%%52L!=0) stop("Projection grid must include all 52 integer phases")
  x <- (0:(grid_n-1))*52/grid_n
  evaluate <- p5_beta_evaluator(annual_semiannual);target <- evaluate(x)
  parameters <- annual_semiannual$evaluation$parameters
  K <- problem$basis$K;C <- rbind(diag(K-1L),rep(-1,K-1L))
  B <- cyclic_basis_matrix(x,problem$basis)
  # Explicit intercept absorbs the target normalizer; c=Cz enforces sum(c)=0.
  design <- cbind(intercept=1,B%*%C)
  least_squares <- qr.solve(design,log(target/parameters$beta0))
  parameters$free_coordinates <- unname(least_squares[-1])
  coefficients <- constrained_coefficients(parameters$free_coordinates,problem$basis)
  reconstructed <- evaluate_spline_beta(x,coefficients,parameters$beta0,problem$basis)
  weekly_target <- evaluate(0:51)
  weekly_spline <- spline_beta_by_week(parameters$free_coordinates,parameters$beta0,problem$basis)
  if(!p6r_close(mean(weekly_target),parameters$beta0,1e-10) || !p6r_close(mean(weekly_spline),parameters$beta0,1e-10))
    stop("Fourier/SPL weekly mean normalization incompatible; beta0 cannot be copied")
  # The existing encoder enforces dimensions, support and the initial simplex.
  q <- tryCatch(encode_parameters(parameters,problem$map),error=function(e) e)
  evaluation <- if(inherits(q,"error")) NULL else evaluate_candidate(q,problem)
  row <- data.frame(beta_RMSE=sqrt(mean((reconstructed-target)^2)),beta_max_absolute_error=max(abs(reconstructed-target)),
    log_beta_RMSE=sqrt(mean((log(reconstructed)-log(target))^2)),log_beta_max_absolute_error=max(abs(log(reconstructed)-log(target))),
    AS_NLL=-annual_semiannual$evaluation$loglik,
    SPL_NLL=if(is.null(evaluation)) NA_real_ else -evaluation$loglik,
    SPL_penalty=if(is.null(evaluation)) NA_real_ else evaluation$penalty,
    SPL_objective=if(is.null(evaluation)) NA_real_ else evaluation$objective,
    AS_beta0=parameters$beta0,SPL_beta0=parameters$beta0,
    least_squares_intercept=least_squares[1],AS_log_normalizer=annual_semiannual$evaluation$log_normalizer,
    SPL_log_normalizer=weekly_log_normalizer(coefficients,problem$basis),
    normalization_intercept_difference= -weekly_log_normalizer(coefficients,problem$basis)-least_squares[1],
    status=if(inherits(q,"error") || !isTRUE(evaluation$valid)) "INVALID_WARM_START" else "VALID_APPROXIMATION_NOT_EXACT_NESTING")
  list(parameters=parameters,coordinates=if(inherits(q,"error")) NULL else q,
    reason=if(inherits(q,"error")) conditionMessage(q) else evaluation$failure_reason,
    report=row,curve=data.frame(x=x,week_coordinate=x+1,AS_beta=target,SPL_beta=reconstructed),
    evaluation=evaluation,coefficients=coefficients)
}
p6r_step <- function(problem,q,indices,optimizer=stats::optim,tolerance=1e-7) {
  start <- evaluate_candidate(q,problem);tracker <- make_objective(problem)
  answer <- list(initial_coordinates=q,initial_evaluation=start,indices=indices,accepted=FALSE,
    coordinates=q,evaluation=start,convergence=NA_integer_,warnings=character(),message="",status="INVALID")
  if(!start$valid) {answer$message <- start$failure_reason;return(answer)}
  fn <- function(v) {full <- q;full[indices] <- v;tracker$fn(full)}
  control <- c(problem$estimation$control,list(ndeps=problem$ndeps[indices],parscale=problem$parscale[indices]))
  args <- list(par=q[indices],fn=fn,method=problem$estimation$method,control=control,hessian=FALSE)
  if(problem$estimation$method=="L-BFGS-B") {args$lower <- problem$map$lower[indices];args$upper <- problem$map$upper[indices]}
  warnings <- character()
  opt <- tryCatch(withCallingHandlers(do.call(optimizer,args),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e) e)
  answer$warnings <- warnings;answer$failure_diagnostics <- tracker$diagnostics();answer$effective_control <- control
  if(inherits(opt,"error")) {answer$message <- conditionMessage(opt);return(answer)}
  if(!is.list(opt)) {answer$message <- "Optimizer returned a non-list result";return(answer)}
  answer$optimizer_return <- opt;answer$convergence <- opt$convergence;answer$message <- opt$message
  if(!is.numeric(opt$par) || length(opt$par)!=length(indices) ||
    !is.numeric(opt$convergence) || length(opt$convergence)!=1L || !is.finite(opt$convergence)) {
    answer$message <- "Malformed optimizer return";return(answer)
  }
  candidate <- q;candidate[indices] <- opt$par
  final <- evaluate_candidate(candidate,problem)
  answer$candidate_coordinates <- candidate;answer$candidate_evaluation <- final
  reported_ok <- is.numeric(opt$value) && length(opt$value)==1L && is.finite(opt$value) && opt$value<problem$estimation$failure_value
  answer$reported_objective_matches <- reported_ok && isTRUE(final$valid) && p6r_close(opt$value,final$objective,tolerance)
  reported_ok <- isTRUE(answer$reported_objective_matches)
  if(!reported_ok) answer$message <- paste(c(opt$message,"Invalid/sentinel/inconsistent reported objective; candidate retained in trace"),collapse="; ")
  answer$accepted <- isTRUE(final$valid) && reported_ok && final$objective<=start$objective
  answer$status <- if(!isTRUE(final$valid) || !reported_ok) "INVALID" else if(opt$convergence!=0) "NONCONVERGED" else "OPTIMIZER_TERMINATED"
  if(answer$accepted) {answer$coordinates <- candidate;answer$evaluation <- final}
  answer
}
p6r_gradient <- function(q,problem,settings,evaluator=evaluate_candidate) {
  # Keep the original ALL-stencils-valid certification rule; retain partial evidence.
  records <- lapply(settings$gradient_steps,function(h) {
    cells <- lapply(seq_along(q),function(i) {
      values <- list();components <- list();reasons <- character()
      for(side in c("plus","minus")) {
        v <- q;v[i] <- q[i]+if(side=="plus") h else -h
        if(v[i]<=problem$map$lower[i] || v[i]>=problem$map$upper[i]) {
          values[[side]] <- NA_real_;reasons <- c(reasons,paste(side,"stencil touches guard"))
        } else {
          e <- tryCatch(evaluator(v,problem),error=function(e) list(valid=FALSE,failure_reason=conditionMessage(e)))
          values[[side]] <- if(isTRUE(e$valid)) e$objective else NA_real_
          components[[side]] <- c(NLL=if(length(e$loglik)==1L) -e$loglik else NA_real_,
            penalty=if(length(e$penalty)==1L) e$penalty else NA_real_)
          if(!isTRUE(e$valid)) reasons <- c(reasons,paste(side,e$failure_stage,e$failure_reason))
        }
      }
      g <- if(length(reasons)) NA_real_ else (values$plus-values$minus)/(2*h)
      if(!length(reasons) && !is.finite(g)) reasons <- "Nonfinite central difference"
      component <- function(side,n) if(is.null(components[[side]])) NA_real_ else unname(components[[side]][n])
      data.frame(coordinate=problem$map$order[i],step=h,plus_objective=values$plus,minus_objective=values$minus,
        plus_NLL=component("plus","NLL"),minus_NLL=component("minus","NLL"),
        plus_penalty=component("plus","penalty"),minus_penalty=component("minus","penalty"),
        gradient=g,status=if(length(reasons) || !is.finite(g)) "UNAVAILABLE" else "COMPUTED",reason=paste(reasons,collapse="; "))
    })
    cells <- do.call(rbind,cells);complete <- all(cells$status=="COMPUTED")
    g <- setNames(cells$gradient,problem$map$order)
    list(status=if(complete) "COMPUTED" else "UNAVAILABLE",step=h,gradient=g,
      norm=if(complete) sqrt(sum(g^2)) else NA_real_,coordinates=cells,
      reason=if(complete) "" else paste(unique(cells$reason[cells$status!="COMPUTED"]),collapse="; "))
  })
  ok <- all(vapply(records,function(r) r$status=="COMPUTED",logical(1)))
  reliable <- FALSE;norm <- discrepancy <- NA_real_
  gs <- do.call(rbind,lapply(records,`[[`,"gradient"));last <- gs[nrow(gs),]
  if(ok) {
    discrepancy <- max(abs(sweep(gs,2,last)))/(1+max(abs(last)))
    reliable <- is.finite(discrepancy) && discrepancy<=settings$gradient_agreement
    if(reliable) norm <- sqrt(sum(last^2))
  }
  denominator <- if(all(is.finite(last))) 1+max(abs(last)) else NA_real_
  details <- do.call(rbind,lapply(records,function(r) {
    t <- r$coordinates;reference <- unname(last[t$coordinate])
    t$reference_step <- tail(settings$gradient_steps,1);t$reference_gradient <- reference
    t$absolute_disagreement <- abs(t$gradient-reference)
    t$coordinate_relative_disagreement <- t$absolute_disagreement/(1+abs(reference))
    t$normalized_disagreement <- t$absolute_disagreement/denominator
    t
  }))
  list(reliable=reliable,norm=norm,records=records,discrepancy=discrepancy,details=details)
}

p6r_run <- function(problem,start,settings,optimizer=stats::optim) {
  q <- start$coordinates
  if(is.null(q)) return(list(label=start$label,status="INVALID",reason=start$reason,trace=list(),initial_coordinates=NULL))
  initial <- evaluate_candidate(q,problem)
  if(!initial$valid) return(list(label=start$label,status="INVALID",reason=initial$failure_reason,trace=list(),initial_coordinates=q,initial_evaluation=initial))
  trace <- list();current <- q
  if(settings$blockwise) for(cycle in seq_len(settings$max_cycles)) {
    before <- evaluate_candidate(current,problem)$objective
    shape <- which(problem$map$order %in% c("log_beta0",problem$map$shape_order))
    for(block in list(shape,setdiff(seq_along(q),shape))) {
      step <- p6r_step(problem,current,block,optimizer,settings$audit_tolerance);step$cycle <- cycle;step$stage <- "block"
      trace[[length(trace)+1L]] <- step;current <- step$coordinates
    }
    after <- evaluate_candidate(current,problem)$objective
    if(before-after<=max(settings$absolute_stop,settings$relative_stop*abs(before))) break
  }
  full <- p6r_step(problem,current,seq_along(q),optimizer,settings$audit_tolerance);full$stage <- "full"
  trace[[length(trace)+1L]] <- full;current <- full$coordinates
  final <- evaluate_candidate(current,problem);gradient <- p6r_gradient(current,problem,settings)
  distance <- pmin(current-problem$map$lower,problem$map$upper-current)
  boundary <- names(current)[distance<=settings$boundary_distance]
  terminated <- identical(full$status,"OPTIMIZER_TERMINATED")
  stable <- terminated && full$accepted && !length(boundary) && gradient$reliable && gradient$norm<=settings$gradient_tolerance
  status <- if(length(boundary)) "BOUNDARY_SOLUTION" else if(stable) "NUMERICALLY_STABLE" else full$status
  dense <- tryCatch({
    x <- (0:(settings$grid_n-1))*52/settings$grid_n
    data.frame(x=x,beta=evaluate_spline_beta(x,final$spline_coefficients,final$parameters$beta0,problem$basis))
  },error=function(e) list(status="UNAVAILABLE",reason=conditionMessage(e)))
  list(label=start$label,status=status,optimizer_terminated=terminated,numerically_stable=stable,
    initial_coordinates=q,initial_evaluation=initial,coordinates=current,evaluation=final,
    gradient=gradient,boundary=boundary,boundary_distance=distance,trace=trace,
    full_step_accepted=full$accepted,parameter_record=parameter_record(problem,final$parameters),
    parameter_changes=p6_physical(current,problem$map)-p6_physical(q,problem$map),dense_beta=dense)
}
p6_rescue <- function(cfg,source,data,optimizer=stats::optim) {
  settings <- p6r_settings(cfg);audit <- p6r_audit(cfg,source,data)
  problem <- p6_problem(cfg,data,"SPL")
  warm <- tryCatch(p6r_project(source$models$AS$fit$representative,problem,settings$grid_n),error=function(e) {
    list(coordinates=NULL,reason=conditionMessage(e),evaluation=NULL,
      report=data.frame(beta_RMSE=NA_real_,beta_max_absolute_error=NA_real_,log_beta_RMSE=NA_real_,
        log_beta_max_absolute_error=NA_real_,AS_NLL= -source$models$AS$fit$representative$evaluation$loglik,
        SPL_NLL=NA_real_,SPL_penalty=NA_real_,SPL_objective=NA_real_,status="PROJECTION_FAILED",reason=conditionMessage(e)),
      curve=data.frame(x=numeric(),week_coordinate=numeric(),AS_beta=numeric(),SPL_beta=numeric()))
  })
  original <- source$models$SPL$fit$representative
  starts <- list(list(label="AS-projected",coordinates=warm$coordinates,reason=warm$reason),
    list(label="saved-SPL-endpoint",coordinates=original$coordinates))
  local <- with_phase2_seed(settings$seed,function() lapply(seq_len(settings$n_local),function(i) {
    if(is.null(warm$coordinates)) return(list(label=paste0("AS-local-",i),coordinates=NULL,reason="Projection unavailable"))
    list(label=paste0("AS-local-",i),coordinates=warm$coordinates+stats::rnorm(problem$map$dimension,sd=settings$local_sd))
  }))
  starts <- c(starts,local)
  runs <- lapply(starts,function(start) p6r_run(problem,start,settings,optimizer))
  valid <- which(vapply(runs,function(r) isTRUE(r$evaluation$valid),logical(1)))
  best <- if(length(valid)) valid[which.min(vapply(runs[valid],function(r) r$evaluation$objective,numeric(1)))] else NA_integer_
  list(type="phase6_rescue",status=if(is.na(best)) "NO_VALID_ENDPOINT" else "RESCUE_ENDPOINTS_STORED_NOT_CERTIFIED",
    source_hash=source$artifact_hash,source_identity=source$identity,source_config=source$config,
    config=cfg,model=problem$model,map=problem$map,basis=problem$basis,settings=settings,
    audit=audit,warm=warm,starts=starts,runs=runs,best_index=best,
    selected=if(is.na(best)) NULL else runs[[best]],source_AS=source$models$AS$fit$representative,
    source_SPL=original,selection_rule="Lowest valid accepted endpoint objective (including unchanged baselines); convergence/quality remain separate")
}
p6r_tables <- function(result) {
  scalar <- function(x) if(is.null(x) || length(x)!=1L) NA_real_ else x
  rows <- lapply(result$runs,function(r) {
    e <- r$evaluation;initial <- r$initial_evaluation;valid <- isTRUE(e$valid)
    trace <- r$trace
    invalid <- sum(vapply(trace,function(t) if(is.null(t$failure_diagnostics$invalid)) 0 else t$failure_diagnostics$invalid,numeric(1)))
    solver <- sum(vapply(trace,function(t) {
      reasons <- t$failure_diagnostics$reasons
      sum(reasons[grepl("^state_transition:",names(reasons))])
    },numeric(1)))
    last <- if(length(trace)) tail(trace,1)[[1]] else NULL
    absolute <- if(valid) initial$objective-e$objective else NA_real_
    data.frame(start=r$label,selected=isTRUE(r$label==result$selected$label),status=r$status,optimizer_terminated=isTRUE(r$optimizer_terminated),
      numerically_stable=isTRUE(r$numerically_stable),full_step_accepted=isTRUE(r$full_step_accepted),
      initial_objective=scalar(initial$objective),final_objective=scalar(e$objective),
      absolute_improvement=absolute,relative_improvement=if(valid) absolute/max(1,abs(initial$objective)) else NA_real_,
      NLL= -scalar(e$loglik),penalty=scalar(e$penalty),gradient_reliable=isTRUE(r$gradient$reliable),gradient_norm=scalar(r$gradient$norm),
      optimizer_code=scalar(last$convergence),boundary=paste(r$boundary,collapse=";"),
      min_coordinate_boundary_distance=if(is.null(r$boundary_distance)) NA_real_ else min(r$boundary_distance),
      invalid_objective_calls=invalid,SIRS_failure_calls=solver,
      weekly_beta_min=if(valid) min(e$beta_by_week) else NA_real_,weekly_beta_max=if(valid) max(e$beta_by_week) else NA_real_,
      dense_beta_min=if(is.data.frame(r$dense_beta)) min(r$dense_beta$beta) else NA_real_,
      dense_beta_max=if(is.data.frame(r$dense_beta)) max(r$dense_beta$beta) else NA_real_,
      S0=scalar(e$parameters$S0),I0=scalar(e$parameters$I0),R0=scalar(e$parameters$R0),rho=scalar(e$parameters$rho),
      coordinate_change_norm=if(valid) sqrt(sum((r$coordinates-r$initial_coordinates)^2)) else NA_real_,
      reason=paste(c(r$reason,last$message),collapse="; "),stringsAsFactors=FALSE)
  })
  table <- do.call(rbind,rows)
  pairs <- list()
  if(length(result$runs)>1) for(i in 2:length(result$runs)) for(j in seq_len(i-1L)) {
    a <- result$runs[[i]];b <- result$runs[[j]]
    if(!isTRUE(a$evaluation$valid) || !isTRUE(b$evaluation$valid)) next
    difference <- abs(a$evaluation$objective-b$evaluation$objective)
    coordinate_distance <- sqrt(sum((a$coordinates-b$coordinates)^2))
    pairs[[length(pairs)+1L]] <- data.frame(start_a=a$label,start_b=b$label,objective_difference=difference,
      coordinate_distance=coordinate_distance,same_objective=difference<=result$settings$agreement_objective,
      same_coordinates=coordinate_distance<=result$settings$agreement_coordinates,
      interpretation="Finite-start agreement only; parameter compensation and global optimality unresolved")
  }
  list(starts=table,agreement=p3_bind_rows(pairs))
}
p6_rescue_report <- function(result,out) {
  if(!identical(result$type,"phase6_rescue")) stop("Saved phase6_rescue result required")
  csv <- function(x,n) utils::write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE)
  tables <- p6r_tables(result);csv(tables$starts,"rescue_start_diagnostics");csv(tables$agreement,"repeated_start_agreement")
  csv(result$warm$report,"warm_start_approximation");csv(result$warm$curve,"warm_start_beta_approximation")
  csv(result$audit,"source_objective_audit")
  all <- list(source_AS=result$source_AS$evaluation,source_SPL=result$source_SPL$evaluation,
    AS_projected=result$warm$evaluation,rescued_SPL=result$selected$evaluation)
  comparison <- do.call(rbind,lapply(names(all),function(n) {
    e <- all[[n]]
    data.frame(endpoint=n,valid=isTRUE(e$valid),NLL=if(is.null(e)) NA_real_ else -e$loglik,
      penalty=if(is.null(e)) NA_real_ else e$penalty,objective=if(is.null(e)) NA_real_ else e$objective,
      definition="AS unpenalized reference; SPL objective includes exact campaign penalty; no model superiority claim")
  }))
  csv(comparison,"objective_comparison")
  parameters <- curves <- changes <- dense <- steps <- list()
  for(r in result$runs) if(isTRUE(r$evaluation$valid)) {
    p <- r$parameter_record;p$start <- r$label;p$selected <- isTRUE(r$label==result$selected$label);parameters[[r$label]] <- p
    curves[[r$label]] <- data.frame(start=r$label,selected=isTRUE(r$label==result$selected$label),week=1:52,beta=r$evaluation$beta_by_week)
    changes[[r$label]] <- data.frame(start=r$label,parameter=names(r$parameter_changes),change=unname(r$parameter_changes))
    if(is.data.frame(r$dense_beta)) {d <- r$dense_beta;d$start <- r$label;dense[[r$label]] <- d}
    for(i in seq_along(r$trace)) {
      t <- r$trace[[i]]
      steps[[length(steps)+1L]] <- data.frame(start=r$label,step=i,stage=t$stage,status=t$status,
        accepted=t$accepted,initial_objective=t$initial_evaluation$objective,
        retained_objective=t$evaluation$objective,candidate_objective=if(is.null(t$candidate_evaluation)) NA_real_ else t$candidate_evaluation$objective)
    }
  }
  csv(p3_bind_rows(parameters),"refined_parameter_estimates");csv(p3_bind_rows(curves),"refined_beta_curves")
  csv(p3_bind_rows(changes),"parameter_changes")
  csv(p3_bind_rows(dense),"refined_beta_dense");csv(p3_bind_rows(steps),"optimization_steps")
  csv(p6r_gradient_table(result$runs),"gradient_coordinates")
  writeLines(c("# Joint Spline-SIRS optimization rescue","",
    paste("Source fit hash:",result$source_hash),paste("K:",result$basis$K,"lambda:",result$model$lambda,
      "seed:",result$settings$seed,"blockwise:",result$settings$blockwise),
    paste("Selected start:",if(is.null(result$selected)) "NONE" else result$selected$label),result$selection_rule,"",
    "The original source artifact is unchanged. A rejected/failed full step retains the last valid baseline; inspect full_step_accepted and every trace.",
    "Code 0 means OPTIMIZER_TERMINATED, not scientific validity. NUMERICALLY_STABLE additionally requires accepted full step, interior guards, reproducible finite differences and small gradient.",
    "BOUNDARY_SOLUTION, NONCONVERGED and INVALID attempts remain explicit. A valid retained baseline can come from an unsuccessful optimizer attempt.",
    "Warm projection is constrained least squares, not exact AS-in-SPL nesting. Both beta0 values use the 52-week arithmetic mean, not a continuous mean.",
    "Inspect NLL, penalty and objective separately; penalized SPL objective need not beat AS NLL. No naive SPL AIC/BIC or AS-SPL chi-square LR is computed.",
    "Tiny beta0 compensated by huge omega, extreme initial fractions/rho, weak waning and different parameters with similar objectives warrant investigation; numerical stability does not establish practical identifiability.",
    "Finite starts, gradients and basis approximation cannot identify a unique cause, global optimum, causal mechanism or calibrated uncertainty. Compare blockwise on/off in separate new outputs."),file.path(out,"rescue_summary.md"))
  invisible(result)
}

# Cross-lambda analysis is deliberately separate from the SAME-lambda rescue contract.
p6l_settings <- function(cfg) {
  s <- cfg$lambda_sensitivity
  finite_numeric(s$values,"lambda grid",3L)
  if(!identical(unname(s$values),c(.1,1,10))) stop("Sensitivity values must be exactly 0.1, 1, 10 in that order")
  integer_scalar(s$seed,"lambda start seed",0L);integer_scalar(s$n_local,"lambda local starts",0L)
  scalar_number(s$local_sd,"lambda perturbation SD",0)
  scalar_number(s$min_prominence,"peak prominence",0)
  scalar_number(s$min_relative_prominence,"relative prominence",0,1)
  scalar_number(s$match_max_shift,"circular peak matching shift",0,26)
  p6r_settings(cfg);s
}
p6l_validate <- function(cfg,source,data,lambda) {
  p6l_settings(cfg)
  if(!lambda %in% c(.1,1,10)) stop("Unsupported target lambda")
  if(!identical(source$type,"phase6_rescue") || !isTRUE(source$config$independent) ||
    !isTRUE(source$selected$evaluation$valid)) stop("Valid saved phase6_rescue selected endpoint required")
  # This is NOT p6r_validate with a falsified source lambda. No source object is edited.
  for(n in c("data","expected_window","fixed","fixed_sources","independent","estimation","spline","rescue","peak_detection"))
    if(!identical(cfg[[n]],source$config[[n]])) stop("Cross-lambda source contract mismatch: ",n)
  if(!identical(as.integer(cfg$spline$K),8L) || source$basis$K!=8L) stop("Cross-lambda experiment requires K=8")
  base <- p6_problem(cfg,data,"SPL")
  if(!identical(source$model,base$model) || !identical(source$basis,base$basis) || !identical(source$map,base$map))
    stop("Source model/basis/map incompatible")
  context <- source$source_SPL$comparison_context
  for(pair in list(c("cases","cases"),c("weeks","seasonal_week"),c("npieff","npieff"),
    c("mask","likelihood_include"),c("markers","time_markers"),c("fixed","fixed_parameters")))
    if(!identical(context[[pair[1]]],data[[pair[2]]])) stop("Cross-lambda data/window mismatch: ",pair[1])
  # Rescue artifacts predate a top-level data field; original contexts and input identity
  # provide the observation and provenance contracts without inventing missing data.
  paths <- p3_inputs(cfg)
  if(length(paths)) {
    paths <- normalizePath(paths,mustWork=TRUE)
    hashes <- source$identity[paths]
    if(length(hashes)!=length(paths) || anyNA(hashes) || !identical(unname(hashes),unname(tools::md5sum(paths))))
      stop("Cross-lambda source input identity mismatch")
  }
  check_optimizer_coordinates(source$selected$coordinates,base$map)
  target <- cfg;target$spline$lambda <- lambda
  problem <- p6_problem(target,data,"SPL")
  expected <- base$model;expected$lambda <- lambda
  if(!identical(problem$model,expected) || !identical(problem$map,base$map) || !identical(problem$basis,base$basis))
    stop("Only lambda may differ in target model")
  invisible(problem)
}
p6_lambda_refit <- function(cfg,source,data,optimizer=stats::optim) {
  settings <- p6l_settings(cfg)
  for(lambda in settings$values) p6l_validate(cfg,source,data,lambda)
  base <- p6_problem(cfg,data,"SPL");q <- source$selected$coordinates
  replay <- evaluate_candidate(q,base);saved <- source$selected$evaluation;tol <- cfg$rescue$audit_tolerance
  if(!replay$valid || !p6r_close(c(replay$loglik,replay$penalty,replay$objective),c(saved$loglik,saved$penalty,saved$objective),tol) ||
    !p6r_close(unlist(replay$parameters),unlist(saved$parameters),tol) ||
    !p6r_close(replay$beta_by_week,saved$beta_by_week,tol) ||
    !p6r_close(as.numeric(as.matrix(replay$states)),as.numeric(as.matrix(saved$states)),tol)) stop("Selected source replay mismatch")
  starts <- c(list(list(label="saved-selected",coordinates=q)),with_phase2_seed(settings$seed,function()
    lapply(seq_len(settings$n_local),function(i) list(label=paste0("selected-local-",i),
      coordinates=q+stats::rnorm(length(q),sd=settings$local_sd)))))
  results <- lapply(settings$values,function(lambda) {
    problem <- p6l_validate(cfg,source,data,lambda)
    runs <- lapply(starts,function(s) p6r_run(problem,s,cfg$rescue,optimizer))
    valid <- which(vapply(runs,function(r) isTRUE(r$evaluation$valid),logical(1)))
    index <- if(length(valid)) valid[which.min(vapply(runs[valid],function(r) r$evaluation$objective,numeric(1)))] else NA_integer_
    selected <- if(is.na(index)) NULL else runs[[index]]
    peaks <- if(is.null(selected)) NULL else tryCatch(p5_beta_periodic_peaks(
      list(evaluation=selected$evaluation,basis=problem$basis,map=problem$map),cfg$peak_detection,
      settings$min_prominence,settings$min_relative_prominence),error=function(e) {
        empty <- data.frame(x=numeric(),week_coordinate=numeric(),beta=numeric(),prominence=numeric(),relative_prominence=numeric())
        list(status="PEAK_DIAGNOSTIC_FAILED",count=NA_integer_,table=empty,filtered=empty,reason=conditionMessage(e))
      })
    list(lambda=lambda,model=problem$model,map=problem$map,basis=problem$basis,runs=runs,selected=selected,
      best_index=index,settings=cfg$rescue,peaks=peaks,status=if(is.null(selected)) "NO_VALID_ENDPOINT" else "ENDPOINT_STORED_NOT_CERTIFIED")
  })
  names(results) <- as.character(settings$values)
  list(type="phase6_lambda_refit",status="LAMBDA_REFITS_STORED_NOT_CERTIFIED",config=cfg,
    source_hash=source$artifact_hash,source_identity=source$identity,source_config=source$config,
    source_coordinates=q,source_lambda=source$model$lambda,source_replay=replay,
    starts=starts,settings=settings,results=results,
    matching_reference="lambda=1; filtered periodic peaks, maximum-cardinality then minimum circular distance")
}
p6r_gradient_table <- function(runs) {
  p3_bind_rows(lapply(runs,function(r) {
    t <- r$gradient$details
    if(is.null(t)) return(NULL) # Older saved rescue artifacts remain reportable.
    t$start <- r$label;t
  }))
}
p6_lambda_report <- function(result,out) {
  if(!identical(result$type,"phase6_lambda_refit")) stop("Saved cross-lambda result required")
  csv <- function(t,name) utils::write.csv(t,file.path(out,paste0(name,".csv")),row.names=FALSE)
  summary <- diagnostics <- starts <- gradients <- peaks <- curves <- parameters <- list()
  for(r in result$results) {
    lambda <- r$lambda;selected <- r$selected;e <- selected$evaluation;available <- isTRUE(e$valid)
    peaks_available <- available && !identical(r$peaks$status,"PEAK_DIAGNOSTIC_FAILED")
    roughness <- if(available) cyclic_penalty(e$spline_coefficients,1) else NA_real_
    summary[[as.character(lambda)]] <- data.frame(lambda=lambda,K=8L,status=r$status,
      selected_start=if(available) selected$label else "",NLL=if(available) -e$loglik else NA_real_,
      roughness=roughness,penalty=if(available) e$penalty else NA_real_,objective=if(available) e$objective else NA_real_,
      raw_peak_count=if(peaks_available) r$peaks$count else NA_integer_,peak_count=if(peaks_available) nrow(r$peaks$filtered) else NA_integer_,
      peak_status=if(available) r$peaks$status else "UNAVAILABLE",convergence_status=if(available) selected$status else "UNAVAILABLE",
      peak_reason=if(is.null(r$peaks$reason)) "" else r$peaks$reason,
      starts=length(r$runs),optimizer_terminated=sum(vapply(r$runs,function(run) isTRUE(run$optimizer_terminated),logical(1))),
      numerically_stable=sum(vapply(r$runs,function(run) isTRUE(run$numerically_stable),logical(1))),
      gradient_reliable=if(available) selected$gradient$reliable else FALSE,
      gradient_norm=if(available) selected$gradient$norm else NA_real_,
      gradient_disagreement=if(available) selected$gradient$discrepancy else NA_real_,
      peak_locations=if(available) paste(r$peaks$filtered$week_coordinate,collapse=";") else "",
      peak_prominences=if(available) paste(r$peaks$filtered$prominence,collapse=";") else "")
    t <- p6r_tables(r)$starts;t$lambda <- lambda;diagnostics[[as.character(lambda)]] <- t
    t <- p6r_gradient_table(r$runs);t$lambda <- rep(lambda,nrow(t));gradients[[as.character(lambda)]] <- t
    for(run in r$runs) {
      initial <- run$initial_evaluation
      starts[[length(starts)+1L]] <- data.frame(lambda=lambda,start=run$label,valid=isTRUE(initial$valid),
        NLL=if(isTRUE(initial$valid)) -initial$loglik else NA_real_,
        roughness=if(isTRUE(initial$valid)) cyclic_penalty(initial$spline_coefficients,1) else NA_real_,
        penalty=if(isTRUE(initial$valid)) initial$penalty else NA_real_,objective=if(isTRUE(initial$valid)) initial$objective else NA_real_,
        coordinates=paste(run$initial_coordinates,collapse=";"))
    }
    if(available) {
      t <- r$peaks$table;t$lambda <- rep(lambda,nrow(t))
      t$retained <- t$prominence>=result$settings$min_prominence & t$relative_prominence>=result$settings$min_relative_prominence
      peaks[[as.character(lambda)]] <- t
      if(is.data.frame(selected$dense_beta)) {
        t <- selected$dense_beta;t$lambda <- rep(lambda,nrow(t));curves[[as.character(lambda)]] <- t
      }
      t <- selected$parameter_record;t$lambda <- lambda;parameters[[as.character(lambda)]] <- t
    }
  }
  matches <- list();reference <- result$results[["1"]]$peaks$filtered$x
  reference_available <- !is.null(result$results[["1"]]$selected) && result$results[["1"]]$peaks$status!="PEAK_DIAGNOSTIC_FAILED"
  for(r in result$results) if(r$lambda!=1) {
    if(!reference_available || is.null(r$selected) || r$peaks$status=="PEAK_DIAGNOSTIC_FAILED") {
      t <- data.frame(reference_x=NA_real_,target_x=NA_real_,circular_distance=NA_real_,signed_shift=NA_real_,status="MATCH_UNAVAILABLE")
    } else t <- p5_circular_peak_match(reference,r$peaks$filtered$x,result$settings$match_max_shift)
    t$lambda <- rep(r$lambda,nrow(t));matches[[as.character(r$lambda)]] <- t
  }
  csv(p3_bind_rows(summary),"lambda_comparison");csv(p3_bind_rows(starts),"starting_objectives")
  csv(p3_bind_rows(diagnostics),"lambda_start_diagnostics");csv(p3_bind_rows(gradients),"gradient_coordinates")
  csv(p3_bind_rows(peaks),"periodic_peaks");csv(p3_bind_rows(matches),"circular_peak_matches")
  beta <- p3_bind_rows(curves);csv(beta,"selected_beta_dense");csv(p3_bind_rows(parameters),"selected_parameters")
  # Same base-graphics device/line plotting pattern as saved Phase 5 reports.
  grDevices::png(file.path(out,"lambda_beta_curves.png"),width=1200,height=750,res=120)
  tryCatch({
    if(!"beta" %in% names(beta) || !nrow(beta)) {graphics::plot.new();graphics::title("No saved dense lambda curves available")}
    else {
      graphics::plot(range(beta$x+1),range(beta$beta),type="n",xlab="Seasonal week coordinate",ylab="beta (week^-1)")
      for(i in seq_along(curves)) graphics::lines(curves[[i]]$x+1,curves[[i]]$beta,col=i,lwd=2)
      graphics::legend("topright",paste("lambda",names(curves)),col=seq_along(curves),lty=1,bty="n")
    }
  },finally=grDevices::dev.off())
  writeLines(c("# Three-value spline smoothing sensitivity","",paste("Source rescue hash:",result$source_hash),
    "K=8; identical source-centered starts; every start re-evaluated at each lambda; no chained warm starts.",
    "Roughness=(1/2)*sum(cyclic second differences squared); penalty=lambda*roughness; objective=NLL+penalty.",
    "Objectives at different lambdas are different criteria, not model-selection evidence. No optimal lambda, ordinary SPL AIC/BIC or LR claim.",
    "Derivative-root peak detection is unchanged. Prominence is sampled periodic topographic prominence; equal highest peaks traverse a full turn.",
    result$matching_reference,"Unmatched peaks and empty peak sets are explicit; matching is descriptive, not proof of a persistent biological peak.",
    "NUMERICALLY_STABLE rules are unchanged. Any coordinate/step failure prevents complete-gradient certification; inspect per-cell reasons/disagreement.",
    "Finite starts, prominence mesh and matching tolerances require sensitivity review. No identifiability/global optimum/causal claim."),file.path(out,"lambda_summary.md"))
  invisible(result)
}
