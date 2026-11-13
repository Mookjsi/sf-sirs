# One production-objective, production-tangent L-BFGS-B attempt per saved lambda.
p7c_defaults <- function() list(
  gate_steps=c(1e-3,3e-4,1e-4),gate_absolute=1e-4,gate_relative=5e-4,
  verification_steps=p7b_defaults()$steps,verification_absolute=1e-4,verification_relative=5e-4,
  replay_tolerance=1e-7,minimum_improvement=1e-8,optimizer_control=list(),
  persistent_weeks=c(32,51),peak_window=4)
p7c_settings <- function(s) {
  for(n in c("gate_steps","verification_steps")) {
    finite_numeric(s[[n]],n)
    if(length(s[[n]])<2L || any(s[[n]]<=0) || any(diff(s[[n]])>=0)) stop("Use multiple decreasing positive steps")
  }
  for(n in c("gate_absolute","gate_relative","verification_absolute","verification_relative","minimum_improvement")) scalar_number(s[[n]],n,0)
  scalar_number(s$replay_tolerance,"replay tolerance",0,strict_lower=TRUE)
  finite_numeric(s$persistent_weeks,"prespecified peak weeks",2L)
  if(any(s$persistent_weeks<1 | s$persistent_weeks>52)) stop("Peak weeks must be in 1..52")
  scalar_number(s$peak_window,"peak window",0,26)
  if(!is.list(s$optimizer_control) || (length(s$optimizer_control) &&
    (is.null(names(s$optimizer_control)) || anyDuplicated(names(s$optimizer_control)) ||
      any(!names(s$optimizer_control) %in% c("maxit","factr","pgtol","lmm"))))) stop("Only explicit L-BFGS-B numerical overrides allowed")
  for(n in names(s$optimizer_control)) scalar_number(s$optimizer_control[[n]],n,0)
  for(n in intersect(names(s$optimizer_control),c("maxit","lmm"))) integer_scalar(s$optimizer_control[[n]],n,1L)
  s
}
p7c_gate <- function(cfg,s,sensitivity=p7b_sensitivity) {
  # Explicit artificial observations and interior design points, not saved fit values.
  n <- 26L
  data <- adapt_preprocessed_island(data.frame(island="Honshu",year=2013L,week=seq_len(n),cases=rep(1:4,length.out=n)),
    data.frame(island="Honshu",pop=1e4),"Honshu")
  data$fixed_parameters$mu <- cfg$fixed$mu;data$fixed_parameters$gamma <- cfg$fixed$gamma
  engine <- p7b_engine("production");rows <- list()
  for(lambda in c(.1,1,10)) for(point in 1:2) {
    local <- cfg;local$spline$lambda <- lambda;p <- p6_problem(local,data,"SPL")
    record <- tryCatch({
      v <- list(beta0=if(point==1) 3 else 2.4,phi=if(point==1) 12 else 8,S0=if(point==1) .3 else .22,
        I0=if(point==1) .002 else .001,rho=if(point==1) .2 else .12,omega=if(point==1) 2 else 3,
        delta=if(point==1) 1/52 else .03,mu=if(cfg$estimation$joint$mu) 1/(80*52) else cfg$fixed$mu,
        gamma=if(cfg$estimation$joint$gamma) 1 else cfg$fixed$gamma,
        free_coordinates=if(point==1) c(.1,-.08,.04,-.03,.02,-.05,.01) else c(-.08,.12,-.03,.02,-.04,.01,.05))
      q <- encode_parameters(v,p$map);e <- evaluate_candidate(q,p);cross <- sensitivity(q,p,e,engine)
      if(!all(cross$status=="COMPUTED") || !identical(cross$coordinate,p$map$order)) stop("Synthetic tangent unavailable/order mismatch")
      g <- cross$objective_gradient;d <- length(q)
      directions <- cbind(diag(d),rep(1,d)/sqrt(d),(-1)^seq_len(d)/sqrt(d),sin(seq_len(d))/sqrt(sum(sin(seq_len(d))^2)))
      labels <- c(p$map$order,"direction_constant","direction_alternating","direction_sine")
      out <- list()
      for(i in seq_len(ncol(directions))) for(h in s$gate_steps) {
        plus <- evaluate_candidate(q+h*directions[,i],p);minus <- evaluate_candidate(q-h*directions[,i],p)
        analytic <- sum(g*directions[,i])
        fd <- if(plus$valid && minus$valid) (plus$objective-minus$objective)/(2*h) else NA_real_
        error <- abs(fd-analytic);limit <- s$gate_absolute+s$gate_relative*abs(analytic)
        out[[length(out)+1L]] <- data.frame(lambda=lambda,point=point,target=labels[i],step=h,
          analytic=analytic,finite_difference=fd,error=error,tolerance=limit,passed=is.finite(error) && error<=limit,reason="")
      }
      do.call(rbind,out)
    },error=function(e) data.frame(lambda=lambda,point=point,target="GATE_ERROR",step=NA_real_,analytic=NA_real_,
      finite_difference=NA_real_,error=NA_real_,tolerance=NA_real_,passed=FALSE,reason=conditionMessage(e)))
    rows[[length(rows)+1L]] <- record
  }
  table <- p3_bind_rows(rows)
  list(passed=all(table$passed),table=table,
    criterion="ALL coordinate and three deterministic directional checks at ALL steps/points/lambdas pass abs+rel tolerance; no empirical optimization in gate")
}
p7c_metrics <- function(q,e,g,problem,boundary_distance=problem$estimation$boundary_tolerance) {
  valid <- !is.null(g) && all(is.finite(g))
  projected <- if(valid) q-pmin(problem$map$upper,pmax(problem$map$lower,q-g)) else rep(NA_real_,length(q))
  distance <- pmin(q-problem$map$lower,problem$map$upper-q)
  list(NLL=-e$loglik,penalty=e$penalty,objective=e$objective,
    gradient_norm=if(valid) sqrt(sum(g^2)) else NA_real_,gradient_max=if(valid) max(abs(g)) else NA_real_,
    projected_gradient_norm=if(valid) max(abs(projected)) else NA_real_,
    boundary=names(q)[distance<=boundary_distance])
}
p7c_optimize <- function(problem,q0,s,source_rules,optimizer=stats::optim) {
  engine <- p7b_engine("production")
  initial <- evaluate_candidate(q0,problem)
  counters <- c(objective_requests=0L,gradient_requests=0L,production_evaluations=0L,tangent_evaluations=0L)
  cache <- NULL;best <- list(coordinates=q0,evaluation=initial);failures <- list()
  get <- function(q) {
    if(is.null(cache) || !identical(q,cache$coordinates)) {
      counters["production_evaluations"] <<- counters["production_evaluations"]+1L
      e <- evaluate_candidate(q,problem);cache <<- list(coordinates=q,evaluation=e,gradient=NULL)
      if(e$valid && e$objective<best$evaluation$objective) best <<- list(coordinates=q,evaluation=e)
      if(!e$valid) failures[[length(failures)+1L]] <<- list(stage="objective",coordinates=q,reason=e$failure_reason)
    }
    cache
  }
  gradient <- function(q) {
    current <- get(q)
    if(!current$evaluation$valid) {
      failures[[length(failures)+1L]] <<- list(stage="gradient",coordinates=q,reason="Invalid objective; no fake gradient supplied")
      stop("Tangent requested at invalid objective; no fake gradient supplied")
    }
    if(is.null(current$gradient)) {
      counters["tangent_evaluations"] <<- counters["tangent_evaluations"]+1L
      t <- p7b_sensitivity(q,problem,current$evaluation,engine)
      if(!all(t$status=="COMPUTED") || !identical(t$coordinate,problem$map$order)) {
        failures[[length(failures)+1L]] <<- list(stage="gradient",coordinates=q,reason=paste(unique(t$reason),collapse="; "))
        stop("Production tangent unavailable/order mismatch")
      }
      cache$gradient <<- setNames(t$objective_gradient,problem$map$order)
    }
    cache$gradient
  }
  fn <- function(q) {counters["objective_requests"] <<- counters["objective_requests"]+1L;e <- get(q)$evaluation
    if(e$valid) e$objective else problem$estimation$failure_value}
  gr <- function(q) {counters["gradient_requests"] <<- counters["gradient_requests"]+1L;gradient(q)}
  initial_g <- tryCatch(gradient(q0),error=function(e) NULL)
  control <- utils::modifyList(problem$estimation$control,s$optimizer_control)
  # REPORT=1 + trace=1 expose actual L-BFGS-B iterations; these are logging controls.
  control$trace <- 1L;control$REPORT <- 1L;control$parscale <- problem$parscale;control$ndeps <- problem$ndeps
  warnings <- character();opt <- NULL
  trace <- capture.output(opt <- tryCatch(withCallingHandlers(optimizer(par=q0,fn=fn,gr=gr,method="L-BFGS-B",
    lower=problem$map$lower,upper=problem$map$upper,control=control),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    }),error=function(e) e))
  if(length(trace)) cat(paste(trace,collapse="\n"),"\n")
  if(inherits(opt,"error") || !is.list(opt)) failures[[length(failures)+1L]] <- list(stage="optimizer",
    coordinates=if(is.null(cache)) q0 else cache$coordinates,
    reason=if(inherits(opt,"error")) conditionMessage(opt) else "Malformed optimizer return")
  normal <- is.list(opt) && !inherits(opt,"error") && is.numeric(opt$convergence) && length(opt$convergence)==1L &&
    is.finite(opt$convergence) && opt$convergence==0L && is.numeric(opt$par) && length(opt$par)==length(q0)
  returned <- if(normal) evaluate_candidate(opt$par,problem) else NULL
  reported <- normal && isTRUE(returned$valid) && p6r_close(opt$value,returned$objective,s$replay_tolerance)
  if(normal && !reported) failures[[length(failures)+1L]] <- list(stage="optimizer_return",coordinates=opt$par,
    reason="Invalid endpoint or inconsistent reported objective")
  candidate <- if(reported) list(coordinates=opt$par,evaluation=returned) else best
  # Fresh call OUTSIDE cache is authoritative, for normal and abnormal candidates.
  independent <- evaluate_candidate(candidate$coordinates,problem)
  replay <- independent$valid && p6r_close(c(-candidate$evaluation$loglik,candidate$evaluation$penalty,candidate$evaluation$objective),
    c(-independent$loglik,independent$penalty,independent$objective),s$replay_tolerance)
  cross <- p7b_sensitivity(candidate$coordinates,problem,independent,engine)
  candidate_g <- if(all(cross$status=="COMPUTED") && identical(cross$coordinate,problem$map$order)) setNames(cross$objective_gradient,problem$map$order) else NULL
  rules <- source_rules;rules$gradient_steps <- s$verification_steps
  fd <- p6r_gradient(candidate$coordinates,problem,rules)
  t <- fd$details;t$tangent <- cross$objective_gradient[match(t$coordinate,cross$coordinate)]
  t$error <- abs(t$gradient-t$tangent);t$tolerance <- s$verification_absolute+s$verification_relative*abs(t$tangent)
  t$passed <- t$status=="COMPUTED" & is.finite(t$error) & t$error<=t$tolerance
  by_step <- vapply(s$verification_steps,function(h) all(t$passed[t$step==h]),logical(1))
  adjacent <- any(head(by_step,-1L) & tail(by_step,-1L))
  # Preserve the original source-grid reliability rule as an additional screen.
  original_fd <- p6r_gradient(candidate$coordinates,problem,source_rules)
  metrics0 <- p7c_metrics(q0,initial,initial_g,problem,source_rules$boundary_distance)
  metrics1 <- p7c_metrics(candidate$coordinates,independent,candidate_g,problem,source_rules$boundary_distance)
  improvement <- initial$objective-independent$objective
  accepted <- reported && replay && !is.null(candidate_g) && improvement>s$minimum_improvement
  supported <- accepted && !length(metrics1$boundary) && original_fd$reliable && adjacent &&
    original_fd$norm<=source_rules$gradient_tolerance &&
    metrics1$gradient_norm<=source_rules$gradient_tolerance && metrics1$projected_gradient_norm<=source_rules$gradient_tolerance
  iterations <- grep("^iter[[:space:]]+[0-9]+",trace,value=TRUE)
  count <- if(length(iterations)) max(as.integer(sub("^iter[[:space:]]+([0-9]+).*","\\1",iterations))) else NA_integer_
  list(initial_coordinates=q0,initial_evaluation=initial,initial_metrics=metrics0,
    candidate_coordinates=candidate$coordinates,candidate_evaluation=independent,candidate_metrics=metrics1,
    retained_coordinates=if(accepted) candidate$coordinates else q0,retained_evaluation=if(accepted) independent else initial,
    accepted=accepted,optimizer_terminated=normal,objective_improved=replay && improvement>s$minimum_improvement,
    numerical_convergence_supported=isTRUE(supported),independent_replay=replay,
    candidate_kind=if(reported) "OPTIMIZER_ENDPOINT" else "BEST_VALID_EVALUATED_DIAGNOSTIC_CANDIDATE",
    status=if(accepted) "IMPROVED_OPTIMIZER_ENDPOINT_ACCEPTED" else if(!reported) "ABNORMAL_DIAGNOSTIC_CANDIDATE_ONLY" else "BASELINE_RETAINED",
    optimizer_return=opt,returned_evaluation=returned,reported_objective_matches=reported,iterations=count,
    iteration_status=if(is.na(count)) "NOT_REPORTED_NOT_INFERRED_FROM_EVALUATIONS" else "L-BFGS-B_TRACE",
    optimizer_trace=trace,effective_control=control,counters=counters,failures=failures,warnings=warnings,best_valid_evaluated=best,
    certification=t,source_grid_diagnostic=original_fd,sensitivity=cross,verification_adjacent_steps=adjacent,
    parameter_displacement=p6_physical(candidate$coordinates,problem$map)-p6_physical(q0,problem$map),
    initial_parameters=p6_physical(q0,problem$map),candidate_parameters=p6_physical(candidate$coordinates,problem$map),
    coordinate_displacement=sqrt(sum((candidate$coordinates-q0)^2)),objective_improvement=improvement)
}
p7c_refit <- function(source,data,s=p7c_defaults(),optimizer=stats::optim,gate_function=p7c_gate) {
  p7c_settings(s)
  validation <- p7b_defaults();validation$replay_tolerance <- s$replay_tolerance
  p7b_validate(source,data,validation)
  if(source$config$spline$K!=8L || source$config$estimation$method!="L-BFGS-B" ||
    any(vapply(source$results,function(r) is.null(r$selected),logical(1)))) stop("K=8, L-BFGS-B and three saved selected endpoints required")
  gate <- gate_function(source$config,s)
  base <- list(type="phase7c_tangent_refit",source_hash=source$artifact_hash,source_identity=source$identity,
    source_config=source$config,diagnostic_settings=s,gate=gate,records=list(),config=source$config)
  if(!isTRUE(gate$passed)) {base$status <- "GRADIENT_GATE_FAILED_NO_OPTIMIZATION";return(base)}
  for(r in source$results) {
    cfg <- source$config;cfg$spline$lambda <- r$lambda;problem <- p6_problem(cfg,data,"SPL")
    run <- p7c_optimize(problem,r$selected$coordinates,s,source$config$rescue,optimizer)
    peak_settings <- source$config$lambda_sensitivity
    detect <- function(e) p5_beta_periodic_peaks(list(evaluation=e,map=problem$map,basis=problem$basis),
      source$config$peak_detection,peak_settings$min_prominence,peak_settings$min_relative_prominence)
    old <- detect(run$initial_evaluation);new <- detect(run$candidate_evaluation)
    matching <- p5_circular_peak_match(old$filtered$x,new$filtered$x,peak_settings$match_max_shift)
    anchor_old <- p5_circular_peak_match(s$persistent_weeks-1,old$filtered$x,s$peak_window)
    anchor_new <- p5_circular_peak_match(s$persistent_weeks-1,new$filtered$x,s$peak_window)
    has_two <- function(t) sum(t$status=="MATCHED")==2L
    x <- (0:(source$config$rescue$grid_n-1L))*52/source$config$rescue$grid_n
    curve <- function(e,x) p5_beta_evaluator(list(evaluation=e,map=problem$map,basis=problem$basis))(x)
    curves <- rbind(data.frame(grid="weekly",x=0:51,original=curve(run$initial_evaluation,0:51),candidate=curve(run$candidate_evaluation,0:51)),
      data.frame(grid="dense",x=x,original=curve(run$initial_evaluation,x),candidate=curve(run$candidate_evaluation,x)))
    curves$absolute_difference <- abs(curves$candidate-curves$original)
    curves$relative_difference <- curves$absolute_difference/curves$original
    run$lambda <- r$lambda;run$original_peaks <- old;run$candidate_peaks <- new;run$peak_matches <- matching
    run$anchor_original <- anchor_old;run$anchor_candidate <- anchor_new
    run$two_prespecified_peaks_preserved <- has_two(anchor_old) && has_two(anchor_new)
    run$peak_structure_preserved <- nrow(old$filtered)==nrow(new$filtered) && all(matching$status=="MATCHED")
    run$curves <- curves;base$records[[as.character(r$lambda)]] <- run
  }
  base$status <- "THREE_SINGLE_ATTEMPTS_STORED_NOT_SCIENTIFIC_CERTIFICATION"
  base
}
p7c_report <- function(result,out) {
  if(!identical(result$type,"phase7c_tangent_refit")) stop("Saved Phase 7C result required")
  csv <- function(t,n) utils::write.csv(t,file.path(out,paste0(n,".csv")),row.names=FALSE)
  summaries <- changes <- peaks <- matching <- anchors <- curves <- certificates <- counts <- list()
  gate <- result$gate$table;gate$stage <- "synthetic_validation_gate";certificates[[1L]] <- gate
  for(r in result$records) {
    a <- r$initial_metrics;b <- r$candidate_metrics;opt <- r$optimizer_return
    scalar <- function(x) if(is.null(x) || length(x)!=1L) NA_real_ else x
    summaries[[length(summaries)+1L]] <- data.frame(lambda=r$lambda,status=r$status,candidate_kind=r$candidate_kind,
      optimizer_terminated=r$optimizer_terminated,optimizer_code=if(inherits(opt,"error") || !is.list(opt)) NA_real_ else scalar(opt$convergence),
      optimizer_message=if(inherits(opt,"error")) conditionMessage(opt) else if(!is.list(opt)) "Malformed optimizer return" else paste(opt$message,collapse="; "),
      objective_improved=r$objective_improved,numerical_convergence_supported=r$numerical_convergence_supported,
      numerical_peak_structure_preserved=r$peak_structure_preserved,two_peaks_near_32_51_preserved=r$two_prespecified_peaks_preserved,
      accepted_for_new_result=r$accepted,initial_NLL=a$NLL,final_NLL=b$NLL,initial_penalty=a$penalty,final_penalty=b$penalty,
      retained_objective=r$retained_evaluation$objective,source_grid_reliable=r$source_grid_diagnostic$reliable,
      source_grid_gradient_norm=r$source_grid_diagnostic$norm,multistep_tangent_agreement=r$verification_adjacent_steps,
      initial_objective=a$objective,final_objective=b$objective,objective_improvement=r$objective_improvement,
      initial_gradient_norm=a$gradient_norm,final_gradient_norm=b$gradient_norm,initial_max_gradient=a$gradient_max,
      final_max_gradient=b$gradient_max,initial_projected_gradient_norm=a$projected_gradient_norm,
      final_projected_gradient_norm=b$projected_gradient_norm,coordinate_displacement=r$coordinate_displacement,
      boundary=paste(b$boundary,collapse=";"),independent_replay=r$independent_replay,iterations=r$iterations,
      iteration_status=r$iteration_status,objective_requests=unname(r$counters["objective_requests"]),
      optimizer_function_count=if(is.list(opt)) scalar(opt$counts["function"]) else NA_real_,
      optimizer_gradient_count=if(is.list(opt)) scalar(opt$counts["gradient"]) else NA_real_,
      gradient_requests=unname(r$counters["gradient_requests"]),production_evaluations=unname(r$counters["production_evaluations"]),
      tangent_evaluations=unname(r$counters["tangent_evaluations"]),failure_count=length(r$failures))
    changes[[length(changes)+1L]] <- data.frame(lambda=r$lambda,parameter=names(r$parameter_displacement),
      original=unname(r$initial_parameters),candidate=unname(r$candidate_parameters),change=unname(r$parameter_displacement))
    for(endpoint in c("original","candidate")) {
      p <- if(endpoint=="original") r$original_peaks else r$candidate_peaks
      t <- p$table;t$lambda <- rep(r$lambda,nrow(t));t$endpoint <- rep(endpoint,nrow(t));peaks[[length(peaks)+1L]] <- t
      t <- if(endpoint=="original") r$anchor_original else r$anchor_candidate
      t$lambda <- rep(r$lambda,nrow(t));t$endpoint <- rep(endpoint,nrow(t));anchors[[length(anchors)+1L]] <- t
      counts[[length(counts)+1L]] <- data.frame(lambda=r$lambda,endpoint=endpoint,status=p$status,
        raw_peak_count=p$count,filtered_peak_count=nrow(p$filtered),
        locations=paste(p$filtered$week_coordinate,collapse=";"),prominences=paste(p$filtered$prominence,collapse=";"))
    }
    t <- r$peak_matches;t$lambda <- rep(r$lambda,nrow(t));matching[[length(matching)+1L]] <- t
    t <- r$curves;t$lambda <- r$lambda;t$candidate_kind <- r$candidate_kind;t$accepted <- r$accepted;curves[[length(curves)+1L]] <- t
    t <- r$certification;t$lambda <- r$lambda;t$stage <- "final_multistep_verification";certificates[[length(certificates)+1L]] <- t
    t <- r$source_grid_diagnostic$details;t$lambda <- r$lambda;t$stage <- "unchanged_source_grid_diagnostic";certificates[[length(certificates)+1L]] <- t
  }
  summary <- if(length(summaries)) p3_bind_rows(summaries) else data.frame(lambda=c(.1,1,10),
    status="NOT_RUN_GRADIENT_GATE_FAILED",optimizer_terminated=FALSE,objective_improved=FALSE,numerical_convergence_supported=FALSE)
  csv(summary,"reoptimization_summary");csv(p3_bind_rows(certificates),"gradient_certification")
  csv(p3_bind_rows(changes),"parameter_changes");csv(p3_bind_rows(peaks),"peak_comparison")
  csv(p3_bind_rows(counts),"peak_counts");csv(p3_bind_rows(matching),"circular_peak_matches");csv(p3_bind_rows(anchors),"prespecified_peak_matches")
  beta <- p3_bind_rows(curves);csv(beta,"beta_curve_comparison")
  grDevices::pdf(file.path(out,"reoptimization_diagnostics.pdf"),width=10,height=7)
  tryCatch({
    if(!length(result$records)) {graphics::plot.new();graphics::title("Gradient gate failed: no optimization attempted")}
    for(r in result$records) {
      t <- r$curves[r$curves$grid=="dense",]
      graphics::plot(t$x+1,t$original,type="l",ylim=range(c(t$original,t$candidate)),xlab="Seasonal week",ylab="beta",
        main=paste("lambda",r$lambda,r$candidate_kind));graphics::lines(t$x+1,t$candidate,col=2)
      graphics::legend("topright",c("original","candidate"),col=1:2,lty=1,bty="n")
      graphics::plot(t$x+1,t$candidate-t$original,type="l",xlab="Seasonal week",ylab="Candidate minus original beta")
      z <- r$certification
      for(coordinate in unique(z$coordinate)) {
        a <- z[z$coordinate==coordinate,];ok <- is.finite(a$gradient)
        if(!any(ok)) {graphics::plot.new();graphics::title(paste("FD unavailable",r$lambda,coordinate));next}
        graphics::plot(log10(a$step),a$gradient,type="b",ylim=range(c(a$gradient[ok],a$tangent[is.finite(a$tangent)])),
          xlab="log10 FD step",ylab="Objective gradient",main=paste("lambda",r$lambda,coordinate))
        if(any(is.finite(a$tangent))) graphics::abline(h=a$tangent[which(is.finite(a$tangent))[1L]],col=2,lty=2)
      }
    }
  },finally=grDevices::dev.off())
  writeLines(c("# Phase 7C controlled tangent-gradient reoptimization","",paste("Source hash:",result$source_hash),
    paste("Gate passed:",result$gate$passed),result$status,"",
    "Exactly one production-objective / production-tangent L-BFGS-B attempt per lambda after the synthetic gate. No arithmetic replacement, restart, block search or hyperparameter search.",
    "Acceptance in THIS NEW result requires normal code 0, finite matching independent production evaluation, a valid final tangent and objective improvement beyond the explicit minimum. Original artifacts are never replaced.",
    "An abnormal return retains the best valid objective-evaluated point as a DIAGNOSTIC CANDIDATE ONLY. Comparisons use that candidate and are labeled; retained coordinates remain the original unless accepted.",
    "Optimization terminated, objective improved, numerical convergence supported, and numerical peak structure preserved are separate fields. None implies biological seasonality or statistical identifiability.",
    "Numerical support additionally requires interior bounds, small raw/projected tangent norms, the UNCHANGED source-grid FD reliability rule and agreement at two adjacent verification steps. No agreement at a single step is promoted to convergence.",
    "The projected norm is ||q-proj_bounds(q-g)||_infinity (unit step). It does not override raw gradients or boundary warnings.",
    "The two prespecified peaks near weeks 32/51 use one-to-one circular matching within the declared window, for BOTH original and candidate curves. Extra/lost peaks are separately recorded; matching is descriptive, not causal.",
    "Iteration count comes from REPORT=1 trace, not function/gradient counts; missing counts are explicitly NA. Optimizer return, trace, warnings, cache counters and failures remain in RDS.",
    "Remaining limitations: synthetic gate is not an empirical gradient proof; tangent arithmetic and FD comparisons can be ill-conditioned. One local attempt is not a global optimization or uncertainty analysis.",
    "", "## Recorded outcomes", "", "```",capture.output(print(summary,row.names=FALSE)),"```"),file.path(out,"phase7c_summary.md"))
  invisible(result)
}
