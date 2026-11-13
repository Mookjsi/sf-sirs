# Diagnostic only: immutable endpoint coordinates; no optimizer, no production patch.
p7b_defaults <- function() list(repetitions=10L,
  steps=c(1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,1e-6),
  accuracy=c("production","stable_exits","stable_exits_nb2"),replay_tolerance=1e-7)
p7b_settings <- function(s) {
  integer_scalar(s$repetitions,"repeat count",2L);finite_numeric(s$steps,"diagnostic steps")
  if(length(s$steps)<2L || any(s$steps<=0) || any(diff(s$steps)>=0)) stop("At least two decreasing positive diagnostic steps required")
  if(!identical(s$accuracy,c("production","stable_exits","stable_exits_nb2"))) stop("Use the predefined three arithmetic settings")
  scalar_number(s$replay_tolerance,"replay tolerance",0,strict_lower=TRUE)
  s
}
p7b_engine <- function(accuracy) {
  if(!accuracy %in% c("production","stable_exits","stable_exits_nb2")) stop("Unknown diagnostic arithmetic setting")
  if(accuracy=="production") return(list(evaluate=evaluate_candidate,exit=function(a) 1-exp(-a),
    replacements=0L,description="Exact production weekly map and NB2; no solver tolerance exists"))
  # Replace ONLY the three cancellation-prone exit expressions in a PRIVATE function
  # clone. Fail closed on code drift; never assign into a production/global binding.
  count <- 0L
  rewrite <- function(expr) {
    if(is.call(expr) && identical(expr[[1]],as.name("-")) && length(expr)==3L &&
      identical(expr[[2]],1) && is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name("exp"))) {
      count <<- count+1L
      return(call("-",call("expm1",expr[[3]][[2]])))
    }
    if(is.call(expr)) for(i in seq_along(expr)[-1L]) expr[[i]] <- rewrite(expr[[i]])
    expr
  }
  solver <- simulate_weekly_sirs;body(solver) <- rewrite(body(solver))
  if(count!=3L) stop("Production solver structure drift: expected exactly three exit expressions")
  env <- new.env(parent=environment(evaluate_candidate));env$simulate_weekly_sirs <- solver
  if(accuracy=="stable_exits_nb2") {
    original <- nb2_loglik
    env$nb2_loglik <- function(y,mean,phi,include=NULL) {
      result <- original(y,mean,phi,include) # Reuse existing validation/schema/mask.
      logp <- stats::dnbinom(y,mu=mean,size=phi,log=TRUE)
      if(anyNA(logp) || any(logp==Inf)) stop("Diagnostic dnbinom arithmetic undefined")
      result$per_week <- logp;result$included_per_week <- logp
      result$included_per_week[!include] <- 0;result$total <- sum(result$included_per_week)
      result
    }
  }
  evaluator <- evaluate_candidate;environment(evaluator) <- env
  list(evaluate=evaluator,exit=function(a) -expm1(-a),replacements=count,
    description="Same discrete map with expm1 exit probabilities; optional equivalent R NB2 PMF arithmetic, NOT a timestep refinement")
}
p7b_sensitivity <- function(q,problem,evaluation,engine) {
  tryCatch({
    if(!isTRUE(evaluation$valid)) stop("Invalid center evaluation")
    if(isTRUE(problem$map$fourier) || !all(unlist(problem$map$free))) stop("Cross-check requires joint SPL simplex/logit chart")
    p <- decode_parameters(q,problem$map);d <- length(q);N <- problem$data$fixed_parameters$pop
    unit <- function(name) as.numeric(problem$map$order==name)
    rate <- function(name) p[[name]]*unit(paste0("log_",name))
    dmu <- rate("mu");dgamma <- rate("gamma");ddelta <- rate("delta");domega <- rate("omega")
    dphi <- rate("phi");drho <- p$rho*(1-p$rho)*unit("logit_rho")
    # Jacobian of the existing joint simplex chart, not another initial-state model.
    DS <- N*(p$S0*(1-p$S0)*unit("log_S_over_R")-p$S0*p$I0*unit("log_I_over_R"))
    DI <- N*(-p$S0*p$I0*unit("log_S_over_R")+p$I0*(1-p$I0)*unit("log_I_over_R"))
    DR <- -DS-DI
    K <- problem$basis$K;C <- vapply(seq_len(K-1L),function(i) constrained_coefficients(diag(K-1L)[,i],problem$basis),numeric(K))
    H <- cyclic_basis_matrix(0:51,problem$basis)%*%C;beta <- evaluation$beta_by_week
    dB <- matrix(0,52L,d);dB[,match("log_beta0",problem$map$order)] <- beta
    shape <- match(problem$map$shape_order,problem$map$order)
    dB[,shape] <- beta*sweep(H,2,colSums(H*(beta/sum(beta))))
    gradient_ll <- numeric(d);states <- evaluation$states
    for(i in 2:nrow(states)) {
      S <- states$S[i-1L];I <- states$I[i-1L];R <- states$R[i-1L]
      w <- problem$data$seasonal_week[i];npi <- problem$data$npieff[i]
      f <- beta[w]*(I+p$omega)/N*npi
      DF <- npi/N*(dB[w,]*(I+p$omega)+beta[w]*(DI+domega))
      # Differentiate existing competing-risk flows along the SHARED saved trajectory.
      a <- f+p$mu;b <- p$gamma+p$mu;c <- p$delta+p$mu
      outS <- engine$exit(a)*S;outI <- engine$exit(b)*I;outR <- engine$exit(c)*R
      DoutS <- engine$exit(a)*DS+S*exp(-a)*(DF+dmu)
      DoutI <- engine$exit(b)*DI+I*exp(-b)*(dgamma+dmu)
      DoutR <- engine$exit(c)*DR+R*exp(-c)*(ddelta+dmu)
      DSI <- (p$mu*DF-f*dmu)/a^2*outS+f/a*DoutS
      DIR <- (p$mu*dgamma-p$gamma*dmu)/b^2*outI+p$gamma/b*DoutI
      DRS <- (p$mu*ddelta-p$delta*dmu)/c^2*outR+p$delta/c*DoutR
      DM <- p$rho*DSI+states$incidence[i]*drho
      if(problem$data$likelihood_include[i]) {
        m <- states$expected_cases[i];y <- problem$data$cases[i]
        if(m<=0) stop("Zero scored mean: cross-check does not invent a boundary derivative")
        score_m <- p$phi*(y-m)/(m*(m+p$phi))
        score_phi <- digamma(y+p$phi)-digamma(p$phi)-log1p(m/p$phi)+(m-y)/(m+p$phi)
        gradient_ll <- gradient_ll+score_m*DM+score_phi*dphi
      }
      DS <- DS+DRS-DoutS+N*dmu;DI <- DI+DSI-DoutI;DR <- DR+DIR-DoutR
    }
    coefficients <- constrained_coefficients(p$free_coordinates,problem$basis)
    DC <- apply(C,2,cyclic_second_difference);gradient_penalty <- numeric(d)
    gradient_penalty[shape] <- problem$model$lambda*colSums(DC*cyclic_second_difference(coefficients))
    g <- -gradient_ll+gradient_penalty
    if(any(!is.finite(g))) stop("Nonfinite tangent cross-check")
    data.frame(coordinate=problem$map$order,status="COMPUTED",NLL_gradient=-gradient_ll,
      penalty_gradient=gradient_penalty,objective_gradient=g,reason="")
  },error=function(e) data.frame(coordinate=problem$map$order,status="UNAVAILABLE",NLL_gradient=NA_real_,
    penalty_gradient=NA_real_,objective_gradient=NA_real_,reason=conditionMessage(e)))
}
p7b_validate <- function(source,data,settings) {
  p7b_settings(settings)
  if(!identical(source$type,"phase6_lambda_refit") || !identical(names(source$results),c("0.1","1","10"))) stop("Saved three-value phase6_lambda_refit required")
  paths <- normalizePath(p3_inputs(source$config),mustWork=TRUE)
  if(length(paths) && (anyNA(source$identity[paths]) || !identical(unname(source$identity[paths]),unname(tools::md5sum(paths))))) stop("Saved input identity mismatch")
  for(r in source$results) {
    cfg <- source$config;cfg$spline$lambda <- r$lambda;p <- p6_problem(cfg,data,"SPL")
    if(!identical(p$model,r$model) || !identical(p$map,r$map) || !identical(p$basis,r$basis)) stop("Saved model/map/basis mismatch")
    if(is.null(r$selected)) next # Explicitly unavailable endpoint, not silently excluded.
    q <- r$selected$coordinates;check_optimizer_coordinates(q,p$map)
    e <- evaluate_candidate(q,p);old <- r$selected$evaluation
    if(!e$valid || !p6r_close(c(e$loglik,e$penalty,e$objective),c(old$loglik,old$penalty,old$objective),settings$replay_tolerance) ||
      !p6r_close(e$beta_by_week,old$beta_by_week,settings$replay_tolerance) ||
      !p6r_close(as.numeric(as.matrix(e$states)),as.numeric(as.matrix(old$states)),settings$replay_tolerance)) stop("Production endpoint replay mismatch")
  }
  invisible(TRUE)
}
p7b_diagnose <- function(source,data,settings=p7b_defaults(),engine_factory=p7b_engine) {
  p7b_validate(source,data,settings)
  repeats <- variability <- fd <- comparison <- sensitivities <- list();coordinates <- list()
  for(r in source$results) {
    if(is.null(r$selected)) {
      comparison[[length(comparison)+1L]] <- data.frame(lambda=r$lambda,accuracy="production",status="NO_SELECTED_ENDPOINT")
      next
    }
    q <- r$selected$coordinates;cfg <- source$config;cfg$spline$lambda <- r$lambda
    problem <- p6_problem(cfg,data,"SPL");production <- evaluate_candidate(q,problem)
    coordinates[[length(coordinates)+1L]] <- data.frame(lambda=r$lambda,coordinate=names(q),value=unname(q),
      lower=unname(problem$map$lower),upper=unname(problem$map$upper),
      saved_numerically_stable=isTRUE(r$selected$numerically_stable),saved_gradient_reliable=isTRUE(r$selected$gradient$reliable))
    production_fd <- NULL
    for(accuracy in settings$accuracy) {
      engine <- engine_factory(accuracy)
      evaluations <- lapply(seq_len(settings$repetitions),function(i) engine$evaluate(q,problem))
      table <- do.call(rbind,lapply(seq_along(evaluations),function(i) {
        e <- evaluations[[i]]
        data.frame(repetition=i,status=if(isTRUE(e$valid)) "VALID" else "INVALID",NLL=-e$loglik,
          penalty=e$penalty,objective=e$objective,reason=paste(c(e$failure_stage,e$failure_reason),collapse=": "))
      }))
      for(component in c("NLL","penalty","objective")) {
        good <- table$status=="VALID" & is.finite(table[[component]])
        baseline <- if(all(good)) table[[component]][1L] else NA_real_
        table[[paste0(component,"_absolute_deviation")]] <- abs(table[[component]]-baseline)
        table[[paste0(component,"_relative_deviation")]] <- abs(table[[component]]-baseline)/max(1,abs(baseline))
        spread <- if(all(good)) diff(range(table[[component]])) else NA_real_
        variability[[length(variability)+1L]] <- data.frame(lambda=r$lambda,accuracy=accuracy,component=component,
          absolute_range=spread,relative_range=spread/max(1,abs(baseline)),
          status=if(!all(good)) "INCONCLUSIVE_INVALID_REPEATS" else if(spread==0) "DETERMINISTIC_AT_TESTED_POINT" else "OBSERVED_REPEAT_VARIABILITY")
      }
      table$lambda <- r$lambda;table$accuracy <- accuracy;repeats[[length(repeats)+1L]] <- table
      center <- evaluations[[1L]]
      cross <- p7b_sensitivity(q,problem,center,engine);cross$lambda <- r$lambda;cross$accuracy <- accuracy
      sensitivities[[length(sensitivities)+1L]] <- cross
      local <- source$config$rescue;local$gradient_steps <- settings$steps
      gradient <- p6r_gradient(q,problem,local,evaluator=engine$evaluate)
      t <- gradient$details;t$lambda <- r$lambda;t$accuracy <- accuracy
      t$center_valid <- isTRUE(center$valid);t$center_objective <- center$objective
      t$coordinate_value <- unname(q[match(t$coordinate,names(q))]);t$lower <- unname(problem$map$lower[match(t$coordinate,names(q))])
      t$upper <- unname(problem$map$upper[match(t$coordinate,names(q))])
      t$NLL_gradient <- (t$plus_NLL-t$minus_NLL)/(2*t$step)
      t$penalty_gradient <- (t$plus_penalty-t$minus_penalty)/(2*t$step)
      t$tangent_gradient <- cross$objective_gradient[match(t$coordinate,cross$coordinate)]
      t$tangent_NLL_gradient <- cross$NLL_gradient[match(t$coordinate,cross$coordinate)]
      t$tangent_penalty_gradient <- cross$penalty_gradient[match(t$coordinate,cross$coordinate)]
      t$tangent_absolute_discrepancy <- abs(t$gradient-t$tangent_gradient)
      t$tangent_relative_discrepancy <- t$tangent_absolute_discrepancy/(1+abs(t$tangent_gradient))
      t$objective_difference <- t$plus_objective-t$minus_objective
      scale <- pmax(1,abs(t$plus_objective),abs(t$minus_objective))
      t$relative_objective_difference <- abs(t$objective_difference)/scale
      t$difference_in_epsilon_scale_units <- abs(t$objective_difference)/(.Machine$double.eps*scale)
      if(accuracy=="production") production_fd <- t
      t$production_gradient <- production_fd$gradient
      t$accuracy_gradient_difference <- t$gradient-t$production_gradient
      fd[[length(fd)+1L]] <- t
      comparison[[length(comparison)+1L]] <- data.frame(lambda=r$lambda,accuracy=accuracy,
        status=if(isTRUE(center$valid)) "EVALUATED_NOT_CERTIFIED" else "INCONCLUSIVE_INVALID_CENTER",
        NLL=-center$loglik,penalty=center$penalty,objective=center$objective,
        NLL_difference=-center$loglik+production$loglik,penalty_difference=center$penalty-production$penalty,
        objective_difference=center$objective-production$objective,
        max_state_absolute_difference=if(isTRUE(center$valid)) max(abs(as.matrix(center$states)-as.matrix(production$states))) else NA_real_,
        max_gradient_difference=if(all(is.finite(t$accuracy_gradient_difference))) max(abs(t$accuracy_gradient_difference)) else NA_real_,
        step_disagreement=gradient$discrepancy,exit_expression_replacements=engine$replacements)
    }
  }
  list(type="phase7b_gradient_diagnosis",status="DIAGNOSTICS_STORED_INFERENCE_INCONCLUSIVE",source_hash=source$artifact_hash,
    source_identity=source$identity,source_config=source$config,diagnostic_settings=settings,
    repeats=p3_bind_rows(repeats),variability=p3_bind_rows(variability),finite_differences=p3_bind_rows(fd),
    accuracy=p3_bind_rows(comparison),sensitivity=p3_bind_rows(sensitivities),coordinates=p3_bind_rows(coordinates),
    solver_contract="Fixed weekly discrete map; no rtol/atol/internal step. Arithmetic variants only; no temporal refinement or production default change")
}
p7b_report <- function(result,out) {
  if(!identical(result$type,"phase7b_gradient_diagnosis")) stop("Saved Phase 7B diagnosis required")
  csv <- function(t,n) utils::write.csv(t,file.path(out,paste0(n,".csv")),row.names=FALSE)
  csv(result$repeats,"objective_repeatability");csv(result$variability,"repeatability_summary")
  csv(result$finite_differences,"finite_difference_grid");csv(result$accuracy,"solver_accuracy_comparison")
  csv(result$sensitivity,"sensitivity_crosscheck");csv(result$coordinates,"endpoint_coordinates")
  grDevices::pdf(file.path(out,"gradient_step_diagnostics.pdf"),width=10,height=7)
  tryCatch({
    t <- result$finite_differences
    if(!all(c("lambda","coordinate") %in% names(t))) {graphics::plot.new();graphics::title("No selected endpoints available")}
    else for(lambda in unique(t$lambda)) for(coordinate in unique(t$coordinate)) {
      rows <- t[t$lambda==lambda & t$coordinate==coordinate,,drop=FALSE]
      good <- is.finite(rows$gradient)
      if(!any(good)) {graphics::plot.new();graphics::title(paste("Unavailable",lambda,coordinate));next}
      tangent <- rows$tangent_gradient[rows$accuracy=="production" & is.finite(rows$tangent_gradient)]
      graphics::plot(range(log10(rows$step)),range(c(rows$gradient[good],tangent)),type="n",xlab="log10 central-difference step",
        ylab="Objective gradient",main=paste("lambda",lambda,"coordinate",coordinate))
      for(i in seq_along(result$diagnostic_settings$accuracy)) {
        s <- rows[rows$accuracy==result$diagnostic_settings$accuracy[i],];s <- s[order(s$step),]
        graphics::lines(log10(s$step),s$gradient,col=i,lty=i,type="b")
      }
      if(length(tangent)) graphics::abline(h=tangent[1L],col=1,lty=3)
      graphics::legend("topright",c(result$diagnostic_settings$accuracy,if(length(tangent)) "production tangent"),
        col=c(1:3,if(length(tangent)) 1),lty=c(1:3,if(length(tangent)) 3),bty="n",cex=.8)
    }
  },finally=grDevices::dev.off())
  observed <- result$variability$status
  writeLines(c("# Numerical gradient reliability diagnosis","",paste("Source hash:",result$source_hash),result$solver_contract,"",
    paste("Repeatability statuses:",paste(unique(observed),collapse="; ")),
    "Exact repeated values only show determinism at tested points. They do not exclude deterministic cancellation, finite-difference truncation, roundoff or strong nonlinear sensitivity.",
    "Compare production vs stable_exits to isolate equivalent exit-probability arithmetic; compare stable_exits vs stable_exits_nb2 to isolate equivalent NB2 arithmetic. This is NOT solver time-step convergence.",
    "Compare per-coordinate step curves and NLL/penalty slopes with tangent sensitivities across MULTIPLE steps. A large coarse-step discrepancy with fine-step stabilization may reflect truncation; fine-step deterioration may reflect cancellation/roundoff. Neither is automatically diagnosed as the cause.",
    "Tangent derivatives use the mathematical discrete-map chain rule along shared states and analytic NB2/penalty scores. They have their own floating-point and conditioning limitations and are not certified truth.",
    "No selected coordinate, optimizer, production equation, basis, likelihood or saved certification flag was changed. No new certification is issued; optimizer code 0 is not gradient convergence.",
    "Inconclusive: no automatic mechanism ranking, accuracy certification, identifiability claim or empirical refit. Invalid stencils/cross-checks remain explicit. User interpretation of the diagnostic evidence is required.",
    "", "## Recorded repeatability evidence", "", "```",capture.output(print(result$variability,row.names=FALSE)),"```",
    "", "## Recorded arithmetic comparison", "", "```",capture.output(print(result$accuracy,row.names=FALSE)),"```"),
    file.path(out,"gradient_diagnosis_summary.md"))
  invisible(result)
}
