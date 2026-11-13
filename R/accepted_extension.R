# The accepted 2013-2019 fit is read-only; all estimation is post-2019.
ae_reference <- function(saved,path) {
  x<-saved$result
  if(!identical(x$type,"phase7c_tangent_refit")||!isTRUE(x$gate$passed)||
     !identical(x$config,x$source_config)||!isTRUE(x$config$independent))
    stop("Accepted Phase 7C reference contract missing")
  r<-x$records[["1"]]
  if(is.null(r)||!identical(r$status,"IMPROVED_OPTIMIZER_ENDPOINT_ACCEPTED")||
     !isTRUE(r$accepted)||!isTRUE(r$optimizer_terminated)||
     !isTRUE(r$independent_replay)||
     !identical(r$retained_coordinates,r$candidate_coordinates)||
     !identical(r$candidate_kind,"OPTIMIZER_ENDPOINT"))
    stop("Lambda=1 endpoint is not the recorded accepted fit")
  e<-r$retained_evaluation
  if(!isTRUE(e$valid)||nrow(e$states)!=364L||length(e$beta_by_week)!=52L||
     !identical(x$config$data$region,"Honshu")||x$config$data$analysis_end!=364L||
     !identical(x$config$data$include_npi,FALSE)) stop("Accepted fit support mismatch")
  p<-e$parameters
  fixed<-c("beta0","phi","S0","I0","rho","omega","delta","mu","gamma")
  if(!all(fixed %in% names(p))) stop("Accepted parameter fields missing")
  if(!all(is.finite(unlist(p[fixed])))) stop("Accepted parameter values invalid")
  final<-unlist(e$states[nrow(e$states),c("S","I","R")],use.names=TRUE)
  names(final)<-c("S","I","R")
  list(source=x,record=r,evaluation=e,parameters=p,beta=e$beta_by_week,
    final_state=final,artifact_path=normalizePath(path,mustWork=TRUE),
    artifact_hash=unname(tools::md5sum(path)),producer_commit=saved$metadata$code_commit)
}
ae_manifest <- function(ref) {
  p<-ref$parameters;e<-ref$evaluation;x<-ref$source
  list(reference_manifest=data.frame(
    artifact=ref$artifact_path,md5=ref$artifact_hash,
    producer_commit=ref$producer_commit,source_hash=x$source_hash,
    accepted_status=ref$record$status,numerical_convergence_supported=
      isTRUE(ref$record$numerical_convergence_supported),
    years="2013-2019",rows=nrow(e$states),likelihood_rows=e$n_observations,
    first_row="2013-W1 initialization; excluded from likelihood",
    final_row="2019-W52 after transition 364",final_S=ref$final_state["S"],
    final_I=ref$final_state["I"],final_R=ref$final_state["R"],
    state_mass=sum(ref$final_state),population_denominator=
      sum(unlist(e$states[1L,c("S","I","R")],use.names=FALSE)),beta_period=52L,
    beta_normalization="discrete weekly arithmetic mean via cyclic spline log normalizer",
    observation_model="NB2 reported case counts; mean=rho*incidence; fixed phi",
    historical_week53="2015 week 53 merged into legacy week 52 before fitting",
    historical_external_multiplier="one throughout accepted 364-row no-NPI fit",
    limitation="accepted optimizer endpoint; numerical convergence support not certified; conditional model"),
    frozen_parameters=data.frame(parameter=names(p)[vapply(p,function(v) is.numeric(v)&&length(v)==1L,logical(1))],
      value=unlist(p[vapply(p,function(v) is.numeric(v)&&length(v)==1L,logical(1))]),
      status="FROZEN_FROM_ACCEPTED_ARTIFACT"),
    seasonal_coefficients=data.frame(coefficient=seq_along(e$spline_coefficients),
      value=e$spline_coefficients,normalization=e$log_normalizer),
    beta_by_week=data.frame(week=1:52,beta=ref$beta),
    reference_candidates=do.call(rbind,lapply(names(x$records),function(key) {
      r<-x$records[[key]]
      data.frame(lambda=as.numeric(key),status=r$status,accepted=isTRUE(r$accepted),
        selected_primary=identical(key,"1"),
        rationale=if(identical(key,"1")) "Phase 8 source explicitly uses accepted lambda=1" else "preserved alternative")
    })),
    final_state=data.frame(year=2019L,week=52L,historical_row=364L,
      S=ref$final_state["S"],I=ref$final_state["I"],R=ref$final_state["R"],
      state_mass=sum(ref$final_state)))
}
ae_replay <- function(ref,tolerance) {
  cfg<-ref$source$config
  data<-p6_data(cfg)
  if(length(data$cases)!=364L||!identical(data$likelihood_include,c(FALSE,rep(TRUE,363L)))||
     any(data$npieff!=1)) stop("Archived data or likelihood mask differs")
  problem<-p6_problem(cfg,data,"SPL")
  live<-evaluate_candidate(ref$record$retained_coordinates,problem)
  old<-ref$evaluation
  if(!isTRUE(live$valid)) stop("Historical replay is invalid")
  delta<-function(a,b) max(abs(as.numeric(a)-as.numeric(b)))
  checks<-data.frame(component=c("parameters","beta","states","incidence",
    "expected_cases","loglik","penalty","objective","final_state"),
    maximum_absolute_difference=c(
      delta(unlist(live$parameters),unlist(old$parameters)),
      delta(live$beta_by_week,old$beta_by_week),
      delta(as.matrix(live$states),as.matrix(old$states)),
      delta(live$states$incidence,old$states$incidence),
      delta(live$states$expected_cases,old$states$expected_cases),
      delta(live$loglik,old$loglik),delta(live$penalty,old$penalty),
      delta(live$objective,old$objective),
      delta(unlist(tail(live$states[c("S","I","R")],1L)),ref$final_state)),
    tolerance=tolerance)
  checks$passed<-checks$maximum_absolute_difference<=tolerance*
    pmax(1,c(max(abs(unlist(old$parameters))),max(abs(old$beta_by_week)),
      max(abs(as.matrix(old$states))),max(abs(old$states$incidence)),
      max(abs(old$states$expected_cases)),abs(old$loglik),abs(old$penalty),
      abs(old$objective),max(abs(ref$final_state))))
  list(status=if(all(checks$passed)) "HISTORICAL_REPLAY_PASSED" else "HISTORICAL_REPLAY_FAILED",
    checks=checks,live=live,data=data)
}
ae_calendar <- function(panel,settings,policy=settings$week53_policy) {
  z<-panel[panel$island=="Honshu"&panel$year %in% settings$future_years,,drop=FALSE]
  z<-z[order(z$chronological_index),,drop=FALSE]
  expected<-do.call(rbind,lapply(settings$future_years,function(y) data.frame(year=y,
    week=seq_len(settings$expected_weeks[[as.character(y)]]))))
  if(nrow(z)!=nrow(expected)||anyDuplicated(z[c("year","week")])||
     !identical(as.integer(z$year),as.integer(expected$year))||
     !identical(as.integer(z$week),as.integer(expected$week))||
     !identical(as.integer(z$chronological_index),seq.int(366L,365L+nrow(z)))||
     anyNA(z$complete)||!all(z$complete)||any(!is.finite(z$complete_reported_cases))||
     any(z$complete_reported_cases<0|z$complete_reported_cases!=floor(z$complete_reported_cases)))
    stop("Future case/calendar support invalid; no imputation")
  w<-as.integer(z$week)
  if(policy=="repeat_week52") w[w==53L]<-52L else if(policy=="map_to_week1")
    w[w==53L]<-1L else stop("Unknown future week-53 policy")
  data.frame(year=z$year,source_week=z$week,chronological_index=z$chronological_index,
    transition_after_historical=seq_len(nrow(z)),seasonal_week=w,
    observed=z$complete_reported_cases,
    week53_policy=ifelse(z$week==53L,policy,"source_week_1_to_52"))
}
ae_q <- function(calendar,model,log_q=numeric()) {
  if(model=="E0") return(rep(1,nrow(calendar)))
  if(model!="E1"||length(log_q)!=3L||any(!is.finite(log_q))) stop("Invalid extension model or q")
  q<-rep(1,nrow(calendar));q[calendar$year==2020L]<-exp(log_q[1]);
  q[calendar$year==2021L]<-exp(log_q[2]);q[calendar$year>=2022L]<-exp(log_q[3]);q
}
ae_continue <- function(ref,calendar,model,log_q=numeric(),future_external=1) {
  q<-ae_q(calendar,model,log_q)
  if(length(future_external)!=1L||!is.finite(future_external)||future_external<=0)
    stop("Future external multiplier must be explicitly positive")
  pop<-sum(unlist(ref$evaluation$states[1L,c("S","I","R")],use.names=FALSE))
  fixed<-list(pop=pop,mu=ref$parameters$mu,gamma=ref$parameters$gamma,
    omega=ref$parameters$omega,delta=ref$parameters$delta,rho=ref$parameters$rho)
  states<-simulate_weekly_sirs(ref$beta,ref$final_state,fixed,
    c(52L,calendar$seasonal_week),c(1,rep(future_external,nrow(calendar))*q),
    initial_mass=sum(ref$final_state))[-1L,,drop=FALSE]
  states$year<-calendar$year;states$source_week<-calendar$source_week
  states$chronological_index<-calendar$chronological_index
  states$seasonal_beta<-ref$beta[calendar$seasonal_week]
  states$q<-q;states$future_external_multiplier<-future_external
  states$effective_beta<-states$seasonal_beta*q*future_external
  states
}
ae_score <- function(ref,calendar,model,log_q=numeric(),future_external=1) {
  states<-ae_continue(ref,calendar,model,log_q,future_external)
  likelihood<-nb2_loglik(calendar$observed,states$expected_cases,ref$parameters$phi,
    rep(TRUE,nrow(calendar)))
  list(states=states,loglik=likelihood$total,per_week_loglik=likelihood$per_week,
    penalty=0,objective=-likelihood$total)
}
ae_fit <- function(ref,calendar,settings) {
  external<-settings$future_external_multiplier
  baseline<-tryCatch(ae_score(ref,calendar,"E0",future_external=external),error=function(e)e)
  E0<-if(inherits(baseline,"error")) list(status="FAILED",model="E0",message=conditionMessage(baseline)) else
    list(status="SCORED",model="E0",log_q=numeric(),score=baseline,
      attempts=data.frame(start=0L,status="FIXED_BASELINE",objective=baseline$objective))
  set.seed(settings$starts$seed)
  attempts<-lapply(seq_along(settings$starts$log_q),function(i) {
    start<-settings$starts$log_q[[i]]
    if(length(start)!=3L||any(!is.finite(start))||any(start<settings$q_bounds[1]|start>settings$q_bounds[2]))
      return(list(start=i,status="INVALID_START",message="q start outside declared bounds"))
    initial<-tryCatch(ae_score(ref,calendar,"E1",start,external)$objective,error=function(e)e)
    if(inherits(initial,"error")) return(list(start=i,status="INVALID_START",message=conditionMessage(initial)))
    fn<-function(p) {
      v<-ae_score(ref,calendar,"E1",p,external)$objective
      if(!is.finite(v)) stop("Invalid extension objective")
      v
    }
    opt<-tryCatch(stats::optim(start,fn,method="L-BFGS-B",
      lower=rep(settings$q_bounds[1],3L),upper=rep(settings$q_bounds[2],3L),
      control=settings$optimizer),error=function(e)e)
    if(inherits(opt,"error")) return(list(start=i,status="OPTIMIZER_ERROR",
      initial_objective=initial,message=conditionMessage(opt)))
    score<-tryCatch(ae_score(ref,calendar,"E1",opt$par,external),error=function(e)e)
    valid<-opt$convergence==0L&&!inherits(score,"error")&&is.finite(score$objective)
    list(start=i,status=if(valid) "NUMERICALLY_VALID" else "FAILED",
      initial=start,initial_objective=initial,coordinates=opt$par,
      convergence=opt$convergence,message=opt$message,counts=opt$counts,
      parameter_movement=sqrt(sum((opt$par-start)^2)),
      objective=if(valid) score$objective else NA_real_,score=if(valid) score else NULL)
  })
  valid<-which(vapply(attempts,function(a) identical(a$status,"NUMERICALLY_VALID"),logical(1)))
  best<-if(length(valid)) valid[which.min(vapply(attempts[valid],function(a) a$objective,numeric(1)))] else NA_integer_
  E1<-list(status=if(is.finite(best)) "NUMERICALLY_VALID" else "FAILED",model="E1",
    attempts=attempts,best=best,
    log_q=if(is.finite(best)) attempts[[best]]$coordinates else NULL,
    score=if(is.finite(best)) attempts[[best]]$score else NULL)
  if(is.finite(best)) {
    E1$near_best_starts<-sum(vapply(attempts[valid],function(a)
      a$objective-E1$score$objective<=settings$near_best_objective_gap,logical(1)))
    E1$minimum_log_bound_distance<-min(E1$log_q-settings$q_bounds[1],
      settings$q_bounds[2]-E1$log_q)
    replay<-ae_score(ref,calendar,"E1",E1$log_q,external)
    E1$independent_replay_difference<-replay$objective-E1$score$objective
    h<-settings$stationarity_step
    E1$objective_neighborhood<-do.call(rbind,lapply(1:3,function(j) do.call(rbind,
      lapply(c(-1,1),function(direction) {
        p<-E1$log_q;p[j]<-max(settings$q_bounds[1],min(settings$q_bounds[2],p[j]+direction*h))
        value<-tryCatch(ae_score(ref,calendar,"E1",p,external)$objective,error=function(e) NA_real_)
        data.frame(coordinate=paste0("log_q_",j),direction=direction,
          actual_step=p[j]-E1$log_q[j],valid=is.finite(value),
          objective_change=value-E1$score$objective)
      }))))
    E1$feasible_descent_found<-any(E1$objective_neighborhood$valid&
      E1$objective_neighborhood$objective_change< -1e-6)
    E1$stationarity_screen<-if(E1$feasible_descent_found)
      "FEASIBLE_DESCENT_FOUND" else "NO_DESCENT_AT_DECLARED_STEP_ONLY"
  }
  list(E0=E0,E1=E1)
}
ae_metrics <- function(fit,calendar,ref) {
  if(is.null(fit$score)) return(list(status=fit$status,annual=data.frame(),weekly=data.frame(),
    boundary=data.frame(),summary=data.frame(model=fit$model,status=fit$status,
      loglik=NA_real_,penalty=NA_real_,objective=NA_real_,MAE=NA_real_,RMSE=NA_real_,
      lag1_residual_correlation=NA_real_)))
  s<-fit$score$states;y<-calendar$observed;e<-s$expected_cases
  annual<-do.call(rbind,lapply(sort(unique(calendar$year)),function(year) {
    i<-which(calendar$year==year);a<-y[i];b<-e[i]
    timing<-function(v,q) if(sum(v)>0) calendar$source_week[i][which(cumsum(v)>=q*sum(v))[1]] else NA_integer_
    data.frame(model=fit$model,year=year,observed_total=sum(a),expected_total=sum(b),
      signed_burden_error=sum(b)-sum(a),MAE=mean(abs(a-b)),RMSE=sqrt(mean((a-b)^2)),
      observed_peak_weeks=paste(calendar$source_week[i][a==max(a)],collapse=";"),
      expected_peak_weeks=paste(calendar$source_week[i][b==max(b)],collapse=";"),
      observed_25=timing(a,.25),observed_50=timing(a,.5),observed_75=timing(a,.75),
      expected_25=timing(b,.25),expected_50=timing(b,.5),expected_75=timing(b,.75))
  }))
  weekly<-data.frame(model=fit$model,year=calendar$year,source_week=calendar$source_week,
    chronological_index=calendar$chronological_index,observed=y,
    incidence=s$incidence,expected_cases=e,residual=y-e,S=s$S,I=s$I,R=s$R,
    seasonal_beta=s$seasonal_beta,q=s$q,effective_beta=s$effective_beta,
    per_week_loglik=fit$score$per_week_loglik)
  boundary<-weekly[weekly$source_week==1L&weekly$year>=2020L,,drop=FALSE]
  boundary$prior_expected<-vapply(boundary$chronological_index,function(ix)
    if(ix==calendar$chronological_index[1]) tail(ref$evaluation$states$expected_cases,1L) else
      weekly$expected_cases[match(ix-1L,weekly$chronological_index)],numeric(1))
  boundary$expected_jump<-boundary$expected_cases-boundary$prior_expected
  boundary$prior_q<-vapply(boundary$chronological_index,function(ix)
    if(ix==calendar$chronological_index[1]) 1 else weekly$q[match(ix-1L,weekly$chronological_index)],numeric(1))
  boundary$imposed_q_step<-abs(log(boundary$q/boundary$prior_q))>1e-12
  list(status=fit$status,annual=annual,weekly=weekly,boundary=boundary,
    summary=data.frame(model=fit$model,status=fit$status,loglik=fit$score$loglik,
      penalty=fit$score$penalty,objective=fit$score$objective,
      MAE=mean(abs(y-e)),RMSE=sqrt(mean((y-e)^2)),
      lag1_residual_correlation=if(stats::sd(y-e)>0) stats::cor((y-e)[-1L],(y-e)[-length(y)]) else NA_real_))
}
ae_figures <- function(ref,historical_observed,comparison,annual,weekly,out) {
  grDevices::pdf(file.path(out,"frozen_history_and_extension.pdf"),width=11,height=6)
  old<-ref$evaluation$states$expected_cases
  h<-historical_observed$observed
  future<-weekly[weekly$model=="E1",,drop=FALSE]
  if(!nrow(future)) future<-weekly[weekly$model=="E0",,drop=FALSE]
  graphics::plot(seq_along(old),h,type="l",col="gray40",xlab="Model transition row",
    ylab="Reported cases",main="Accepted historical fit and conditional extension",
    xlim=c(1,length(old)+nrow(future)),ylim=range(c(h,old,future$observed,future$expected_cases),finite=TRUE))
  graphics::lines(seq_along(old),old,col="blue")
  if(nrow(future)) {
    ix<-length(old)+seq_len(nrow(future))
    graphics::lines(ix,future$observed,col="gray40")
    graphics::lines(ix,future$expected_cases,col="red")
  }
  graphics::abline(v=length(old)+.5,lty=2)
  grDevices::dev.off()
  if(nrow(future)) {
    grDevices::pdf(file.path(out,"future_observed_expected.pdf"),width=11,height=5)
    graphics::plot(future$chronological_index,future$observed,type="l",col="gray40",
      xlab="Source chronological week",ylab="Reported cases",
      main="Future observed and conditional expected cases")
    graphics::lines(future$chronological_index,future$expected_cases,col="red")
    graphics::abline(v=future$chronological_index[future$source_week==1L],lty=3)
    grDevices::dev.off()
    grDevices::pdf(file.path(out,"future_beta_q_states.pdf"),width=10,height=9)
    op<-graphics::par(mfrow=c(3,1))
    graphics::plot(future$chronological_index,future$seasonal_beta,type="l",
      xlab="Source chronological week",ylab="Accepted seasonal beta")
    graphics::plot(future$chronological_index,future$q,type="s",
      xlab="Source chronological week",ylab="Additional q")
    graphics::matplot(future$chronological_index,future[c("S","I","R")],type="l",lty=1,
      xlab="Source chronological week",ylab="Compartment counts")
    graphics::par(op);grDevices::dev.off()
    grDevices::pdf(file.path(out,"future_boundaries.pdf"),width=10,height=7)
    op<-graphics::par(mfrow=c(2,1))
    graphics::plot(future$chronological_index,future$expected_cases,type="l",
      xlab="Source chronological week",ylab="Expected reported cases")
    graphics::abline(v=future$chronological_index[future$source_week==1L],lty=3)
    graphics::plot(future$chronological_index,future$q,type="s",
      xlab="Source chronological week",ylab="Post-period q")
    graphics::abline(v=future$chronological_index[future$source_week==1L],lty=3)
    graphics::par(op);grDevices::dev.off()
  }
  if(nrow(annual)) {
    grDevices::pdf(file.path(out,"future_annual_burden_timing.pdf"),width=10,height=7)
    op<-graphics::par(mfrow=c(2,1))
    for(model in unique(annual$model)) {
      q<-annual[annual$model==model,,drop=FALSE]
      if(model==unique(annual$model)[1]) graphics::plot(q$year,q$observed_total,type="b",
        xlab="Source year",ylab="Reported cases",main="Future annual burden")
      graphics::lines(q$year,q$expected_total,type="b",col=if(model=="E0") "blue" else "red")
    }
    graphics::plot(annual$year,annual$expected_50-annual$observed_50,type="b",
      xlab="Source year",ylab="Median cumulative timing error (weeks)")
    graphics::par(op);grDevices::dev.off()
  }
  invisible(out)
}
