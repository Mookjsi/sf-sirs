# Phase 21 uses the accepted end state and the production weekly recurrence.
p21_registry <- function(cfg) {
  rows<-list(data.frame(id="E0",family="E0",spacing=NA_integer_,lambda=0),
    data.frame(id="E1",family="E1",spacing=NA_integer_,lambda=0))
  for(lambda in cfg$smoothing_grid) rows[[length(rows)+1L]]<-
    data.frame(id=paste0("E2-BLOCK-L",lambda),family="E2-BLOCK",
      spacing=cfg$block_width,lambda=lambda)
  for(spacing in cfg$spline_spacings) for(lambda in cfg$smoothing_grid)
    rows[[length(rows)+1L]]<-data.frame(id=paste0("E2-SPLINE-K",spacing,"-L",lambda),
      family="E2-SPLINE",spacing=spacing,lambda=lambda)
  do.call(rbind,rows)
}
p21_basis <- function(family,spacing,n=209L) {
  if(family=="E1") return(list(family=family,n=n,dimension=3L))
  if(family=="E0") return(list(family=family,n=n,dimension=0L))
  if(family=="E2-BLOCK") {
    block<-((seq_len(n)-1L)%/%spacing)+1L
    return(list(family=family,n=n,spacing=spacing,block=block,
      boundaries=data.frame(block=seq_len(max(block)),first=seq(1,n,by=spacing),
        last=pmin(seq(1,n,by=spacing)+spacing-1L,n)),dimension=max(block)))
  }
  if(family!="E2-SPLINE") stop("Unknown candidate family")
  knots<-seq.int(1L+spacing,n-1L,by=spacing)
  raw<-splines::bs(seq_len(n),degree=3,knots=knots,
    Boundary.knots=c(1,n),intercept=TRUE)
  center<-colMeans(raw)
  centered<-sweep(raw,2,center)
  rank<-qr(centered,tol=1e-10)
  if(rank$rank!=ncol(raw)-1L) stop("Centered spline basis rank mismatch")
  keep<-rank$pivot[seq_len(rank$rank)]
  shape<-centered[,keep,drop=FALSE]
  X<-cbind(intercept=1,shape)
  stops<-sort(unique(c(1,knots,n)))
  u<-c(0,1/3,2/3,1)
  bernstein<-cbind((1-u)^3,3*u*(1-u)^2,3*u^2*(1-u),u^3)
  control<-do.call(rbind,lapply(seq_len(length(stops)-1L),function(i) {
    x<-stops[i]+u*(stops[i+1L]-stops[i])
    b<-predict(raw,x)
    solve(bernstein,cbind(1,sweep(b,2,center)[,keep,drop=FALSE]))
  }))
  # Four-point Gauss-Legendre quadrature integrates squared cubic-spline curvature.
  nodes<-c(-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526)
  weights<-c(.3478548451374539,.6521451548625461,.6521451548625461,.3478548451374539)
  full_knots<-c(rep(1,4),knots,rep(n,4))
  D<-matrix(0,0,ncol(X));w<-numeric()
  for(i in seq_len(length(stops)-1L)) {
    width<-stops[i+1L]-stops[i]
    x<-(stops[i+1L]+stops[i])/2+nodes*width/2
    deriv<-splines::splineDesign(full_knots,x,ord=4,
      derivs=rep(2L,length(x)),outer.ok=TRUE)
    D<-rbind(D,cbind(0,deriv[,keep,drop=FALSE]))
    w<-c(w,weights*width/2)
  }
  P<-crossprod(D,D*w)
  P<-(P+t(P))/2
  scale<-sum(diag(P)[-1L])/(ncol(P)-1L)
  if(!is.finite(scale)||scale<=0) stop("Invalid roughness scale")
  list(family=family,n=n,spacing=spacing,knots=knots,boundary=c(1,n),
    raw_basis=raw,center=center,keep=keep,shape=shape,X=X,
    control=control,penalty_matrix=P/scale,penalty_scale=scale,
    quadrature_nodes=nodes,quadrature_weights=weights,dimension=ncol(X),
    coefficient_order=c("intercept",paste0("shape_",seq_len(ncol(shape)))))
}
p21_logq <- function(candidate,basis,theta,calendar) {
  if(candidate$family=="E0") return(rep(0,nrow(calendar)))
  if(length(theta)!=basis$dimension||any(!is.finite(theta))) stop("Parameter dimension mismatch")
  if(candidate$family=="E1") {
    g<-rep(theta[3L],nrow(calendar));g[calendar$year==2020L]<-theta[1L]
    g[calendar$year==2021L]<-theta[2L];return(g)
  }
  if(candidate$family=="E2-BLOCK") return(theta[basis$block])
  as.vector(basis$X%*%theta)
}
p21_feasible <- function(candidate,basis,theta,cfg,strict=FALSE) {
  if(candidate$family=="E0") return(TRUE)
  if(length(theta)!=basis$dimension||any(!is.finite(theta))) return(FALSE)
  values<-if(candidate$family=="E2-SPLINE") as.vector(basis$control%*%theta) else theta
  margin<-if(strict) 1e-8 else -1e-8
  all(values>cfg$q_log_bounds[1]+margin&values<cfg$q_log_bounds[2]-margin)
}
p21_penalty <- function(candidate,basis,theta) {
  if(candidate$family %in% c("E0","E1")) return(0)
  rough<-if(candidate$family=="E2-BLOCK")
    sum(diff(theta)^2)/(basis$dimension-1L) else
      as.numeric(crossprod(theta,basis$penalty_matrix%*%theta))
  candidate$lambda*max(0,rough)
}
p21_score <- function(ref,calendar,candidate,basis,theta,cfg,mask=rep(TRUE,nrow(calendar)),
                      forecast_origin=NULL) {
  if(length(mask)!=nrow(calendar)||anyNA(mask)||!is.logical(mask)||!any(mask))
    stop("Invalid training mask")
  if(!p21_feasible(candidate,basis,theta,cfg)) stop("Continuous q bound violated")
  g<-p21_logq(candidate,basis,theta,calendar)
  if(!is.null(forecast_origin)) {
    if(forecast_origin<1L||forecast_origin>=nrow(calendar)) stop("Invalid forecast origin")
    g[(forecast_origin+1L):nrow(calendar)]<-g[forecast_origin]
  }
  q<-exp(g)
  pop<-sum(unlist(ref$evaluation$states[1L,c("S","I","R")],use.names=FALSE))
  fixed<-list(pop=pop,mu=ref$parameters$mu,gamma=ref$parameters$gamma,
    omega=ref$parameters$omega,delta=ref$parameters$delta,rho=ref$parameters$rho)
  state<-simulate_weekly_sirs(ref$beta,ref$final_state,fixed,
    c(52L,calendar$seasonal_week),c(1,q*cfg$future_external_multiplier),
    initial_mass=sum(ref$final_state))[-1L,,drop=FALSE]
  state$year<-calendar$year;state$source_week<-calendar$source_week
  state$chronological_index<-calendar$chronological_index
  state$reference_beta<-ref$beta[calendar$seasonal_week]
  state$q<-q;state$effective_beta<-state$reference_beta*q*cfg$future_external_multiplier
  # Hidden counts are replaced before computing the training likelihood.
  train_y<-calendar$observed;train_y[!mask]<-0
  ll<-nb2_loglik(train_y,state$expected_cases,ref$parameters$phi,mask)
  penalty<-p21_penalty(candidate,basis,theta)
  list(states=state,loglik=ll$total,per_week_loglik=ll$included_per_week,
    penalty=penalty,objective=-ll$total+penalty,training_rows=sum(mask))
}
p21_starts <- function(candidate,basis,cfg,calendar,e1=NULL,allow_e1=FALSE,seed=cfg$optimizer$seed) {
  p<-basis$dimension
  if(!p) return(list(list(label="zero",theta=numeric())))
  zero<-rep(0,p);constant<-zero;constant[1]<-log(.7)
  if(candidate$family=="E1") constant[]<-log(.7)
  starts<-list(list(label="zero",theta=zero),list(label="constant",theta=constant))
  if(allow_e1&&!is.null(e1)&&length(e1)==3L) {
    target<-p21_logq(data.frame(family="E1"),list(dimension=3L),e1,calendar)
    theta<-if(candidate$family=="E1") e1 else if(candidate$family=="E2-BLOCK")
      vapply(seq_len(p),function(j) mean(target[basis$block==j]),numeric(1)) else
        as.vector(qr.solve(basis$X,target))
    while(!p21_feasible(candidate,basis,theta,cfg,strict=TRUE)) theta<-theta*.5
    starts[[length(starts)+1L]]<-list(label="training_E1_approximation",theta=theta)
  }
  set.seed(seed)
  for(i in seq_len(cfg$optimizer$perturbed_starts)) {
    theta<-zero+stats::rnorm(p,sd=cfg$optimizer$perturb_sd)
    while(!p21_feasible(candidate,basis,theta,cfg,strict=TRUE)) theta<-theta*.5
    starts[[length(starts)+1L]]<-list(label=paste0("seeded_",i),theta=theta)
  }
  starts
}
p21_numeric_gradient <- function(theta,fn,feasible,step,adaptive=FALSE) {
  gradient<-numeric(length(theta))
  for(j in seq_along(theta)) {
    h<-step*max(1,abs(theta[j]))
    plus<-minus<-theta
    plus[j]<-plus[j]+h;minus[j]<-minus[j]-h
    can_plus<-feasible(plus);can_minus<-feasible(minus)
    if(adaptive) while(!can_plus&&!can_minus&&h>1e-10) {
      h<-h/2
      plus[j]<-theta[j]+h;minus[j]<-theta[j]-h
      can_plus<-feasible(plus);can_minus<-feasible(minus)
    }
    if(can_plus&&can_minus) gradient[j]<-(fn(plus)-fn(minus))/(2*h)
    else if(can_plus) gradient[j]<-(fn(plus)-fn(theta))/h
    else if(can_minus) gradient[j]<-(fn(theta)-fn(minus))/h
    else stop("No feasible finite-difference step for parameter ",j)
  }
  if(any(!is.finite(gradient))) stop("Nonfinite numerical gradient")
  gradient
}
p21_neighborhood <- function(ref,calendar,candidate,basis,theta,cfg,mask,base) {
  directions<-c(lapply(seq_along(theta),function(i) {v<-rep(0,length(theta));v[i]<-1;v}),
    list(rep(1/sqrt(length(theta)),length(theta))))
  rows<-list()
  for(i in seq_along(directions)) for(step in cfg$diagnostic_steps) for(sign in c(-1,1)) {
      trial<-theta+sign*step*directions[[i]]
      feasible<-p21_feasible(candidate,basis,trial,cfg)
      value<-if(feasible) tryCatch(p21_score(ref,calendar,candidate,basis,trial,cfg,mask)$objective,
        error=function(e) NA_real_) else NA_real_
      rows[[length(rows)+1L]]<-data.frame(direction=if(i<=length(theta)) paste0("coordinate_",i) else "joint_all",
        step=sign*step,feasible=feasible,objective_change=value-base,
        one_sided_slope=(value-base)/(sign*step))
  }
  do.call(rbind,rows)
}
p21_fit <- function(ref,calendar,candidate,basis,cfg,mask=rep(TRUE,nrow(calendar)),
                    e1=NULL,allow_e1=FALSE,seed=cfg$optimizer$seed) {
  if(candidate$family=="E0") {
    score<-p21_score(ref,calendar,candidate,basis,numeric(),cfg,mask)
    return(list(id=candidate$id,status="SCORED",theta=numeric(),score=score,
      seed=seed,attempts=list(list(label="fixed",status="SCORED",start=numeric(),
        final=numeric(),final_objective=score$objective))))
  }
  starts<-p21_starts(candidate,basis,cfg,calendar,e1,allow_e1,seed)
  fn<-function(theta) p21_score(ref,calendar,candidate,basis,theta,cfg,mask)$objective
  attempts<-lapply(seq_along(starts),function(i) {
    start<-starts[[i]];initial<-tryCatch(fn(start$theta),error=function(e)e)
    if(inherits(initial,"error")) return(list(label=start$label,status="INVALID_START",
      seed=seed,start=start$theta,message=conditionMessage(initial)))
    opt<-tryCatch({
      if(candidate$family=="E2-SPLINE") {
        A<-basis$control;lo<-cfg$q_log_bounds[1];hi<-cfg$q_log_bounds[2]
        grad<-function(p) p21_numeric_gradient(p,fn,
          function(v) p21_feasible(candidate,basis,v,cfg,strict=TRUE),
          cfg$optimizer$gradient_step)
        stats::constrOptim(start$theta,fn,grad=grad,ui=rbind(A,-A),
          ci=c(rep(lo,nrow(A)),rep(-hi,nrow(A))),method="BFGS",
          control=list(maxit=cfg$optimizer$maxit,reltol=cfg$optimizer$reltol),
          outer.iterations=cfg$optimizer$outer_iterations,mu=cfg$optimizer$barrier_mu)
      } else stats::optim(start$theta,fn,method="L-BFGS-B",
        lower=rep(cfg$q_log_bounds[1],basis$dimension),
        upper=rep(cfg$q_log_bounds[2],basis$dimension),
        control=list(maxit=cfg$optimizer$maxit,factr=1e7,pgtol=0))
    },error=function(e)e)
    if(inherits(opt,"error")) return(list(label=start$label,status="OPTIMIZER_ERROR",
      seed=seed,start=start$theta,initial_objective=initial,message=conditionMessage(opt)))
    score<-tryCatch(p21_score(ref,calendar,candidate,basis,opt$par,cfg,mask),error=function(e)e)
    valid<-isTRUE(opt$convergence==0L)&&!inherits(score,"error")&&is.finite(score$objective)
    list(label=start$label,seed=seed,status=if(valid) "NUMERICALLY_VALID" else "FAILED",
      start=start$theta,final=opt$par,initial_objective=initial,
      final_objective=if(inherits(score,"error")) NA_real_ else score$objective,
      convergence=opt$convergence,message=if(inherits(score,"error"))
        conditionMessage(score) else if(is.null(opt$message)) "" else opt$message,
      counts=opt$counts,outer_iterations=opt$outer.iterations,
      movement=sqrt(sum((opt$par-start$theta)^2)),
      bound_distance=if(candidate$family=="E2-SPLINE")
        min(as.vector(basis$control%*%opt$par)-cfg$q_log_bounds[1],
          cfg$q_log_bounds[2]-as.vector(basis$control%*%opt$par)) else
        min(opt$par-cfg$q_log_bounds[1],cfg$q_log_bounds[2]-opt$par),
      score=if(valid) score else NULL)
  })
  good<-which(vapply(attempts,function(a) identical(a$status,"NUMERICALLY_VALID"),logical(1)))
  if(!length(good)) return(list(id=candidate$id,status="FAILED",seed=seed,attempts=attempts))
  winner<-good[which.min(vapply(attempts[good],function(a) a$final_objective,numeric(1)))]
  chosen<-attempts[[winner]];replay<-p21_score(ref,calendar,candidate,basis,chosen$final,cfg,mask)
  near<-sum(vapply(attempts[good],function(a)
    a$final_objective-chosen$final_objective<=cfg$near_best_gap,logical(1)))
  winner_logq<-p21_logq(candidate,basis,chosen$final,calendar)
  clusters<-do.call(rbind,lapply(good,function(i) {
    a<-attempts[[i]]
    distance<-max(abs(p21_logq(candidate,basis,a$final,calendar)-winner_logq))
    data.frame(start=i,objective_gap=a$final_objective-chosen$final_objective,
      maximum_logq_difference=distance,near_best=
        a$final_objective-chosen$final_objective<=cfg$near_best_gap,
      shape_cluster=if(distance<=cfg$shape_cluster_max_logq_difference)
        "WINNER_NEIGHBORHOOD" else "DISTINCT_SHAPE")
  }))
  neighborhood<-p21_neighborhood(ref,calendar,candidate,basis,chosen$final,cfg,mask,
    chosen$final_objective)
  list(id=candidate$id,status="NUMERICALLY_VALID",seed=seed,
    theta=chosen$final,score=chosen$score,
    selected_start=winner,attempts=attempts,near_best_starts=near,
    near_best_clusters=clusters,
    replay_difference=replay$objective-chosen$score$objective,
    neighborhood=neighborhood,feasible_descent=any(neighborhood$feasible&
      neighborhood$objective_change< -1e-6,na.rm=TRUE))
}
# Evaluation and selection keep validation score separate from full-data objective.
p21_candidate <- function(registry,i) as.list(registry[i,,drop=FALSE])
p21_score_rows <- function(observed,expected,phi,index) {
  ll<-nb2_loglik(observed[index],expected[index],phi,rep(TRUE,length(index)))
  data.frame(rows=length(index),log_score=ll$total,
    MAE=mean(abs(observed[index]-expected[index])),
    RMSE=sqrt(mean((observed[index]-expected[index])^2)),
    observed_total=sum(observed[index]),expected_total=sum(expected[index]),
    observed_peak=index[which.max(observed[index])],
    expected_peak=index[which.max(expected[index])])
}
p21_fit_table <- function(fits) {
  do.call(rbind,lapply(names(fits),function(id) {
    f<-fits[[id]]
    do.call(rbind,lapply(seq_along(f$attempts),function(i) {
      a<-f$attempts[[i]]
      data.frame(candidate=id,start=i,label=a$label,status=a$status,
        seed=if(is.null(a$seed)) NA_integer_ else a$seed,
        initial_objective=if(is.null(a$initial_objective)) NA_real_ else a$initial_objective,
        final_objective=if(is.null(a$final_objective)) NA_real_ else a$final_objective,
        convergence=if(is.null(a$convergence)) NA_integer_ else a$convergence,
        function_evaluations=if(is.null(a$counts)) NA_integer_ else unname(a$counts[1]),
        gradient_evaluations=if(is.null(a$counts)) NA_integer_ else unname(a$counts[2]),
        outer_iterations=if(is.null(a$outer_iterations)) NA_integer_ else a$outer_iterations,
        movement=if(is.null(a$movement)) NA_real_ else a$movement,
        bound_distance=if(is.null(a$bound_distance)) NA_real_ else a$bound_distance,
        start_parameters=if(is.null(a$start)) "" else paste(signif(a$start,10),collapse=";"),
        final_parameters=if(is.null(a$final)) "" else paste(signif(a$final,10),collapse=";"),
        message=if(is.null(a$message)) "" else a$message)
    }))
  }))
}
p21_fit_all <- function(ref,calendar,registry,bases,cfg,mask,baseline=NULL,
                        allow_saved=FALSE,seed_offset=0L) {
  fits<-vector("list",nrow(registry));names(fits)<-registry$id
  for(i in seq_len(nrow(registry))) {
    c<-p21_candidate(registry,i);basis<-bases[[c$id]]
    if(allow_saved&&c$family %in% c("E0","E1")) {
      old<-baseline[[c$family]]
      if(is.null(old$score)||!all(mask)) stop("Saved baseline unavailable")
      theta<-if(c$family=="E0") numeric() else old$log_q
      replay<-p21_score(ref,calendar,c,basis,theta,cfg,mask)
      if(abs(replay$objective-old$score$objective)>1e-6*
         max(1,abs(old$score$objective))) stop("Saved baseline objective mismatch")
      fits[[i]]<-list(id=c$id,status="REUSED_BASELINE",theta=theta,score=replay,
        attempts=list(list(label="saved_baseline",status="REUSED_BASELINE",
          start=theta,final=theta,final_objective=replay$objective,
          message="replayed from accepted-extension baseline")))
    } else {
      training_e1<-if(c$family %in% c("E2-BLOCK","E2-SPLINE")&&
        !is.null(fits[["E1"]]$theta)) fits[["E1"]]$theta else NULL
      fits[[i]]<-tryCatch(p21_fit(ref,calendar,c,basis,cfg,mask,
        e1=training_e1,allow_e1=!is.null(training_e1),
        seed=cfg$optimizer$seed+seed_offset+i),error=function(e)
          list(id=c$id,status="FAILED",message=conditionMessage(e),
            attempts=list(list(label="action",status="ERROR",message=conditionMessage(e)))))
    }
  }
  fits
}
p21_validate <- function(ref,calendar,registry,bases,cfg,origin=nrow(calendar)) {
  available<-which(vapply(cfg$validation_blocks,function(b) max(b)<=origin,logical(1)))
  if(!length(available)) stop("No declared held-out block before origin")
  rows<-list();fit_detail<-list()
  for(k in available) {
    block<-cfg$validation_blocks[[k]]
    mask<-seq_len(nrow(calendar))<=origin
    mask[block]<-FALSE
    fits<-p21_fit_all(ref,calendar,registry,bases,cfg,mask,
      seed_offset=1000L*k+origin)
    fit_detail[[paste0("block_",k)]]<-fits
    for(i in seq_len(nrow(registry))) {
      c<-p21_candidate(registry,i);f<-fits[[c$id]]
      if(is.null(f$score)) rows[[length(rows)+1L]]<-data.frame(
        candidate=c$id,block=k,first=min(block),last=max(block),
        status=f$status,rows=length(block),log_score=NA_real_,MAE=NA_real_,
        RMSE=NA_real_,observed_total=NA_real_,expected_total=NA_real_,
        observed_peak=NA_integer_,expected_peak=NA_integer_)
      else rows[[length(rows)+1L]]<-cbind(data.frame(candidate=c$id,block=k,
        first=min(block),last=max(block),status=f$status),
        p21_score_rows(calendar$observed,f$score$states$expected_cases,
          ref$parameters$phi,block))
    }
  }
  table<-do.call(rbind,rows)
  aggregate<-do.call(rbind,lapply(registry$id,function(id) {
    z<-table[table$candidate==id,,drop=FALSE]
    data.frame(candidate=id,completed=sum(is.finite(z$log_score)),
      required=length(available),heldout_rows=sum(z$rows),
      heldout_log_score=if(all(is.finite(z$log_score))) sum(z$log_score) else NA_real_,
      heldout_MAE=if(all(is.finite(z$MAE))) weighted.mean(z$MAE,z$rows) else NA_real_,
      heldout_RMSE=if(all(is.finite(z$RMSE))) sqrt(weighted.mean(z$RMSE^2,z$rows)) else NA_real_)
  }))
  list(blocks=table,aggregate=aggregate,fits=fit_detail,origin=origin)
}
p21_select <- function(aggregate,registry,cfg) {
  eligible<-aggregate[is.finite(aggregate$heldout_log_score)&
    aggregate$completed==aggregate$required,,drop=FALSE]
  if(!nrow(eligible)) return(list(candidate=NA_character_,status="NO_VALID_CANDIDATE",
    rule="all declared folds required"))
  best<-max(eligible$heldout_log_score)
  tied<-eligible[best-eligible$heldout_log_score<=
    cfg$selection$score_tie_per_week*eligible$heldout_rows,,drop=FALSE]
  family_order<-c("E0","E1","E2-BLOCK","E2-SPLINE")
  reg<-registry[match(tied$candidate,registry$id),,drop=FALSE]
  order_index<-order(match(reg$family,family_order),reg$spacing,reg$lambda,na.last=TRUE)
  list(candidate=tied$candidate[order_index[1]],status="VALIDATION_SELECTED",
    rule="maximum held-out NB log score; within declared per-week tie choose simpler family, coarser spacing, lower lambda",
    best_score=best,tie_count=nrow(tied))
}
p21_metrics <- function(ref,calendar,candidate,fit,cfg,period=rep("reconstruction",nrow(calendar))) {
  if(is.null(fit$score)) return(list(status=fit$status,weekly=data.frame(),
    annual=data.frame(),summary=data.frame(candidate=candidate$id,status=fit$status)))
  state<-fit$score$states;y<-calendar$observed;e<-state$expected_cases
  ll<-nb2_loglik(y,e,ref$parameters$phi,rep(TRUE,length(y)))
  weekly<-data.frame(candidate=candidate$id,year=calendar$year,week=calendar$source_week,
    index=calendar$transition_after_historical,period=period,observed=y,
    expected=e,incidence=state$incidence,S=state$S,I=state$I,R=state$R,
    reference_beta=state$reference_beta,q=state$q,effective_beta=state$effective_beta,
    residual=y-e,log_score=ll$per_week)
  previous_mass<-c(sum(ref$final_state),head(rowSums(state[c("S","I","R")]),-1L))
  pop<-sum(unlist(ref$evaluation$states[1L,c("S","I","R")],use.names=FALSE))
  weekly$accounting_error<-rowSums(state[c("S","I","R")])-previous_mass-
    ref$parameters$mu*pop+
    state$deaths_S+state$deaths_I+state$deaths_R
  weekly$nonnegative_state<-apply(state[c("S","I","R")],1,function(v) all(v>=0))
  annual<-do.call(rbind,lapply(sort(unique(calendar$year)),function(year) {
    z<-weekly[weekly$year==year,,drop=FALSE]
    timing<-function(v,p) if(sum(v)<=0) NA_integer_ else z$week[which(cumsum(v)>=p*sum(v))[1]]
    local<-function(v) {
      i<-which(v>=c(-Inf,head(v,-1L))&v>=c(tail(v,-1L),-Inf)&
        (v>c(-Inf,head(v,-1L))|v>c(tail(v,-1L),-Inf)))
      paste(paste0(z$week[i],":",signif(v[i],5)),collapse=";")
    }
    data.frame(candidate=candidate$id,year=year,observed_total=sum(z$observed),
      expected_total=sum(z$expected),relative_burden_error=
        if(sum(z$observed)>0) (sum(z$expected)-sum(z$observed))/sum(z$observed) else NA_real_,
      log_score=sum(z$log_score),MAE=mean(abs(z$residual)),RMSE=sqrt(mean(z$residual^2)),
      observed_global_peaks=paste(z$week[z$observed==max(z$observed)],collapse=";"),
      expected_global_peaks=paste(z$week[z$expected==max(z$expected)],collapse=";"),
      observed_local_peaks=local(z$observed),expected_local_peaks=local(z$expected),
      observed_25=timing(z$observed,.25),expected_25=timing(z$expected,.25),
      observed_50=timing(z$observed,.5),expected_50=timing(z$expected,.5),
      observed_75=timing(z$observed,.75),expected_75=timing(z$expected,.75))
  }))
  boundary<-weekly[weekly$week==1L,,drop=FALSE]
  prev<-c(1,weekly$q[match(boundary$index[-1L]-1L,weekly$index)])
  boundary$prior_q<-prev;boundary$q_change<-boundary$q-prev
  boundary$prior_beta<-c(tail(ref$beta,1L),
    weekly$effective_beta[match(boundary$index[-1L]-1L,weekly$index)])
  boundary$beta_change<-boundary$effective_beta-boundary$prior_beta
  boundary$prior_expected<-c(tail(ref$evaluation$states$expected_cases,1L),
    weekly$expected[match(boundary$index[-1L]-1L,weekly$index)])
  boundary$expected_change<-boundary$expected-boundary$prior_expected
  error_rank<-order(abs(weekly$residual),decreasing=TRUE)
  largest<-weekly[head(error_rank,10L),,drop=FALSE]
  ac<-function(k) if(length(y)>k&&stats::sd(weekly$residual)>0)
    stats::cor(weekly$residual[1:(length(y)-k)],weekly$residual[(k+1):length(y)]) else NA_real_
  summary<-data.frame(candidate=candidate$id,status=fit$status,
    loglik=ll$total,training_loglik=fit$score$loglik,penalty=fit$score$penalty,
    objective=fit$score$objective,MAE=mean(abs(y-e)),RMSE=sqrt(mean((y-e)^2)),
    residual_lag1=ac(1),residual_lag4=ac(4),
    maximum_accounting_error=max(abs(weekly$accounting_error)),
    all_states_nonnegative=all(weekly$nonnegative_state),
    minimum_q=min(weekly$q),maximum_q=max(weekly$q))
  list(status=fit$status,weekly=weekly,annual=annual,boundary=boundary,
    largest_errors=largest,summary=summary)
}
p21_adequacy <- function(annual,weekly,cfg) {
  if(!nrow(annual)) return(list(status="NO_ADEQUATE_MODEL",reason="no annual metrics"))
  burden<-all(abs(annual$relative_burden_error)<=cfg$selection$annual_relative_burden,
    na.rm=FALSE)
  timing<-all(abs(annual$expected_50-annual$observed_50)<=
    cfg$selection$cumulative_median_timing_weeks,na.rm=FALSE)
  low<-weekly$observed<=pmax(1,stats::median(weekly$observed))
  false_peak<-any(weekly$expected[low]>
    cfg$selection$large_false_peak_ratio*pmax(1,weekly$observed[low]))
  passed<-isTRUE(burden)&&isTRUE(timing)&&!false_peak
  list(status=if(passed) "ADEQUATE_UNDER_OPERATIONAL_TOLERANCES" else "NO_ADEQUATE_MODEL",
    reason=paste("burden",burden,"timing",timing,"false_large_prediction",false_peak))
}
p21_forecast <- function(ref,calendar,registry,bases,cfg) {
  output<-list();rows<-list();weekly_rows<-list()
  for(origin in cfg$forecast_origins) {
    tuning<-p21_validate(ref,calendar,registry,bases,cfg,origin)
    choice<-p21_select(tuning$aggregate,registry,cfg)
    if(is.na(choice$candidate)) {
      rows[[length(rows)+1L]]<-data.frame(origin=origin,status="NO_VALID_CANDIDATE",
        candidate=NA_character_,rows=nrow(calendar)-origin,log_score=NA_real_,
        MAE=NA_real_,RMSE=NA_real_,observed_total=NA_real_,expected_total=NA_real_,
        observed_peak=NA_integer_,expected_peak=NA_integer_)
      output[[as.character(origin)]]<-list(tuning=tuning,choice=choice)
      next
    }
    c<-p21_candidate(registry,match(choice$candidate,registry$id))
    mask<-seq_len(nrow(calendar))<=origin
    train_e1<-NULL
    if(c$family %in% c("E2-BLOCK","E2-SPLINE")) {
      e1c<-p21_candidate(registry,match("E1",registry$id))
      e1fit<-p21_fit(ref,calendar,e1c,bases[["E1"]],cfg,mask,
        seed=cfg$optimizer$seed+origin)
      if(!is.null(e1fit$theta)) train_e1<-e1fit$theta
    }
    fit<-tryCatch(p21_fit(ref,calendar,c,bases[[c$id]],cfg,mask,
      e1=train_e1,allow_e1=!is.null(train_e1),
      seed=cfg$optimizer$seed+origin+100L),error=function(e)
        list(status="FAILED",message=conditionMessage(e)))
    if(is.null(fit$score)) {
      rows[[length(rows)+1L]]<-data.frame(origin=origin,status="FIT_FAILED",
        candidate=c$id,rows=nrow(calendar)-origin,log_score=NA_real_,
        MAE=NA_real_,RMSE=NA_real_,observed_total=NA_real_,expected_total=NA_real_,
        observed_peak=NA_integer_,expected_peak=NA_integer_)
      output[[as.character(origin)]]<-list(tuning=tuning,choice=choice,fit=fit)
      next
    }
    forecast<-p21_score(ref,calendar,c,bases[[c$id]],fit$theta,cfg,mask,
      forecast_origin=origin)
    future<-(origin+1L):nrow(calendar)
    s<-forecast$states
    scores<-nb2_loglik(calendar$observed[future],s$expected_cases[future],
      ref$parameters$phi,rep(TRUE,length(future)))$per_week
    weekly_rows[[as.character(origin)]]<-data.frame(origin=origin,candidate=c$id,
      index=future,year=calendar$year[future],week=calendar$source_week[future],
      observed=calendar$observed[future],expected=s$expected_cases[future],
      incidence=s$incidence[future],S=s$S[future],I=s$I[future],R=s$R[future],
      reference_beta=s$reference_beta[future],q=s$q[future],
      effective_beta=s$effective_beta[future],log_score=scores)
    rows[[length(rows)+1L]]<-cbind(data.frame(origin=origin,status="FORECAST_SCORED",
      candidate=c$id),p21_score_rows(calendar$observed,
      forecast$states$expected_cases,ref$parameters$phi,future))
    output[[as.character(origin)]]<-list(tuning=tuning,choice=choice,fit=fit,
      forecast=forecast,training_mask=mask,forecast_index=future)
  }
  list(origins=output,summary=do.call(rbind,rows),
    weekly=if(length(weekly_rows)) do.call(rbind,weekly_rows) else data.frame(),
    rule="constant fitted log-q at each origin; no future outcome used for tuning or fit")
}
p21_bootstrap <- function(ref,calendar,candidate,basis,fit,cfg) {
  if(is.null(fit$score)) stop("Selected fit unavailable for bootstrap")
  set.seed(cfg$bootstrap$seed)
  rows<-list();means<-list();predictive<-list()
  for(i in seq_len(cfg$bootstrap$replicates)) {
    draw<-calendar
    draw$observed<-stats::rnbinom(nrow(calendar),mu=fit$score$states$expected_cases,
      size=ref$parameters$phi)
    trial<-tryCatch(p21_fit(ref,draw,candidate,basis,cfg,
      seed=cfg$bootstrap$seed+i),error=function(e)
        list(status="FAILED",message=conditionMessage(e)))
    rows[[i]]<-data.frame(replicate=i,status=trial$status,
      objective=if(is.null(trial$score)) NA_real_ else trial$score$objective,
      parameters=if(is.null(trial$theta)) "" else paste(signif(trial$theta,9),collapse=";"),
      message=if(is.null(trial$message)) "" else trial$message)
    if(!is.null(trial$score)) {
      means[[as.character(i)]]<-trial$score$states$expected_cases
      predictive[[as.character(i)]]<-stats::rnbinom(nrow(calendar),
        mu=trial$score$states$expected_cases,size=ref$parameters$phi)
    }
  }
  interval<-if(length(means)) {
    q<-cfg$bootstrap$interval
    m<-do.call(cbind,means);p<-do.call(cbind,predictive)
    data.frame(index=seq_len(nrow(calendar)),
      mean_low=apply(m,1,stats::quantile,probs=q[1]),
      mean_high=apply(m,1,stats::quantile,probs=q[2]),
      predictive_low=apply(p,1,stats::quantile,probs=q[1]),
      predictive_high=apply(p,1,stats::quantile,probs=q[2]))
  } else data.frame()
  list(replicates=do.call(rbind,rows),interval=interval,
    selection_uncertainty="excluded; candidate and smoothing fixed",
    historical_uncertainty="excluded; accepted endpoint frozen")
}
p21_figures <- function(ref,calendar,weekly,annual,forecast,out,historical_observed=NULL) {
  selected<-unique(weekly$candidate)
  grDevices::pdf(file.path(out,"phase21_frozen_and_future.pdf"),width=11,height=6)
  h<-ref$evaluation$states$expected_cases
  graphics::plot(seq_along(h),h,type="l",col="gray30",
    xlim=c(1,length(h)+nrow(calendar)),xlab="Transition row",
    ylim=range(c(h,historical_observed$observed,weekly$observed,weekly$expected),finite=TRUE),
    ylab="Expected reported cases",main=paste("Frozen historical and",selected,"future"))
  if(!is.null(historical_observed)&&nrow(historical_observed)==length(h))
    graphics::points(seq_along(h),historical_observed$observed,pch=16,cex=.2,col="gray65")
  graphics::lines(length(h)+seq_len(nrow(weekly)),weekly$expected,col="blue")
  graphics::points(length(h)+seq_len(nrow(weekly)),weekly$observed,pch=16,cex=.25,col="gray60")
  graphics::abline(v=length(h)+.5,lty=2);grDevices::dev.off()
  grDevices::pdf(file.path(out,"phase21_future_diagnostics.pdf"),width=11,height=10)
  op<-graphics::par(mfrow=c(4,1))
  graphics::plot(weekly$index,weekly$observed,type="l",col="gray40",
    xlab="Future transition",ylab="Reported cases")
  graphics::lines(weekly$index,weekly$expected,col="red")
  graphics::plot(weekly$index,weekly$reference_beta,type="l",
    xlab="Future transition",ylab="Accepted seasonal beta")
  graphics::plot(weekly$index,weekly$q,type="l",
    xlab="Future transition",ylab="Effective q")
  graphics::matplot(weekly$index,weekly[c("S","I","R")],type="l",lty=1,
    xlab="Future transition",ylab="State counts")
  graphics::par(op);grDevices::dev.off()
  grDevices::pdf(file.path(out,"phase21_annual_and_forecast.pdf"),width=10,height=7)
  op<-graphics::par(mfrow=c(2,1))
  graphics::plot(annual$year,annual$observed_total,type="b",
    xlab="Year",ylab="Annual reported cases")
  graphics::lines(annual$year,annual$expected_total,type="b",col="red")
  graphics::plot(annual$year,annual$expected_50-annual$observed_50,type="b",
    xlab="Year",ylab="Median cumulative timing error (weeks)")
  graphics::par(op);grDevices::dev.off()
  invisible(out)
}
