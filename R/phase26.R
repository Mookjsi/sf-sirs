# Phase 26 fits effective transmission corrections on a continuous trajectory.
p26_clock <- function(inputs) p24_clock(inputs)
p26_season <- function(year,week) ifelse(week>=27L,year,year-1L)
p26_centered_shape <- function(week) {
  x <- cbind(sin=sin(2*pi*(1:52)/52),cos=cos(2*pi*(1:52)/52))
  x <- sweep(x,2L,colMeans(x))
  x[week,,drop=FALSE]
}
p26_expand <- function(free,seasons) {
  if(length(free)!=length(seasons)-1L) stop("Wrong anchored correction length")
  setNames(c(free,-sum(free)),seasons)
}
p26_split <- function(theta,inputs,setup,cfg,model) {
  d <- length(setup$problem$map$order); q <- inputs$basis$dimension
  k <- length(cfg$historical_seasons); needed <- d+q+k-1L+if(model=="M2") 2L*k else 0L
  if(length(theta)!=needed) stop("Wrong candidate dimension")
  c <- p26_expand(theta[d+q+seq_len(k-1L)],cfg$historical_seasons)
  shape <- matrix(0,k,2L,dimnames=list(as.character(cfg$historical_seasons),c("sin","cos")))
  if(model=="M2") shape[] <- theta[(d+q+k):needed]
  list(coords=theta[seq_len(d)],post=theta[d+seq_len(q)],level=c,shape=shape)
}
p26_joint_theta <- function(inputs,source_cfg,cfg,model) {
  ref <- p22_references(inputs,source_cfg)$JOINT
  c(ref$coords,ref$post,rep(0,length(cfg$historical_seasons)-1L),
    if(model=="M2") rep(0,2L*length(cfg$historical_seasons)))
}
p26_score <- function(theta,inputs,setup,cfg,model,penalty,mask=NULL,
                      projected_seasons=integer(),blend=FALSE) {
  a <- p26_split(theta,inputs,setup,cfg,model)
  if(any(abs(a$level)>cfg$correction_bound+1e-10) ||
     any(abs(a$shape)>cfg$shape_bound+1e-10)) stop("Correction outside bounds")
  if(!p21_feasible(inputs$candidate,inputs$basis,a$post,
     list(q_log_bounds=inputs$q_log_bounds))) stop("Future correction outside bounds")
  problem <- setup$problem;p <- decode_parameters(a$coords,problem$map)
  profile <- if(isTRUE(problem$map$fourier))
    fourier_profile(p$free_coordinates,p$beta0,problem$basis) else
    spline_profile(p$free_coordinates,p$beta0,problem$basis)
  clock <- p26_clock(inputs); n <- length(problem$data$cases); N <- nrow(clock)
  if(is.null(mask)) mask <- c(FALSE,rep(TRUE,N-1L))
  if(length(mask)!=N || mask[1L] || anyNA(mask)) stop("Invalid likelihood mask")
  sy <- p26_season(clock$year,clock$week)
  logc <- numeric(N); shp <- numeric(N)
  historical <- seq_len(n)
  active <- sy[historical] %in% cfg$historical_seasons &
    !sy[historical] %in% projected_seasons
  logc[historical[active]] <- a$level[as.character(sy[historical][active])]
  if(model=="M2") {
    basis <- p26_centered_shape(clock$seasonal_week[historical])
    shp[historical[active]] <- rowSums(basis[active,,drop=FALSE]*
      a$shape[as.character(sy[historical][active]),,drop=FALSE])
  }
  raw <- logc+shp
  if(blend) {
    for(i in which(clock$week==27L & seq_len(N)<=n)) {
      if(i>1L) for(j in seq_len(cfg$boundary_blend_weeks)) {
        at <- i+j-1L
        if(at<=n) raw[at] <- (j/(cfg$boundary_blend_weeks+1L))*raw[at]+
          (1-j/(cfg$boundary_blend_weeks+1L))*raw[i-1L]
      }
    }
  }
  g <- as.vector(inputs$basis$X%*%a$post)
  future <- c(rep(0,n),g)
  transition <- exp(raw+future)
  fixed <- list(pop=problem$data$fixed_parameters$pop,
    mu=ref_or(p$mu,problem$data$fixed_parameters$mu),
    gamma=ref_or(p$gamma,problem$data$fixed_parameters$gamma),
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  states <- simulate_weekly_sirs(profile$beta_by_week,
    make_initial_state(p$S0,p$I0,fixed$pop),fixed,clock$seasonal_week,
    transition)
  obs <- c(problem$data$cases,inputs$calendar$observed)
  ll <- nb2_loglik(obs,states$expected_cases,p$phi,mask)
  seasonal_penalty <- if(isTRUE(problem$map$fourier)) 0 else
    cyclic_penalty(profile$coefficients,problem$model$lambda)
  future_penalty <- p21_penalty(inputs$candidate,inputs$basis,a$post)
  level_penalty <- penalty$level*sum(a$level^2)+penalty$change*sum(diff(a$level)^2)
  shape_penalty <- if(model=="M2") penalty$level*sum(a$shape^2)+
    penalty$change*sum(diff(a$shape)^2) else 0
  objective <- -ll$total+seasonal_penalty+future_penalty+level_penalty+shape_penalty
  list(objective=objective,loglik=ll$total,seasonal_penalty=seasonal_penalty,
    future_penalty=future_penalty,level_penalty=level_penalty,shape_penalty=shape_penalty,
    historical_loglik=sum(ll$included_per_week[seq_len(n)]),
    future_loglik=sum(ll$included_per_week[n+seq_len(N-n)]),
    states=states,observed=obs,per_week=ll$per_week,mask=mask,
    parameters=p,level=a$level,shape=a$shape,
    weekly=data.frame(clock,season=sy,observed=obs,expected=states$expected_cases,
      S=states$S,I=states$I,R=states$R,incidence=states$incidence,
      seasonal_beta=profile$beta_by_week[clock$seasonal_week],
      log_level=logc,log_shape=shp,log_correction=raw,
      future_q=exp(future),effective_beta=profile$beta_by_week[clock$seasonal_week]*transition,
      nb2_loglik=ll$per_week,in_likelihood=mask))
}
p26_constraints <- function(inputs,setup,cfg,model) {
  d <- length(setup$problem$map$order);q <- inputs$basis$dimension
  k <- length(cfg$historical_seasons);n <- d+q+k-1L+if(model=="M2") 2L*k else 0L
  rows <- list();rhs <- numeric()
  add <- function(v,b) {rows[[length(rows)+1L]] <<- v;rhs <<- c(rhs,b)}
  for(i in seq_len(d)) {
    v <- numeric(n);v[i] <- 1;add(v,setup$problem$map$lower[i]);add(-v,-setup$problem$map$upper[i])
  }
  control <- inputs$basis$control
  for(i in seq_len(nrow(control))) {
    v <- numeric(n);v[d+seq_len(q)] <- control[i,];add(v,inputs$q_log_bounds[1]);add(-v,-inputs$q_log_bounds[2])
  }
  for(i in seq_len(k-1L)) {
    v <- numeric(n);v[d+q+i] <- 1;add(v,-cfg$correction_bound);add(-v,-cfg$correction_bound)
  }
  v <- numeric(n);v[d+q+seq_len(k-1L)] <- 1
  add(v,-cfg$correction_bound);add(-v,-cfg$correction_bound)
  if(model=="M2") for(i in (d+q+k):n) {
    v <- numeric(n);v[i] <- 1;add(v,-cfg$shape_bound);add(-v,-cfg$shape_bound)
  }
  list(ui=do.call(rbind,rows),ci=rhs)
}
p26_season_metrics <- function(score,inputs,cfg,model) {
  clock <- p26_clock(inputs);z <- score$weekly
  ys <- 2012:2019
  out <- lapply(ys,function(y) {
    ix <- tryCatch(p24_window(clock,y),error=function(e) integer())
    if(!length(ix)) return(NULL)
    a <- z$observed[ix];b <- z$expected[ix]
    ao <- p24_stats(z$observed,clock,ix,model,y,"OBSERVED")
    bo <- p24_stats(z$expected,clock,ix,model,y,"FITTED")
    data.frame(model=model,season=y,observed_peak=ao$peak_season_week,
      fitted_peak=bo$peak_season_week,peak_error=bo$peak_season_week-ao$peak_season_week,
      observed_ties=ao$tied_season_weeks,fitted_ties=bo$tied_season_weeks,
      observed_height=ao$peak_height,fitted_height=bo$peak_height,
      height_error=bo$peak_height-ao$peak_height,
      observed_cog=ao$center_of_gravity,fitted_cog=bo$center_of_gravity,
      cog_error=bo$center_of_gravity-ao$center_of_gravity,
      observed_total=ao$total,fitted_total=bo$total,total_error=bo$total-ao$total,
      mae=mean(abs(a-b)),rmse=sqrt(mean((a-b)^2)),
      nb2_loglik=sum(z$nb2_loglik[ix]))
  })
  do.call(rbind,Filter(Negate(is.null),out))
}
p26_transition <- function(metrics,cfg) {
  a <- metrics[match(c(2015,2017),metrics$season),]
  op <- diff(a$observed_peak);fp <- diff(a$fitted_peak)
  oc <- diff(a$observed_cog);fc <- diff(a$fitted_cog)
  data.frame(model=a$model[1L],observed_peak_change=op,fitted_peak_change=fp,
    peak_change_error=fp-op,observed_cog_change=oc,fitted_cog_change=fc,
    cog_change_error=fc-oc,
    endpoint_peak_pass=all(abs(a$peak_error)<=cfg$endpoint_peak_tolerance),
    peak_change_pass=abs(fp-op)<=cfg$peak_change_tolerance,
    cog_change_pass=abs(fc-oc)<=cfg$cog_change_tolerance)
}
p26_periods <- function(score,model) {
  z <- score$weekly;period <- ifelse(z$year<=2014,"2013-2014",
    ifelse(z$year<=2017,"2015-2017",ifelse(z$year<=2019,"2018-2019","2020-2023")))
  do.call(rbind,lapply(unique(period),function(p) {
    ix <- which(period==p & z$in_likelihood)
    data.frame(model=model,period=p,rows=length(ix),
      nb2_loglik=sum(z$nb2_loglik[ix]),mae=mean(abs(z$observed[ix]-z$expected[ix])),
      rmse=sqrt(mean((z$observed[ix]-z$expected[ix])^2)))
  }))
}
p26_year_metrics <- function(score,model) {
  z <- score$weekly
  do.call(rbind,lapply(2013:2023,function(y) {
    ix <- which(z$year==y & z$in_likelihood)
    data.frame(model=model,year=y,rows=length(ix),
      mae=mean(abs(z$observed[ix]-z$expected[ix])),
      rmse=sqrt(mean((z$observed[ix]-z$expected[ix])^2)),
      nb2_loglik=sum(z$nb2_loglik[ix]),
      burden_error=sum(z$expected[ix]-z$observed[ix]))
  }))
}
p26_boundary <- function(score,cfg) {
  z <- score$weekly
  ix <- which(z$week==27L & z$row>1L & z$year<=2019)
  data.frame(season=z$season[ix],previous_log_correction=z$log_correction[ix-1L],
    current_log_correction=z$log_correction[ix],
    log_jump=z$log_correction[ix]-z$log_correction[ix-1L],
    flagged=abs(z$log_correction[ix]-z$log_correction[ix-1L])>cfg$boundary_jump_threshold)
}
p26_replay <- function(inputs,source_cfg,cfg) {
  refs <- p22_references(inputs,source_cfg)
  old <- p22_saved_joint(inputs,source_cfg)$fit$primary
  saved <- old$attempts[[old$best]]$score
  now <- p22_score(refs$JOINT$coords,refs$JOINT$post,inputs,refs$JOINT$setup)
  candidate_zero <- p26_score(p26_joint_theta(inputs,source_cfg,cfg,"M1"),inputs,
    refs$JOINT$setup,cfg,"M1",cfg$penalty_grid[2L,])
  columns <- c("S","I","R","incidence","expected_cases")
  err <- max(abs(as.matrix(now$states[,columns])-as.matrix(saved$states[,columns])))
  llerr <- abs(now$loglik-saved$loglik)
  objerr <- abs(now$objective-saved$objective)
  weekly <- utils::read.csv(cfg$phase24_weekly)
  w <- weekly[weekly$reference=="JOINT",]
  clock <- p26_clock(inputs)
  if(nrow(w)!=nrow(clock) || any(w$year!=clock$year | w$week!=clock$week))
    stop("Saved Phase 22 weekly row alignment changed")
  werr <- max(abs(w$expected-now$states$expected_cases))
  candidate_error <- max(abs(candidate_zero$states$expected_cases-now$states$expected_cases))
  candidate_objective_error <- abs(candidate_zero$objective-now$objective)
  peaks <- vapply(2015:2017,function(y) {
    ix <- p24_window(clock,y);which.max(now$states$expected_cases[ix])
  },integer(1))
  observed_peaks <- vapply(2015:2017,function(y) {
    ix <- p24_window(clock,y);which.max(now$observed[ix])
  },integer(1))
  p24 <- read_phase2_result(cfg$phase24)
  p25 <- read_phase2_result(cfg$phase25)
  if(!identical(p24$metadata$phase,24L) || !identical(p25$metadata$phase,25L))
    stop("Phase 24 or 25 artifact identity mismatch")
  checks <- data.frame(check=c("state_and_expected","likelihood","objective","weekly_csv",
    "candidate_zero_expected","candidate_zero_objective",
    "peak_2015","peak_2016","peak_2017",
    "observed_peak_2015","observed_peak_2016","observed_peak_2017"),
    error=c(err,llerr,objerr,werr,candidate_error,candidate_objective_error,
      abs(peaks-c(17L,14L,12L)),abs(observed_peaks-c(24L,14L,11L))),
    tolerance=c(rep(cfg$replay_tolerance,6),rep(0,6)))
  checks$passed <- checks$error<=checks$tolerance
  list(checks=checks,passed=all(checks$passed),score=now,source_hashes=tools::md5sum(c(
    unlist(source_cfg[c("reference","replay","plan","fit","joint_fit")]),
    cfg$phase24,cfg$phase25,cfg$phase24_weekly)))
}
p26_fit <- function(inputs,source_cfg,cfg,model,penalty,out,mask=NULL,
                    projected_seasons=integer(),holdout=FALSE,
                    setup_override=NULL,baseline_override=NULL,
                    extra_starts=list(),blend=FALSE) {
  log_file <- file.path(out,"fit_progress.log")
  log_progress <- function(message) {
    line <- sprintf("[%s] %s",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),message)
    cat(line,"\n",file=log_file,append=TRUE,sep="")
    cat(line,"\n",sep="")
    flush(stdout())
  }
  evaluations <- 0L
  last_progress <- Sys.time()
  current_start <- NA_integer_
  best_objective <- Inf
  log_progress(sprintf("Starting %s fit (penalty: %s)",model,penalty$id))
  setup <- if(is.null(setup_override)) p22_problem(inputs) else setup_override
  baseline <- if(is.null(baseline_override))
    p26_joint_theta(inputs,source_cfg,cfg,model) else baseline_override
  constraint <- p26_constraints(inputs,setup,cfg,model)
  feasible <- function(x) all(as.vector(constraint$ui%*%x-constraint$ci)>0)
  fn <- function(x) {
    z <- tryCatch(p26_score(x,inputs,setup,cfg,model,penalty,mask,projected_seasons,blend),error=function(e) NULL)
    evaluations <<- evaluations+1L
    if(!is.null(z) && is.finite(z$objective)) best_objective <<- min(best_objective,z$objective)
    now <- Sys.time()
    if(as.numeric(difftime(now,last_progress,units="secs"))>=30) {
      log_progress(sprintf("Start %s: %d objective evaluations; best objective %s",
        current_start,evaluations,if(is.finite(best_objective)) format(best_objective,digits=10) else "none"))
      last_progress <<- now
    }
    if(is.null(z) || !is.finite(z$objective)) 1e20 else z$objective
  }
  d <- length(setup$problem$map$order);q <- inputs$basis$dimension
  independent <- baseline
  design <- setup$cfg$initial
  design$omega <- inputs$ref$parameters$omega
  design$delta <- inputs$ref$parameters$delta
  design$free_coordinates <- rep(0,setup$problem$basis$K-1L)
  independent[seq_len(d)] <- encode_parameters(design,setup$problem$map)
  independent[seq_len(d)] <- pmax(setup$problem$map$lower+1e-5,
    pmin(setup$problem$map$upper-1e-5,independent[seq_len(d)]))
  independent[d+seq_len(q)] <- 0
  starts <- if(holdout) list(independent) else list(baseline,independent)
  starts <- c(starts,extra_starts)
  set.seed(cfg$seed+if(model=="M2") 1L else 0L)
  if(cfg$starts>length(starts)) for(i in (length(starts)+1L):cfg$starts) {
    x <- if(holdout) independent else baseline
    x[d+q+seq_len(length(cfg$historical_seasons)-1L)] <-
      rnorm(length(cfg$historical_seasons)-1L,sd=.025)
    if(model=="M2") x[(d+q+length(cfg$historical_seasons)):length(x)] <-
      rnorm(2L*length(cfg$historical_seasons),sd=.015)
    starts[[i]] <- x
  }
  attempts <- vector("list",length(starts))
  rows <- list()
  for(i in seq_along(starts)) {
    current_start <- i
    evaluations <- 0L
    best_objective <- Inf
    last_progress <- Sys.time()
    log_progress(sprintf("Start %d of %d: optimizing",i,length(starts)))
    x <- starts[[i]]
    result <- tryCatch({
      if(!feasible(x)) stop("Infeasible start")
      opt <- stats::constrOptim(x,fn,
        grad=function(y) p21_numeric_gradient(y,fn,feasible,cfg$gradient_step,
          adaptive=isTRUE(cfg$adaptive_gradient)),
        ui=constraint$ui,ci=constraint$ci,
        control=list(maxit=cfg$maxit,reltol=1e-8),outer.iterations=20L)
      score <- p26_score(opt$par,inputs,setup,cfg,model,penalty,mask,projected_seasons,blend)
      list(status=if(opt$convergence==0L) "CONVERGED" else "NONCONVERGED",
        initial=x,coordinates=opt$par,optimizer=opt,
        objective=score$objective,loglik=score$loglik,
        penalty=score$objective+score$loglik,
        replay_difference=score$objective-opt$value,
        movement=sqrt(sum((opt$par-x)^2)),
        minimum_slack=min(as.vector(constraint$ui%*%opt$par-constraint$ci)))
    },error=function(e) list(status="ERROR",initial=x,message=conditionMessage(e)))
    attempts[[i]] <- result
    rows[[i]] <- data.frame(start=i,status=result$status,
      convergence=if(!is.null(result$optimizer)) result$optimizer$convergence else NA_integer_,
      objective=if(!is.null(result$objective)) result$objective else NA_real_,
      loglik=if(!is.null(result$loglik)) result$loglik else NA_real_,
      penalty=if(!is.null(result$penalty)) result$penalty else NA_real_,
      movement=if(!is.null(result$movement)) result$movement else NA_real_,
      minimum_slack=if(!is.null(result$minimum_slack)) result$minimum_slack else NA_real_,
      replay_difference=if(!is.null(result$replay_difference)) result$replay_difference else NA_real_,
      message=if(!is.null(result$message)) result$message else "")
    utils::write.csv(do.call(rbind,rows),file.path(out,"fit_starts.csv"),row.names=FALSE)
    saveRDS(attempts,file.path(out,"attempts.partial.rds"))
    log_progress(sprintf("Start %d of %d: %s; objective %s; %d evaluations",
      i,length(starts),result$status,
      if(!is.null(result$objective)) format(result$objective,digits=10) else "unavailable",
      evaluations))
  }
  good <- which(vapply(attempts,function(a) identical(a$status,"CONVERGED") &&
    is.finite(a$objective),logical(1)))
  if(!length(good)) stop("No converged candidate; diagnostics retained in ",out)
  best <- good[which.min(vapply(attempts[good],function(a) a$objective,numeric(1)))]
  endpoint <- attempts[[best]]$coordinates
  score <- p26_score(endpoint,inputs,setup,cfg,model,penalty,mask,projected_seasons,blend)
  diagnose <- function(endpoint,score) {
    set.seed(cfg$seed+19L)
    best_delta <- 0
    best_candidate <- NULL
    rows <- do.call(rbind,lapply(cfg$descent_steps,function(step) {
    dirs <- rbind(diag(length(endpoint)),
      -diag(length(endpoint)),matrix(rnorm(10L*length(endpoint)),10L))
    dirs <- dirs/sqrt(rowSums(dirs^2))
    vals <- apply(dirs,1L,function(v) {
      y <- endpoint+step*v
      if(!feasible(y)) return(NA_real_)
      delta <- fn(y)-score$objective
      if(is.finite(delta) && delta<best_delta) {
        best_delta <<- delta
        best_candidate <<- y
      }
      delta
    })
    data.frame(step=step,minimum_delta=if(all(is.na(vals))) NA_real_ else min(vals,na.rm=TRUE),
      feasible_directions=sum(!is.na(vals)))
    }))
    attr(rows,"best_candidate") <- best_candidate
    attr(rows,"best_delta") <- best_delta
    rows
  }
  descent <- diagnose(endpoint,score)
  refined <- FALSE
  if(isTRUE(cfg$refine_descent) &&
     any(descent$minimum_delta < -cfg$descent_tolerance,na.rm=TRUE)) {
    log_progress("Local descent detected; refining the selected endpoint")
    refinement <- tryCatch(stats::constrOptim(endpoint,fn,
      grad=function(y) p21_numeric_gradient(y,fn,feasible,
        cfg$refine_gradient_step,adaptive=TRUE),
      ui=constraint$ui,ci=constraint$ci,
      control=list(maxit=cfg$refine_maxit,reltol=1e-10),
      outer.iterations=20L),error=function(e) e)
    if(inherits(refinement,"error")) {
      log_progress(paste("Refinement failed:",conditionMessage(refinement)))
    } else {
      refined_score <- p26_score(refinement$par,inputs,setup,cfg,model,
        penalty,mask,projected_seasons,blend)
      if(is.finite(refined_score$objective) &&
         refined_score$objective < score$objective) {
        endpoint <- refinement$par
        score <- refined_score
        refined <- TRUE
        descent <- diagnose(endpoint,score)
        log_progress(sprintf("Refined objective %.9f; minimum local delta %.9g",
          score$objective,min(descent$minimum_delta,na.rm=TRUE)))
      }
    }
  }
  if(!is.null(cfg$descent_repair_max_steps) &&
     cfg$descent_repair_max_steps>0L &&
     any(descent$minimum_delta < -cfg$descent_tolerance,na.rm=TRUE)) {
    repair <- data.frame(iteration=integer(),objective=numeric(),
      minimum_delta=numeric())
    for(iteration in seq_len(cfg$descent_repair_max_steps)) {
      candidate <- attr(descent,"best_candidate")
      if(is.null(candidate) ||
         attr(descent,"best_delta") >= -cfg$descent_tolerance) break
      candidate_score <- p26_score(candidate,inputs,setup,cfg,model,
        penalty,mask,projected_seasons,blend)
      if(!is.finite(candidate_score$objective) ||
         candidate_score$objective >= score$objective) break
      endpoint <- candidate
      score <- candidate_score
      refined <- TRUE
      descent <- diagnose(endpoint,score)
      repair <- rbind(repair,data.frame(iteration=iteration,
        objective=score$objective,
        minimum_delta=min(descent$minimum_delta,na.rm=TRUE)))
      log_progress(sprintf("Descent repair %d: objective %.9f; minimum local delta %.9g",
        iteration,score$objective,tail(repair$minimum_delta,1L)))
    }
    utils::write.csv(repair,file.path(out,"descent_repair.csv"),row.names=FALSE)
  }
  utils::write.csv(descent,file.path(out,"numerical_diagnostics.csv"),row.names=FALSE)
  if(any(descent$minimum_delta < -cfg$descent_tolerance,na.rm=TRUE))
    warning("Feasible local descent detected; inspect numerical_diagnostics.csv")
  names(endpoint) <- c(setup$problem$map$order,paste0("future_q_",seq_len(q)),
    paste0("level_",cfg$historical_seasons[-length(cfg$historical_seasons)]),
    if(model=="M2") paste0(rep(c("sin_","cos_"),each=length(cfg$historical_seasons)),
      rep(cfg$historical_seasons,2L)))
  params <- data.frame(parameter=names(endpoint),value=endpoint,
    baseline=baseline,movement=endpoint-baseline)
  params$lower <- NA_real_;params$upper <- NA_real_
  params$lower[seq_len(d)] <- setup$problem$map$lower
  params$upper[seq_len(d)] <- setup$problem$map$upper
  level_ix <- d+q+seq_len(length(cfg$historical_seasons)-1L)
  params$lower[level_ix] <- -cfg$correction_bound
  params$upper[level_ix] <- cfg$correction_bound
  if(model=="M2") {
    shape_ix <- (max(level_ix)+1L):nrow(params)
    params$lower[shape_ix] <- -cfg$shape_bound
    params$upper[shape_ix] <- cfg$shape_bound
  }
  anchored <- p26_expand(endpoint[level_ix],cfg$historical_seasons)
  params <- rbind(params,data.frame(parameter=paste0("level_",tail(cfg$historical_seasons,1L),"_anchored"),
    value=unname(tail(anchored,1L)),baseline=0,movement=unname(tail(anchored,1L)),
    lower=-cfg$correction_bound,upper=cfg$correction_bound))
  params$distance_to_bound <- pmin(params$value-params$lower,params$upper-params$value)
  utils::write.csv(params,file.path(out,"parameter_table.csv"),row.names=FALSE)
  list(status="FIT_SAVED",model=model,
    family=if(isTRUE(setup$problem$map$fourier))
      paste0("FOURIER_",setup$problem$basis$harmonics) else "SPL",
    penalty=penalty,best=best,
    attempts=attempts,endpoint=unname(endpoint),score=score,refined=refined,
    numerical=descent,parameter_table=params,holdout=holdout)
}
p26_evaluate <- function(fit,inputs,source_cfg,cfg,out) {
  setup <- p22_problem(inputs);model <- fit$model
  score <- p26_score(fit$endpoint,inputs,setup,cfg,model,fit$penalty)
  metrics <- p26_season_metrics(score,inputs,cfg,model)
  transition <- p26_transition(metrics,cfg)
  boundary <- p26_boundary(score,cfg)
  base <- p22_references(inputs,source_cfg)$JOINT
  m0theta <- p26_joint_theta(inputs,source_cfg,cfg,"M1")
  m0score <- p26_score(m0theta,inputs,setup,cfg,"M1",fit$penalty)
  m0metrics <- p26_season_metrics(m0score,inputs,cfg,"M0")
  targeted <- mean(metrics$mae[metrics$season %in% 2015:2017]) <
    mean(m0metrics$mae[m0metrics$season %in% 2015:2017])
  years <- p26_year_metrics(score,model)
  baseline_years <- p26_year_metrics(m0score,"M0")
  other <- years$mae[years$year %in% c(2013,2014,2018,2019)] /
    baseline_years$mae[baseline_years$year %in% c(2013,2014,2018,2019)]
  no_material <- all(other<=1+cfg$other_year_mae_deterioration)
  transition$targeted_mae_improved <- targeted
  transition$other_historical_years_preserved <- no_material
  transition$timing_recovered <- transition$endpoint_peak_pass & transition$peak_change_pass
  transition$waveform_adequate <- transition$timing_recovered & transition$cog_change_pass &
    targeted & no_material
  utils::write.csv(score$weekly,file.path(out,"weekly_fit.csv"),row.names=FALSE)
  utils::write.csv(metrics,file.path(out,"season_metrics.csv"),row.names=FALSE)
  utils::write.csv(rbind(baseline_years,years),file.path(out,"year_metrics.csv"),row.names=FALSE)
  utils::write.csv(transition,file.path(out,"transition_metrics.csv"),row.names=FALSE)
  utils::write.csv(rbind(p26_periods(m0score,"M0"),p26_periods(score,model)),
    file.path(out,"period_fit_comparison.csv"),row.names=FALSE)
  if(any(boundary$flagged)) {
    blend <- p26_score(fit$endpoint,inputs,setup,cfg,model,fit$penalty,blend=TRUE)
    bmetrics <- p26_season_metrics(blend,inputs,cfg,model)
    bt <- p26_transition(bmetrics,cfg)
    sensitivity <- data.frame(boundary,blended_objective=blend$objective,
      blended_peak_change=bt$fitted_peak_change,blended_cog_change=bt$fitted_cog_change)
    utils::write.csv(sensitivity,file.path(out,"boundary_sensitivity.csv"),row.names=FALSE)
  }
  p26_figures(score,metrics,out)
  list(status="EVALUATED",model=model,score=score,metrics=metrics,
    transition=transition,boundary=boundary,baseline=m0score,
    baseline_metrics=m0metrics,periods=p26_periods(score,model),years=years,
    baseline_years=baseline_years)
}
p26_figures <- function(score,metrics,out) {
  z <- score$weekly
  grDevices::pdf(file.path(out,"weekly_2015_2017.pdf"),width=11,height=8)
  op <- par(mfrow=c(3,1),mar=c(3,4,2,1))
  for(y in 2015:2017) {
    a <- z[z$season==y,]
    plot(seq_len(nrow(a)),a$observed,type="l",col="black",xlab="Season week",
      ylab="Cases",main=paste0(y,"/",substr(y+1L,3,4)))
    lines(seq_len(nrow(a)),a$expected,col="blue")
    lines(seq_len(nrow(a)),a$observed-a$expected,col="red",lty=3)
    legend("topright",c("Observed","Expected","Residual"),col=c("black","blue","red"),lty=c(1,1,3),bty="n")
  }
  par(op);grDevices::dev.off()
  grDevices::pdf(file.path(out,"timing_and_cog.pdf"))
  op <- par(mfrow=c(2,1))
  matplot(metrics$season,cbind(metrics$observed_peak,metrics$fitted_peak),type="b",
    pch=1:2,xlab="Season",ylab="Peak week")
  matplot(metrics$season,cbind(metrics$observed_cog,metrics$fitted_cog),type="b",
    pch=1:2,xlab="Season",ylab="Center of gravity")
  par(op);grDevices::dev.off()
  grDevices::pdf(file.path(out,"historical_correction.pdf"))
  a <- z[z$year<=2019,]
  plot(a$row,a$log_correction,type="l",xlab="Weekly row",ylab="Log transmission correction")
  abline(h=0,lty=3);grDevices::dev.off()
}
p26_variant <- function(score,inputs,u=0,v=0,beta=NULL) {
  p <- score$parameters;z <- score$weekly
  pop <- inputs$replay$data$fixed_parameters$pop
  fixed <- list(pop=pop,mu=ref_or(p$mu,inputs$replay$data$fixed_parameters$mu),
    gamma=ref_or(p$gamma,inputs$replay$data$fixed_parameters$gamma),
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  if(is.null(beta)) beta <- vapply(1:52,function(w)
    z$seasonal_beta[which(z$seasonal_week==w)[1L]],numeric(1))
  initial <- make_initial_state(p$S0,p$I0,pop)
  initial["S"] <- initial["S"]*exp(u)
  initial["R"] <- pop-initial["S"]-initial["I"]
  if(any(initial<=0)) stop("Infeasible initial-state perturbation")
  transition <- z$effective_beta/z$seasonal_beta*exp(v)
  states <- simulate_weekly_sirs(beta,initial,fixed,z$seasonal_week,transition)
  list(states=states,expected=states$expected_cases,
    loglik=nb2_loglik(z$observed,states$expected_cases,p$phi,z$in_likelihood)$total)
}
p26_candidate_audits <- function(fit,inputs,cfg,out) {
  saveRDS(fit$endpoint,file.path(out,"candidate_endpoint.rds"))
  score <- p26_score(fit$endpoint,inputs,p22_problem(inputs),cfg,fit$model,fit$penalty)
  clock <- p26_clock(inputs);years <- 2015:2019
  # Fixed grid is declared in code before any candidate audit is run.
  grid <- seq(-log(1.2),log(1.2),length.out=17L)
  mechanism <- do.call(rbind,lapply(years,function(y) {
    ix <- p24_window(clock,y)
    do.call(rbind,lapply(grid,function(v) {
      a <- p26_variant(score,inputs,v=v)$expected[ix]
      local <- which(a>=c(-Inf,head(a,-1L)) & a>=c(tail(a,-1L),-Inf))
      data.frame(season=y,log_beta_multiplier=v,peak_week=which.max(a),
        peak_height=max(a),local_peak_weeks=paste(local,collapse=";"),
        local_peak_count=length(local))
    }))
  }))
  utils::write.csv(mechanism,file.path(out,"candidate_mechanism.csv"),row.names=FALSE)
  steps <- log(1.2)/2^(2:6)
  sensitivity <- do.call(rbind,lapply(years,function(y) {
    ix <- p24_window(clock,y)
    do.call(rbind,lapply(steps,function(h) {
      a <- p26_variant(score,inputs,u=h)$expected[ix]
      b <- p26_variant(score,inputs,u=-h)$expected[ix]
      c <- p26_variant(score,inputs,v=h)$expected[ix]
      d <- p26_variant(score,inputs,v=-h)$expected[ix]
      cog <- function(x) sum(seq_along(x)*x)/sum(x)
      data.frame(season=y,step=h,
        S_peak_slope=(which.max(a)-which.max(b))/(2*h),
        beta_peak_slope=(which.max(c)-which.max(d))/(2*h),
        S_cog_slope=(cog(a)-cog(b))/(2*h),
        beta_cog_slope=(cog(c)-cog(d))/(2*h))
    }))
  }))
  utils::write.csv(sensitivity,file.path(out,"candidate_sensitivity.csv"),row.names=FALSE)
  beta <- vapply(1:52,function(w) score$weekly$seasonal_beta[
    which(score$weekly$seasonal_week==w)[1L]],numeric(1))
  harmonics <- do.call(rbind,lapply(1:2,function(h) {
    w <- 1:52
    X <- cbind(1,cos(2*pi*w/52),sin(2*pi*w/52))
    if(h==2L) X <- cbind(X,cos(4*pi*w/52),sin(4*pi*w/52))
    approx <- exp(as.vector(X%*%qr.solve(X,log(beta))))
    candidate <- p26_variant(score,inputs,beta=approx)
    data.frame(harmonics=h,policy="fixed_candidate_harmonic_projection",
      beta_log_rmse=sqrt(mean((log(approx)-log(beta))^2)),
      nb2_loglik=candidate$loglik,
      delta_loglik_from_candidate=candidate$loglik-score$loglik)
  }))
  utils::write.csv(harmonics,file.path(out,"candidate_harmonics.csv"),row.names=FALSE)
  source_setup <- p22_problem(inputs)
  parts <- p26_split(fit$endpoint,inputs,source_setup,cfg,fit$model)
  refits <- lapply(c("A","AS"),function(family) {
    setup <- p22_problem(inputs,family)
    X <- do.call(cbind,lapply(seq_len(setup$problem$basis$harmonics),function(j)
      cbind(cos(2*pi*j*(0:51)/52),sin(2*pi*j*(0:51)/52))))
    parameters <- score$parameters
    parameters$beta0 <- mean(beta)
    parameters$free_coordinates <- as.vector(qr.solve(X,log(beta)))
    coords <- encode_parameters(parameters,setup$problem$map)
    baseline <- c(coords,parts$post,unname(head(parts$level,-1L)),
      if(fit$model=="M2") as.vector(parts$shape))
    subdir <- file.path(out,paste0("refit_",family))
    dir.create(subdir)
    run <- p26_fit(inputs,NULL,cfg,fit$model,fit$penalty,subdir,
      setup_override=setup,baseline_override=baseline)
    saveRDS(run,file.path(subdir,"fit.rds"))
    data.frame(family=family,loglik=run$score$loglik,
      objective=run$score$objective,historical_loglik=run$score$historical_loglik,
      future_loglik=run$score$future_loglik,
      seasonal_maxima=p22_fourier_maxima(run$score$parameters,setup$problem$basis))
  })
  refits <- do.call(rbind,refits)
  utils::write.csv(refits,file.path(out,"candidate_fourier_refits.csv"),row.names=FALSE)
  writeLines(c("# Exact-candidate research-question diagnostics","",
    "The mechanism grid recomputes continuous trajectories under the saved candidate.",
    "Sensitivity slopes are finite differences of the selected weekly peak and center of gravity.",
    "The annual/semiannual comparison jointly refits the same candidate correction and Phase 22 nuisance policy under A and AS Fourier forcing.",
    "The harmonic projection table is a separate fixed-parameter diagnostic.",
    "Peak jumps are discrete and are not automatically classified as dominance exchanges.",
    paste("Candidate endpoint MD5:",unname(tools::md5sum(file.path(out,"candidate_endpoint.rds"))))),
    file.path(out,"candidate_audits_report.md"))
  list(status="CANDIDATE_AUDITS_SAVED",model=fit$model,
    mechanism=mechanism,sensitivity=sensitivity,harmonics=harmonics,
    fourier_refits=refits)
}
