# Paper Table 1 is the Phase 14 interval audit; Table 2 is Phase 12A;
# Table 3 is the saved unpenalized Phase 5/6 Fourier comparison.
p22_historical_replay <- function(cfg) {
  paths <- unlist(cfg[c("paper","bundle","historical_mechanism","historical_sensitivity",
    "historical_fourier_fixed","historical_fourier_joint","annual_profiles",
    "annual_scenarios","annual_windows")])
  if(any(!file.exists(paths))) stop("Paper replay input missing: ",paste(paths[!file.exists(paths)],collapse=", "))
  bundle <- p22_bundle_audit(cfg)
  m <- utils::read.csv(cfg$historical_mechanism)
  classes <- c(DISPLACEMENT_SUPPORTED="within_branch_displacement",
    STABLE_SAME_BRANCH="stable_same_branch",
    WINDOW_SELECTION_SUPPORTED="window_selection",UNRESOLVED="unresolved")
  expected <- matrix(c(14,26,1,2,12,16,1,14),nrow=4,
    dimnames=list(unname(classes),c("2018","2019")))
  actual <- table(factor(classes[m$revised_classification],levels=rownames(expected)),
    factor(m$year,levels=2018:2019))
  if(anyNA(classes[m$revised_classification]) || !all(actual==expected))
    stop("Paper Table 1 differs from saved Phase 14 classification")
  table1 <- data.frame(classification=rep(rownames(expected),2),
    year=rep(2018:2019,each=4),count=as.vector(expected),
    denominator_all=43L)
  s <- utils::read.csv(cfg$historical_sensitivity)
  want <- rbind(data.frame(year=2018,h=cfg$sensitivity_steps[-1L],
    E_S=c(-6.9870,-7.0514,-7.0434,-7.0414),
    E_beta=c(-156.9306,-107.0670,-101.7325,-103.0692)),
    data.frame(year=2019,h=cfg$sensitivity_steps[-1L],
    E_S=c(-3.7701,-3.7291,-3.7189,-3.7163),
    E_beta=c(-107.2845,-75.2268,-69.2034,-67.9179)))
  want$step_key <- round(want$h,10)
  s$step_key <- round(s$h,10)
  got <- merge(want,s[s$method=="quadratic" & s$valid,c("year","step_key","E_S","E_beta")],
    by=c("year","step_key"),suffixes=c("_paper","_saved"))
  if(nrow(got)!=8L || any(abs(got$E_S_paper-got$E_S_saved)>5e-4) ||
     any(abs(got$E_beta_paper-got$E_beta_saved)>5e-4))
    stop("Paper Table 2 differs from saved Phase 12A quadratic estimates")
  largest <- s[s$method=="quadratic" & s$year %in% 2018:2019 &
    abs(s$h-cfg$sensitivity_steps[1L])<1e-10,]
  if(nrow(largest)!=2L || any(largest$valid))
    stop("Paper Table 2 largest-step exclusion differs from saved Phase 12A")
  fixed <- utils::read.csv(cfg$historical_fourier_fixed)
  joint <- utils::read.csv(cfg$historical_fourier_joint)
  paper <- rbind(data.frame(nuisance="fixed",model=c("A","AS"),
    loglik=c(-2785.392264,-2717.820745)),
    data.frame(nuisance="joint",model=c("A","AS"),
      loglik=c(-2761.695839,-2692.888379)))
  saved <- rbind(data.frame(nuisance="fixed",model=fixed$model,
    loglik=fixed$logLik,parameters=fixed$estimated_parameters,
    AIC=fixed$AIC,BIC=fixed$BIC),
    data.frame(nuisance="joint",model=joint$model,
      loglik=joint$logLik,parameters=joint$estimated_parameters,
      AIC=joint$AIC,BIC=joint$BIC))
  table3 <- merge(paper,saved,by=c("nuisance","model"),suffixes=c("_paper","_saved"))
  if(nrow(table3)!=4L || any(abs(table3$loglik_paper-table3$loglik_saved)>5e-4))
    stop("Paper Table 3 differs from saved Phase 5/6 Fourier fits")
  profiles <- utils::read.csv(cfg$annual_profiles)
  pp <- profiles[profiles$profile_id %in% c("profile-004","profile-009"),]
  if(nrow(pp)!=2L || any(abs(pp$b-c(60.9560072,91.4340108))>1e-6) ||
     !identical(pp$continuous_maxima,c(1L,1L)) ||
     any(abs(pp$amplitude-.3)>1e-12) ||
     !identical(as.integer(pp$phase_week),c(1L,14L)))
    stop("Annual forcing examples differ from saved Phase 16 profiles")
  # For A>0, d log(beta)/dw is proportional to -sin(2*pi*(w-phi)/52).
  # Its two stationary points per period have opposite curvature, so one is a maximum.
  pp$analytic_continuous_maxima <- 1L
  scenarios <- utils::read.csv(cfg$annual_scenarios)
  windows <- utils::read.csv(cfg$annual_windows)
  sw <- merge(scenarios,windows,by.x="scenario_id",by.y="scenario")
  sw$a <- sw$a.x;sw$status <- sw$status.y
  sw <- sw[sw$profile_id %in% pp$profile_id & sw$year==4L & sw$axis %in% c("baseline","beta"),]
  expected_annual <- data.frame(profile_id=c(rep("profile-004",4),rep("profile-009",7)),
    a=c(1/1.2,.8721959495,.9128709292,1,1/1.2,.9128709292,
      .9554427922,1,1.0466351394,1.0954451150,1.2),
    week=c(12,7,4,3,23,20,18,15,13,10,5))
  expected_annual$a_key <- round(expected_annual$a,6)
  sw$a_key <- round(sw$a,6)
  annual <- merge(expected_annual,sw[,c("profile_id","a_key","week","status")],
    by=c("profile_id","a_key"),suffixes=c("_paper","_saved"))
  if(nrow(annual)!=nrow(expected_annual) || any(annual$week_paper!=annual$week_saved))
    stop("Paper annual profile paths differ from saved Phase 16 windows")
  list(status="PAPER_TABLES_REPLAYED",table1=table1,table2=got,
    table3=table3,annual=annual,forcing=pp,input_hashes=tools::md5sum(paths),
    bundle_audit=bundle,
    tolerances=c(table1_exact=0,table2_absolute=5e-4,
      table3_absolute=5e-4,annual_week_exact=0))
}

p22_bundle_audit <- function(cfg) {
  zip <- cfg$bundle
  index <- utils::unzip(zip,list=TRUE)
  members <- c("fit__result.rds","plan__result.rds","evaluate__result.rds")
  names <- paste0("artifacts/phase21-bundle-003/",members)
  if(!all(names %in% index$Name)) stop("Phase 21 bundle lacks selected evidence")
  dest <- tempfile("phase22-bundle-");dir.create(dest)
  utils::unzip(zip,files=names,exdir=dest)
  extracted <- file.path(dest,names)
  local <- file.path("artifacts/phase21-bundle-003",members)
  upstream <- c(cfg$fit,cfg$plan,cfg$evaluation)
  hashes <- data.frame(member=members,
    zip_md5=unname(tools::md5sum(extracted)),
    bundle_dir_md5=unname(tools::md5sum(local)),
    upstream_md5=unname(tools::md5sum(upstream)))
  if(anyNA(hashes$zip_md5) || any(hashes$zip_md5!=hashes$bundle_dir_md5) ||
     any(hashes$zip_md5!=hashes$upstream_md5))
    stop("Phase 21 bundle evidence differs from saved fit/plan/evaluation")
  hashes
}

p22_saved_joint <- function(inputs,cfg) {
  saved <- read_phase2_result(cfg$joint_fit)
  fit <- saved$result
  if(!identical(fit$type,"phase22_fit") ||
     !identical(unname(fit$hashes),unname(inputs$hashes)) ||
     !identical(fit$primary$status,"ENDPOINTS_RECORDED") ||
     !length(fit$primary$competitive)) stop("Phase 22 primary fit or ancestry invalid")
  list(fit=fit,hash=unname(tools::md5sum(cfg$joint_fit)))
}

p22_references <- function(inputs,cfg) {
  jf <- p22_saved_joint(inputs,cfg)$fit
  setup <- p22_problem(inputs)
  hist <- inputs$ref$record$retained_coordinates[setup$problem$map$order]
  joint <- jf$primary$attempts[[jf$primary$best]]
  list(HIST=list(coords=hist,post=rep(0,inputs$basis$dimension),
      setup=setup,label="2013-2019 accepted lambda=1"),
    FROZEN=list(coords=hist,post=inputs$frozen$theta,
      setup=setup,label="Phase 21 selected q; adequacy preserved"),
    JOINT=list(coords=joint$coordinates[seq_along(setup$problem$map$order)],
      post=tail(joint$coordinates,inputs$basis$dimension),setup=setup,
      label="2013-2023 continuous joint fit"))
}

p22_fit_diagnostics <- function(inputs,cfg) {
  saved <- p22_saved_joint(inputs,cfg)
  rows <- steps <- list()
  for(policy in c("primary","biology_sensitivity")) {
    fit <- saved$fit[[policy]];setup <- fit$setup
    d <- length(setup$problem$map$order)
    for(i in seq_along(fit$attempts)) {
      a <- fit$attempts[[i]]
      if(is.null(a$score) || is.null(a$coordinates)) {
        rows[[length(rows)+1L]] <- data.frame(policy=policy,start=i,status=a$status,
          objective=NA_real_,replay_difference=NA_real_,optimizer_convergence=NA_integer_,
          feasible=FALSE,competitive=FALSE,message=ref_or(a$message,""))
        next
      }
      x <- a$coordinates
      replay <- p22_score(x[seq_len(d)],x[-seq_len(d)],inputs,setup)
      feasible <- all(x[seq_len(d)]>=setup$problem$map$lower) &&
        all(x[seq_len(d)]<=setup$problem$map$upper) &&
        p21_feasible(inputs$candidate,inputs$basis,x[-seq_len(d)],
          list(q_log_bounds=inputs$q_log_bounds))
      rows[[length(rows)+1L]] <- data.frame(policy=policy,start=i,status=a$status,
        objective=replay$objective,replay_difference=replay$objective-a$score$objective,
        optimizer_convergence=a$convergence,feasible=feasible,
        competitive=i %in% fit$competitive,message=ref_or(a$message,""))
      if(!i %in% fit$competitive) next
      fn <- function(z) p22_score(z[seq_len(d)],z[-seq_len(d)],inputs,setup)$objective
      directions <- list(seasonal_scale=c(1,rep(0,length(x)-1L)),
        post_intercept=c(rep(0,d),1,rep(0,length(x)-d-1L)))
      for(k in names(directions)) for(h in cfg$derivative_steps) {
        v <- directions[[k]];plus <- x+h*v;minus <- x-h*v
        value <- tryCatch((fn(plus)-fn(minus))/(2*h),error=function(e) NA_real_)
        steps[[length(steps)+1L]] <- data.frame(policy=policy,start=i,
          direction=k,step=h,derivative=value)
      }
    }
  }
  attempts <- do.call(rbind,rows);derivatives <- do.call(rbind,steps)
  derivatives$consistent <- NA
  for(k in unique(paste(derivatives$policy,derivatives$start,derivatives$direction))) {
    ix <- which(paste(derivatives$policy,derivatives$start,derivatives$direction)==k)
    v <- derivatives$derivative[ix];t <- cfg$derivative_tolerance
    derivatives$consistent[ix] <- all(is.finite(v)) &&
      max(v)-min(v)<=t$absolute+t$relative*max(abs(v))
  }
  endpoint <- do.call(rbind,lapply(c("primary","biology_sensitivity"),function(policy) {
    fit <- saved$fit[[policy]];best <- fit$attempts[[fit$best]]$score
    do.call(rbind,lapply(fit$competitive,function(i) {
      e <- fit$attempts[[i]]$score
      data.frame(policy=policy,start=i,objective_gap=e$objective-best$objective,
        max_seasonal_log_beta_difference=max(abs(log(e$beta/best$beta))),
        max_post_log_q_difference=max(abs(log(e$q/best$q))),
        max_historical_S_difference=max(abs(e$states$S[1:364]-best$states$S[1:364])))
    }))
  }))
  list(status="DIAGNOSTICS_RECORDED",attempts=attempts,
    derivatives=derivatives,competitive_endpoints=endpoint,
    fit_hash=saved$hash)
}

p22_period_comparison <- function(inputs,cfg) {
  refs <- p22_references(inputs,cfg)
  rows <- weekly <- seasonal <- states <- correction <- alternative <- list()
  old_env <- new.env(parent=baseenv())
  sys.source("config/accepted_extension.R",old_env)
  alternate_calendar <- ae_calendar(inputs$saved$panel$result$island_cases,
    old_env$accepted_extension_config(),"map_to_week1")
  alternate_inputs <- inputs;alternate_inputs$calendar <- alternate_calendar
  hist_time <- inputs$replay$data$time_markers
  years <- c(hist_time$year,inputs$calendar$year)
  weeks <- c(hist_time$week,inputs$calendar$source_week)
  for(name in names(refs)) {
    r <- refs[[name]]
    s <- p22_score(r$coords,r$post,inputs,r$setup)
    rows[[name]] <- data.frame(reference=name,
      historical_loglik=s$historical_loglik,extension_loglik=s$extension_loglik,
      raw_loglik=s$loglik,seasonal_penalty=s$seasonal_penalty,
      post_penalty=s$post_penalty,objective=s$objective,
      frozen_adequacy=if(name=="FROZEN")
        inputs$saved$evaluation$result$adequacy$status else NA_character_)
    weekly[[name]] <- data.frame(reference=name,row=seq_along(years),year=years,
      week=weeks,observed=s$observed,expected=s$states$expected_cases,
      incidence=s$states$incidence,per_week_loglik=s$per_week,
      likelihood_included=s$mask)
    seasonal[[name]] <- data.frame(reference=name,week=1:52,beta=s$beta)
    states[[name]] <- data.frame(reference=name,year=2019,week=52,
      S=s$states$S[364],I=s$states$I[364],R=s$states$R[364])
    correction[[name]] <- data.frame(reference=name,year=inputs$calendar$year,
      week=inputs$calendar$source_week,q=s$q)
    alt <- p22_score(r$coords,r$post,alternate_inputs,r$setup)
    alternative[[name]] <- data.frame(reference=name,policy="map_to_week1",
      historical_loglik=alt$historical_loglik,
      extension_loglik=alt$extension_loglik,
      difference_from_primary=alt$loglik-s$loglik)
  }
  list(status="SCORED",scores=do.call(rbind,rows),weekly=do.call(rbind,weekly),
    seasonal=do.call(rbind,seasonal),terminal_state=do.call(rbind,states),
    correction=do.call(rbind,correction),
    calendar_sensitivity=do.call(rbind,alternative))
}

p22_fourier_fit <- function(inputs,cfg) {
  records <- list();summary <- annual <- list();counter <- 0L
  total <- 2L*length(cfg$q_smoothness)*2L
  for(nuisance in c("fixed","joint")) for(lambda in cfg$q_smoothness)
    for(family in c("A","AS")) {
      key <- paste(nuisance,lambda,family,sep="/")
      counter <- counter+1L
      p22_progress("Fourier fit %d/%d: %s",counter,total,key)
      inp <- inputs;inp$candidate$lambda <- lambda
      setup <- p22_problem(inp,family,fixed_nuisance=nuisance=="fixed")
      fit <- p22_fit(inp,setup,cfg,seed=cfg$seed+counter)
      records[[key]] <- fit
      p22_progress("Fourier fit %d/%d finished: %s",counter,total,fit$status)
      if(identical(fit$status,"ENDPOINTS_RECORDED")) {
        s <- fit$attempts[[fit$best]]$score
        summary[[key]] <- data.frame(nuisance=nuisance,lambda=lambda,family=family,
          status=fit$status,historical_loglik=s$historical_loglik,
          extension_loglik=s$extension_loglik,raw_loglik=s$loglik,
          post_penalty=s$post_penalty,objective=s$objective,
          coordinate_count=length(fit$attempts[[fit$best]]$coordinates),
          AIC=NA_real_,BIC=NA_real_,
          seasonal_maxima=p22_fourier_maxima(s$parameters,setup$problem$basis))
        annual[[key]] <- p22_annual_fit_errors(s,inputs,nuisance,lambda,family)
      } else summary[[key]] <- data.frame(nuisance=nuisance,lambda=lambda,
        family=family,status=fit$status,historical_loglik=NA_real_,
        extension_loglik=NA_real_,raw_loglik=NA_real_,post_penalty=NA_real_,
        objective=NA_real_,coordinate_count=NA_integer_,AIC=NA_real_,
        BIC=NA_real_,seasonal_maxima=NA_integer_)
    }
  list(status="ATTEMPTS_SAVED",fits=records,summary=do.call(rbind,summary),
    annual_errors=if(length(annual)) do.call(rbind,annual) else data.frame(),
    input_hashes=inputs$hashes)
}

p22_annual_fit_errors <- function(score,inputs,nuisance,lambda,family) {
  year <- c(inputs$replay$data$time_markers$year,inputs$calendar$year)
  week <- c(inputs$replay$data$time_markers$week,inputs$calendar$source_week)
  do.call(rbind,lapply(sort(unique(year)),function(y) {
    ix <- which(year==y & score$mask)
    observed <- score$observed[ix];expected <- score$states$expected_cases[ix]
    timing <- function(v) week[ix][which(cumsum(v)>=.5*sum(v))[1L]]
    data.frame(nuisance=nuisance,lambda=lambda,family=family,year=y,
      period=if(y<=2019) "historical" else "extension",
      observed_total=sum(observed),expected_total=sum(expected),
      burden_error=sum(expected)-sum(observed),
      observed_peak_week=week[ix][which.max(observed)],
      expected_peak_week=week[ix][which.max(expected)],
      peak_week_error=week[ix][which.max(expected)]-week[ix][which.max(observed)],
      median_timing_error=timing(expected)-timing(observed))
  }))
}

p22_fourier_maxima <- function(parameters,basis) {
  x <- seq(0,52,length.out=5201L)[-5201L]
  y <- fourier_log_shape(x,parameters$free_coordinates,basis)
  sum(y>c(tail(y,1L),head(y,-1L)) & y>=c(tail(y,-1L),head(y,1L)))
}

p22_block_validation <- function(inputs,cfg,include_fourier=TRUE) {
  rows <- list();attempts <- list()
  p21_env <- new.env(parent=baseenv())
  sys.source("config/phase21.R",p21_env)
  pcfg <- p21_env$phase21_config()
  families <- if(include_fourier) c("FROZEN","JOINT","A","AS",
    "A_FIXED","AS_FIXED") else c("FROZEN","JOINT")
  for(k in seq_along(cfg$validation_blocks)) {
    held <- cfg$validation_blocks[[k]]
    p22_progress("validation block %d/%d: held future weeks %d-%d",
      k,length(cfg$validation_blocks),min(held),max(held))
    mask <- c(FALSE,rep(TRUE,363L),!seq_len(209L) %in% held)
    future_mask <- !seq_len(209L) %in% held
    for(family in families) {
      lambdas <- if(family %in% c("A","AS","A_FIXED","AS_FIXED"))
        cfg$q_smoothness else .1
      for(lambda in lambdas) {
        p22_progress("validation block %d/%d: %s lambda=%s start",
          k,length(cfg$validation_blocks),family,format(lambda))
        inp <- inputs;inp$candidate$lambda <- lambda
        key <- paste(k,family,lambda,sep="/")
        if(family=="FROZEN") {
          c21 <- inp$candidate;c21$lambda <- lambda
          f <- p21_fit(inp$ref,inp$calendar,c21,inp$basis,pcfg,future_mask,
            allow_e1=FALSE,seed=cfg$seed+k)
          if(!identical(f$status,"NUMERICALLY_VALID")) {
            p22_progress("validation block %d: %s lambda=%s failed: %s",
              k,family,format(lambda),f$status)
            attempts[[key]] <- f
            rows[[key]] <- data.frame(block=k,family=family,lambda=lambda,
              status=f$status,held_log_score=NA_real_,held_MAE=NA_real_,
              held_burden_error=NA_real_,held_median_timing_error=NA_real_)
            next
          }
          setup <- p22_problem(inp)
          coords <- inp$ref$record$retained_coordinates[setup$problem$map$order]
          score <- p22_score(coords,f$theta,inp,setup)
          attempts[[key]] <- f
        } else {
          setup <- p22_problem(inp,if(family=="JOINT") "SPL" else
            sub("_FIXED$","",family),fixed_nuisance=grepl("_FIXED$",family))
          f <- p22_fit(inp,setup,cfg,mask,seed=cfg$seed+100L*k+
            match(family,families),allow_full_start=FALSE)
          attempts[[key]] <- f
          if(!identical(f$status,"ENDPOINTS_RECORDED")) {
            p22_progress("validation block %d: %s lambda=%s failed: %s",
              k,family,format(lambda),f$status)
            rows[[key]] <- data.frame(block=k,family=family,lambda=lambda,
              status=f$status,held_log_score=NA_real_,held_MAE=NA_real_,
              held_burden_error=NA_real_,held_median_timing_error=NA_real_)
            next
          }
          score <- f$attempts[[f$best]]$score
        }
        idx <- 364L+held
        observed <- score$observed[idx];expected <- score$states$expected_cases[idx]
        median_at <- function(x) held[which(cumsum(x)>=.5*sum(x))[1L]]
        rows[[key]] <- data.frame(block=k,family=family,lambda=lambda,
          status="SCORED",held_log_score=sum(score$per_week[idx]),
          held_MAE=mean(abs(observed-expected)),
          held_burden_error=sum(expected)-sum(observed),
          held_median_timing_error=median_at(expected)-median_at(observed))
        p22_progress("validation block %d: %s lambda=%s held log score=%0.4f",
          k,family,format(lambda),sum(score$per_week[idx]))
      }
    }
  }
  scores <- do.call(rbind,rows)
  aggregate <- stats::aggregate(cbind(held_log_score,held_MAE,
    held_burden_error,held_median_timing_error)~family+lambda,
    scores,FUN=sum)
  list(status="BLOCKED_VALIDATION_SAVED",fold_scores=scores,
    aggregate=aggregate,attempts=attempts,
    policy="Only future rows held out; historical reference trained on all historical rows; no full-data fit used as a validation start")
}

p22_simulate_reference <- function(reference,inputs,u=0,v=0,
                                   origin="2013") {
  if(!origin %in% as.character(2013:2023)) stop("Unknown intervention origin")
  base <- p22_score(reference$coords,reference$post,inputs,reference$setup)
  p <- base$parameters;N <- inputs$replay$data$fixed_parameters$pop
  fixed <- list(pop=N,mu=ref_or(p$mu,inputs$replay$data$fixed_parameters$mu),
    gamma=ref_or(p$gamma,inputs$replay$data$fixed_parameters$gamma),
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  if(origin=="2013") {
    S <- p$S0*exp(u)
    if(S<=0 || S>=1-p$I0) stop("Susceptibility perturbation violates initial simplex")
    initial <- make_initial_state(S,p$I0,N)
    weeks <- c(inputs$replay$data$seasonal_week,inputs$calendar$seasonal_week)
    multiplier <- c(rep(exp(v),364L),base$q*exp(v))
    z <- simulate_weekly_sirs(base$beta,initial,fixed,weeks,multiplier)
  } else {
    hist_year <- c(inputs$replay$data$time_markers$year,inputs$calendar$year)
    cut <- max(which(hist_year < as.integer(origin)))
    prior <- base$states[seq_len(cut),,drop=FALSE]
    terminal <- unlist(prior[cut,c("S","I","R")],use.names=TRUE)
    terminal["S"] <- terminal["S"]*exp(u)
    terminal["R"] <- sum(unlist(prior[cut,c("S","I","R")],use.names=FALSE))-
      terminal["S"]-terminal["I"]
    if(any(!is.finite(terminal)) || any(terminal<=0)) stop("Pandemic-onset state infeasible")
    future <- simulate_weekly_sirs(base$beta,terminal,fixed,
      c(if(cut<=364L) inputs$replay$data$seasonal_week[cut] else inputs$calendar$seasonal_week[cut-364L],
        c(inputs$replay$data$seasonal_week,inputs$calendar$seasonal_week)[(cut+1L):length(hist_year)]),
      c(1,c(rep(1,364L),base$q)[(cut+1L):length(hist_year)]*exp(v)),
      initial_mass=sum(terminal))[-1L,,drop=FALSE]
    z <- rbind(prior,future)
  }
  hist <- inputs$replay$data$time_markers
  z$year <- c(hist$year,inputs$calendar$year)
  z$source_week <- c(hist$week,inputs$calendar$source_week)
  z$initialization_only <- seq_len(nrow(z))==1L
  z
}

p22_peak_record <- function(reference,inputs,u,v,id,origin="2013") {
  z <- p22_simulate_reference(reference,inputs,u,v,origin)
  z$week <- z$source_week
  list(id=id,a=exp(v),u=u,v=v,trajectory=z,status="VALID")
}

p22_mechanism_grid <- function(cfg) {
  z <- utils::read.csv(cfg$mechanism_grid)
  grid <- sort(unique(c(z$left_v[z$year %in% 2018:2019],
    z$right_v[z$year %in% 2018:2019])))
  if(length(grid)!=44L || abs(grid[1]+log(1.2))>1e-12 || tail(grid,1L)!=0)
    stop("Historical transmission grid differs from paper Table 1")
  grid
}

p22_assign_branches <- function(peaks,matches,grid) {
  if(!nrow(peaks)) return(peaks)
  peaks$branch <- NA_character_
  for(name in unique(peaks$reference)) for(year in unique(peaks$year[peaks$reference==name])) {
    previous <- character();next_id <- 0L
    for(k in seq_along(grid)) {
      ix <- which(peaks$reference==name & peaks$year==year &
        abs(peaks$v-grid[k])<1e-12 & peaks$genuine)
      if(!length(ix)) {previous <- character();next}
      current <- setNames(rep(NA_character_,length(ix)),as.character(peaks$peak_id[ix]))
      if(k>1L && nrow(matches)) {
        edge <- matches[matches$reference==name & matches$year==year &
          abs(matches$v_left-grid[k-1L])<1e-12 &
          abs(matches$v_right-grid[k])<1e-12 &
          matches$status=="MATCHED" & matches$direction_consistent,,drop=FALSE]
        for(j in seq_len(nrow(edge))) {
          value <- previous[as.character(edge$left_peak[j])]
          if(length(value)==1L && !is.na(value))
            current[as.character(edge$right_peak[j])] <- value
        }
      }
      for(j in which(is.na(current))) {
        next_id <- next_id+1L
        current[j] <- paste(name,year,"branch",next_id,sep="/")
      }
      peaks$branch[ix] <- unname(current[as.character(peaks$peak_id[ix])])
      previous <- current
    }
  }
  peaks
}

p22_mechanism <- function(inputs,cfg,reference_names=c("HIST","FROZEN","JOINT"),
                          grid=NULL,refine=TRUE,year_scope=NULL) {
  refs <- p22_references(inputs,cfg)
  fit <- p22_saved_joint(inputs,cfg)$fit$primary
  d <- length(refs$JOINT$coords)
  include_competitive <- identical(reference_names,c("HIST","FROZEN","JOINT"))
  for(i in setdiff(fit$competitive,fit$best)) {
    a <- fit$attempts[[i]]
    key <- paste0("JOINT_COMPETITIVE_",i)
    refs[[key]] <- list(coords=a$coordinates[seq_len(d)],
      post=tail(a$coordinates,inputs$basis$dimension),setup=refs$JOINT$setup)
    if(include_competitive) reference_names <- c(reference_names,key)
  }
  if(is.null(grid)) grid <- p22_mechanism_grid(cfg)
  policy <- p10_defaults()
  rows <- peaks_out <- matches_out <- windows_out <- padded_out <- list()
  for(name in reference_names) {
    p22_progress("mechanism %s: %d transmission values, %d-%d windows",
      name,length(grid),if(name=="HIST") 2018 else 2018,
      if(name=="HIST") 2019 else 2023)
    years <- if(is.null(year_scope)) { if(name=="HIST") 2018:2019 else 2018:2023 } else if(!is.null(year_scope[[name]])) year_scope[[name]] else year_scope$JOINT
    if(!length(years)) next
    records <- lapply(seq_along(grid),function(k) {
      if(k==1L || k%%5L==0L || k==length(grid))
        p22_progress("mechanism %s: trajectory %d/%d",name,k,length(grid))
      p22_peak_record(refs[[name]],inputs,0,grid[k],sprintf("%s-%03d",name,k))
    })
    P <- W <- vector("list",length(records))
    for(k in seq_along(records)) {
      P[[k]] <- W[[k]] <- list()
      for(y in years) {
        record <- records[[k]]
        if(y<=2019L) record$trajectory <- record$trajectory[1:364,,drop=FALSE]
        pk <- p10_detect(record,y,policy$peaks,policy$identity)
        P[[k]][[as.character(y)]] <- pk
        W[[k]][[as.character(y)]] <- p10_window(record,y,pk,policy$peaks)
        ix <- which(record$trajectory$year==y & !record$trajectory$initialization_only)
        pad_left <- max(2L,min(ix)-cfg$mechanism$padding_weeks)
        pad_right <- min(nrow(record$trajectory),max(ix)+cfg$mechanism$padding_weeks)
        raw <- p4_weekly_peaks(record$trajectory$incidence[-1L],
          2:nrow(record$trajectory),policy$peaks)$local
        padded_out[[length(padded_out)+1L]] <- data.frame(reference=name,
          year=y,v=grid[k],calendar_start=min(ix),calendar_end=max(ix),
          padded_start=pad_left,padded_end=pad_right,
          chronological_local_peaks=sum(raw$first>=pad_left &
            raw$first<=pad_right & !raw$plateau),
          missing_right_support=max(ix)+cfg$mechanism$padding_weeks>
            nrow(record$trajectory))
        if(nrow(pk)) {
          pk$reference <- name;pk$v <- grid[k]
          peaks_out[[length(peaks_out)+1L]] <- pk
        }
        w <- W[[k]][[as.character(y)]];w$reference <- name;w$v <- grid[k]
        windows_out[[length(windows_out)+1L]] <- w
      }
    }
    for(y in years) for(k in seq_len(length(grid)-1L)) {
      left <- P[[k]][[as.character(y)]];right <- P[[k+1L]][[as.character(y)]]
      wl <- W[[k]][[as.character(y)]];wr <- W[[k+1L]][[as.character(y)]]
      forward <- p10_assignment(left,right,policy$identity)
      reverse <- p10_assignment(left,right,policy$identity,TRUE)
      consistent <- p10_identity_consistent(forward,reverse)
      if(nrow(forward)) {
        forward$reference <- name;forward$year <- y
        forward$v_left <- grid[k];forward$v_right <- grid[k+1L]
        forward$direction_consistent <- consistent
        matches_out[[length(matches_out)+1L]] <- forward
      }
      matched <- forward[forward$status=="MATCHED",,drop=FALSE]
      left_winner <- wl$peak_id;right_winner <- wr$peak_id
      same <- consistent && !is.na(left_winner) && !is.na(right_winner) &&
        any(matched$left_peak==left_winner & matched$right_peak==right_winner)
      gap_left <- gap_right <- NA_real_;exchange <- FALSE
      if(consistent && nrow(matched)>=2L && !is.na(left_winner) && !is.na(right_winner)) {
        a <- which(matched$left_peak==left_winner)
        b <- which(matched$right_peak==right_winner)
        if(length(a)==1L) {
          rivals <- setdiff(seq_len(nrow(matched)),a)
          rival_height <- left$height[match(matched$left_peak[rivals],left$peak_id)]
          competitor <- rivals[which.max(rival_height)]
          la <- left$height[match(matched$left_peak[c(a,competitor)],left$peak_id)]
          ra <- right$height[match(matched$right_peak[c(a,competitor)],right$peak_id)]
          gap_left <- la[1]-la[2];gap_right <- ra[1]-ra[2]
          exchange <- length(b)==1L && b==competitor &&
            is.finite(gap_left) && is.finite(gap_right) && gap_left*gap_right<0
        }
      }
      boundary <- wl$status %in% c("WINDOW_SELECTION","SUPPORT_BOUNDARY") ||
        wr$status %in% c("WINDOW_SELECTION","SUPPORT_BOUNDARY")
      suppression <- !nrow(left) || !nrow(right) ||
        !any(left$genuine) || !any(right$genuine)
      timing_change <- wl$position!=wr$position
      verified_exchange <- exchange &&
        grid[k+1L]-grid[k]<=cfg$mechanism$min_log_width
      cls <- if(suppression) "NO_ELIGIBLE_PEAK" else if(boundary) "WINDOW_SELECTION" else
        if(verified_exchange) "DOMINANCE_EXCHANGE" else if(exchange) "UNRESOLVED" else if(same && timing_change)
          "WITHIN_BRANCH_DISPLACEMENT" else if(same && !timing_change)
            "STABLE_SAME_BRANCH" else "UNRESOLVED"
      rows[[length(rows)+1L]] <- data.frame(reference=name,year=y,
        v_left=grid[k],v_right=grid[k+1L],left_week=wl$week,
        right_week=wr$week,left_peak=left_winner,right_peak=right_winner,
        left_status=wl$status,right_status=wr$status,
        signed_competitor_gap_left=gap_left,signed_competitor_gap_right=gap_right,
        persistent_winner=same,crossing_candidate=exchange,
        identity_consistent=consistent,
        boundary_selection=boundary,right_support_missing=y==2023,
        timing_change=timing_change,classification=cls)
    }
    p22_progress("mechanism %s: annual intervals classified",name)
  }
  intervals <- do.call(rbind,rows)
  replay_difference <- data.frame()
  if("HIST" %in% reference_names && length(grid)==44L &&
     all(2018:2019 %in% unique(intervals$year[intervals$reference=="HIST"]))) {
    paper <- utils::read.csv(cfg$historical_mechanism)
    translation <- c(DISPLACEMENT_SUPPORTED="WITHIN_BRANCH_DISPLACEMENT",
      STABLE_SAME_BRANCH="STABLE_SAME_BRANCH",
      WINDOW_SELECTION_SUPPORTED="WINDOW_SELECTION",UNRESOLVED="UNRESOLVED")
    h <- intervals[intervals$reference=="HIST" & intervals$year %in% 2018:2019,]
    h$left_key <- round(h$v_left,10);h$right_key <- round(h$v_right,10)
    paper$left_key <- round(paper$left_v,10)
    paper$right_key <- round(paper$right_v,10)
    replay_difference <- merge(h[,c("year","left_key","right_key","classification")],
      paper[,c("year","left_key","right_key","revised_classification")],
      by=c("year","left_key","right_key"))
    if(nrow(replay_difference)!=86L) stop("Historical mechanism interval join is incomplete")
    replay_difference$paper_classification <- unname(translation[
      replay_difference$revised_classification])
    replay_difference$matches <- replay_difference$classification==
      replay_difference$paper_classification
    replay_difference$reason <- ifelse(replay_difference$matches,"MATCHED",
      "Phase 10 chronological classifier versus saved Phase 14 adjudication; compare boundary and identity evidence")
  }
  if(all(c("HIST","FROZEN") %in% reference_names) &&
     is.null(year_scope)) {
    fields <- c("year","v_left","v_right","left_week","right_week",
      "left_status","right_status","classification")
    h <- intervals[intervals$reference=="HIST" & intervals$year<=2019,fields]
    f <- intervals[intervals$reference=="FROZEN" & intervals$year<=2019,fields]
    rownames(h) <- rownames(f) <- NULL
    if(!identical(h,f)) stop("FROZEN historical-window invariance failed")
  }
  # A conservative audit never promotes an unresolved edge to an exchange.
  counts <- as.data.frame(table(intervals$reference,intervals$year,
    intervals$classification),stringsAsFactors=FALSE)
  names(counts) <- c("reference","year","classification","count")
  counts <- counts[counts$count>0,]
  coverage <- do.call(rbind,lapply(split(intervals,list(intervals$reference,
    intervals$year),drop=TRUE),function(x)
      data.frame(reference=x$reference[1],year=x$year[1],all_intervals=nrow(x),
        resolved_timing_changes=sum(x$classification %in%
          c("WITHIN_BRANCH_DISPLACEMENT","WINDOW_SELECTION","DOMINANCE_EXCHANGE")),
        unresolved=sum(x$classification=="UNRESOLVED"),
        no_eligible_peak=sum(x$classification=="NO_ELIGIBLE_PEAK"),
        missing_right_support=any(x$right_support_missing))))
  peak_table <- if(length(peaks_out)) do.call(rbind,peaks_out) else data.frame()
  match_table <- if(length(matches_out)) do.call(rbind,matches_out) else data.frame()
  peak_table <- p22_assign_branches(peak_table,match_table,grid)
  result <- list(status="MECHANISM_AUDITED",intervals=intervals,counts=counts,
    coverage=coverage,peaks=peak_table,
    windows=do.call(rbind,windows_out),
    padded_chronological=do.call(rbind,padded_out),
    matches=match_table,
    historical_classifier_difference=replay_difference,
    refinement_history=data.frame(),
    limitation="Unresolved cells require bracket refinement before exchange verdict")
  if(refine && cfg$mechanism$refine_levels>0L) {
    history <- list();used <- 0L
    initial_replay_difference <- result$historical_classifier_difference
    for(level in seq_len(cfg$mechanism$refine_levels)) {
      x <- result$intervals
      x <- x[x$classification %in% c("UNRESOLVED","DOMINANCE_EXCHANGE",
        "WINDOW_SELECTION") & x$v_right-x$v_left>cfg$mechanism$min_log_width,]
      if(!nrow(x)) break
      x$level <- level;x$midpoint <- (x$v_left+x$v_right)/2
      pair <- unique(x[c("v_left","v_right","midpoint")])
      if(used+nrow(pair)>cfg$mechanism$max_new_points)
        pair <- head(pair,max(0L,cfg$mechanism$max_new_points-used))
      if(!nrow(pair)) break
      p22_progress("mechanism refinement level %d/%d: %d new transmission values",
        level,cfg$mechanism$refine_levels,nrow(pair))
      selected <- paste(pair$v_left,pair$v_right)
      x$selected <- paste(x$v_left,x$v_right) %in% selected
      history[[length(history)+1L]] <- x[,c("level","reference","year",
        "v_left","v_right","midpoint","classification","selected")]
      grid <- sort(unique(c(grid,pair$midpoint)))
      used <- used+nrow(pair)
      result <- p22_mechanism(inputs,cfg,reference_names,grid,refine=FALSE,year_scope=year_scope)
    }
    result$refinement_history <- if(length(history)) do.call(rbind,history) else data.frame()
    result$historical_classifier_difference <- initial_replay_difference
  }
  result
}

p22_branch_peak <- function(base,other,year,policy,base_peak) {
  if(year<=2019L) {
    base$trajectory <- base$trajectory[1:364,,drop=FALSE]
    other$trajectory <- other$trajectory[1:364,,drop=FALSE]
  }
  P <- p10_detect(other,year,policy$peaks,policy$identity)
  B <- p10_detect(base,year,policy$peaks,policy$identity)
  if(!nrow(P) || !nrow(B) || is.na(base_peak))
    return(list(status="NO_ELIGIBLE_PEAK",time=NA_real_,fixed=NA_real_,
      stencil_changed=NA))
  f <- p10_assignment(B,P,policy$identity)
  r <- p10_assignment(B,P,policy$identity,TRUE)
  match <- f[f$status=="MATCHED" & f$left_peak==base_peak,,drop=FALSE]
  if(nrow(match)!=1L || !p10_identity_consistent(f,r))
    return(list(status="UNMATCHED_BRANCH",time=NA_real_,fixed=NA_real_,
      stencil_changed=NA))
  peak <- P[P$peak_id==match$right_peak,]
  if(nrow(peak)!=1L || !peak$genuine || peak$quadratic_status!="OK")
    return(list(status="CURVATURE_OR_BOUNDARY",time=NA_real_,fixed=NA_real_,
      stencil_changed=NA))
  z <- other$trajectory
  fixed <- p10_interpolate(z$incidence,base_peak,base_peak,policy$identity)
  list(status="ELIGIBLE",time=peak$quadratic_position,
    fixed=if(fixed$status=="OK") fixed$position else NA_real_,
    stencil_changed=peak$peak_id!=base_peak,
    peak_id=peak$peak_id,height=peak$height,curvature=peak$curvature)
}

p22_sensitivity <- function(inputs,cfg,year_scope=NULL,origin="2013") {
  refs <- p22_references(inputs,cfg)
  fit <- p22_saved_joint(inputs,cfg)$fit$primary
  d <- length(refs$JOINT$coords)
  for(i in setdiff(fit$competitive,fit$best)) {
    a <- fit$attempts[[i]]
    refs[[paste0("JOINT_COMPETITIVE_",i)]] <- list(
      coords=a$coordinates[seq_len(d)],
      post=tail(a$coordinates,inputs$basis$dimension),setup=refs$JOINT$setup)
  }
  policy <- p10_defaults();rows <- list();details <- list()
  for(name in names(refs)) {
    p22_progress("sensitivity %s: %d symmetric steps",name,
      length(cfg$sensitivity_steps))
    years <- if(is.null(year_scope)) { if(name=="HIST") 2018:2019 else 2018:2023 } else if(!is.null(year_scope[[name]])) year_scope[[name]] else year_scope$JOINT
    if(!length(years)) next
    base <- p22_peak_record(refs[[name]],inputs,0,0,paste0(name,"/base"),origin)
    baseline <- lapply(years,function(y) {
      record <- base
      if(y<=2019L) record$trajectory <- record$trajectory[1:364,,drop=FALSE]
      B <- p10_detect(record,y,policy$peaks,policy$identity)
      W <- p10_window(record,y,B,policy$peaks)
      list(peaks=B,window=W)
    });names(baseline) <- as.character(years)
    for(h in cfg$sensitivity_steps) {
      p22_progress("sensitivity %s: h=%0.8f",name,h)
      variants <- list()
      for(axis in c("S","beta")) for(sign in c(-1,1)) {
        key <- paste(axis,sign,sep="/")
        variants[[key]] <- tryCatch(p22_peak_record(refs[[name]],inputs,
          if(axis=="S") sign*h else 0,if(axis=="beta") sign*h else 0,
          paste(name,axis,sign,h,sep="/"),origin),error=function(e)e)
      }
      for(y in years) {
        b <- baseline[[as.character(y)]]
        ref_peak <- b$window$peak_id
        values <- list()
        for(axis in c("S","beta")) for(sign in c(-1,1)) {
          key <- paste(axis,sign,sep="/")
          values[[key]] <- if(inherits(variants[[key]],"error"))
            list(status=paste0("INFEASIBLE: ",conditionMessage(variants[[key]])),
              time=NA_real_,fixed=NA_real_,stencil_changed=NA) else
            p22_branch_peak(base,variants[[key]],y,policy,ref_peak)
          details[[length(details)+1L]] <- data.frame(reference=name,year=y,h=h,
            axis=axis,sign=sign,status=values[[key]]$status,
            branch_time=values[[key]]$time,fixed_stencil_time=values[[key]]$fixed,
            stencil_changed=values[[key]]$stencil_changed)
        }
        status <- vapply(values,function(v) v$status,character(1))
        paired <- all(status=="ELIGIBLE") &&
          b$window$status=="GENUINE_LOCAL_PEAK"
        deriv <- function(axis,field) {
          plus <- values[[paste0(axis,"/1")]][[field]]
          minus <- values[[paste0(axis,"/-1")]][[field]]
          if(paired && is.finite(plus) && is.finite(minus))
            (plus-minus)/(2*h) else NA_real_
        }
        es <- deriv("S","time");eb <- deriv("beta","time")
        fs <- deriv("S","fixed");fb <- deriv("beta","fixed")
        rows[[length(rows)+1L]] <- data.frame(reference=name,year=y,h=h,
          estimand="2013_INITIAL_S_AND_COMPLETE_BETA",
          paired_eligible=paired,
          exclusion=if(paired) "" else paste(unique(c(b$window$status,status[status!="ELIGIBLE"])),collapse=";"),
          E_S=es,E_beta=eb,fixed_E_S=fs,fixed_E_beta=fb,
          stencil_changed=any(vapply(values,function(v)isTRUE(v$stencil_changed),logical(1))),
          signed_difference=eb-es,absolute_difference=abs(eb)-abs(es),
          ratio=if(is.finite(es) && abs(es)>1e-4) abs(eb)/abs(es) else NA_real_,
          sign_ordering=is.finite(es) && is.finite(eb) && es<0 && eb<0,
          magnitude_ordering=is.finite(es) && is.finite(eb) && abs(eb)>abs(es))
      }
    }
  }
  list(status="PAIRED_SENSITIVITY_SAVED",estimates=do.call(rbind,rows),
    variants=do.call(rbind,details))
}

p22_annual <- function(inputs,cfg) {
  replay <- p22_historical_replay(cfg)
  profiles <- utils::read.csv(cfg$annual_profiles)
  profiles <- profiles[profiles$profile_id %in% c("profile-004","profile-009"),]
  scenarios <- utils::read.csv(cfg$annual_scenarios)
  scenarios <- scenarios[scenarios$profile_id %in% profiles$profile_id &
    scenarios$axis %in% c("baseline","beta") & scenarios$u==0,]
  refs <- p22_references(inputs,cfg)
  joint <- p22_score(refs$JOINT$coords,refs$JOINT$post,inputs,refs$JOINT$setup)
  p <- joint$parameters;N <- inputs$replay$data$fixed_parameters$pop
  fixed <- list(pop=N,mu=ref_or(p$mu,inputs$replay$data$fixed_parameters$mu),
    gamma=ref_or(p$gamma,inputs$replay$data$fixed_parameters$gamma),
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  out <- list();records <- list();peak_list <- list();window_list <- list()
  policy <- p10_defaults()
  for(i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i,];profile <- profiles[profiles$profile_id==sc$profile_id,]
    p22_progress("annual transfer %d/%d: %s a=%0.6f",
      i,nrow(scenarios),sc$profile_id,sc$a)
    beta <- sc$a*profile$b*exp(profile$amplitude*
      cos(2*pi*((1:52)-profile$phase_week)/52))
    n <- 1L+52L*4L;week <- ((seq_len(n)-2L)%%52L)+1L
    z <- simulate_weekly_sirs(beta,make_initial_state(p$S0,p$I0,N),
      fixed,week,rep(1,n))
    z$year <- c(0L,rep(1:4,each=52L))
    z$initialization_only <- seq_len(n)==1L
    record <- list(id=sc$scenario_id,a=sc$a,trajectory=z)
    records[[sc$scenario_id]] <- record
    peak_list[[sc$scenario_id]] <- p10_detect(record,4L,policy$peaks,policy$identity)
    window_list[[sc$scenario_id]] <- p10_window(record,4L,
      peak_list[[sc$scenario_id]],policy$peaks)
    wave <- z$incidence[(n-51L):n]
    peak <- which.max(wave)
    out[[length(out)+1L]] <- data.frame(scenario=sc$scenario_id,
      profile_id=sc$profile_id,a=sc$a,
      year=4L,peak_week=peak,peak_height=wave[peak],
      continuous_forcing_maxima=1L,forcing_q="1 throughout",
      nuisance="JOINT",status=if(peak %in% c(1L,52L))
        "WINDOW_SELECTION" else "INTERIOR_PEAK")
  }
  joint <- do.call(rbind,out)
  pair <- do.call(rbind,lapply(split(joint,joint$profile_id),function(z) {
    z <- z[order(z$a),]
    if(nrow(z)<2L) return(NULL)
    do.call(rbind,lapply(seq_len(nrow(z)-1L),function(k) {
      left <- peak_list[[z$scenario[k]]]
      right <- peak_list[[z$scenario[k+1L]]]
      wl <- window_list[[z$scenario[k]]]
      wr <- window_list[[z$scenario[k+1L]]]
      f <- p10_assignment(left,right,policy$identity)
      rev <- p10_assignment(left,right,policy$identity,TRUE)
      consistent <- p10_identity_consistent(f,rev)
      m <- f[f$status=="MATCHED",,drop=FALSE]
      same <- consistent && !is.na(wl$peak_id) && !is.na(wr$peak_id) &&
        any(m$left_peak==wl$peak_id & m$right_peak==wr$peak_id)
      gap_left <- gap_right <- NA_real_;exchange <- FALSE
      a <- which(m$left_peak==wl$peak_id)
      b <- which(m$right_peak==wr$peak_id)
      if(consistent && nrow(m)>=2L && length(a)==1L && length(b)==1L && a!=b) {
        lh <- left$height[match(m$left_peak[c(a,b)],left$peak_id)]
        rh <- right$height[match(m$right_peak[c(a,b)],right$peak_id)]
        gap_left <- lh[1]-lh[2];gap_right <- rh[1]-rh[2]
        exchange <- is.finite(gap_left) && is.finite(gap_right) &&
          gap_left*gap_right<0
      }
      boundary <- wl$status=="WINDOW_SELECTION" || wr$status=="WINDOW_SELECTION"
      class <- if(boundary) "WINDOW_SELECTION" else if(exchange)
        "DOMINANCE_EXCHANGE" else if(same && wl$week!=wr$week)
          "FINITE_WITHIN_BRANCH_DISPLACEMENT" else if(same)
            "STABLE_SAME_BRANCH" else "UNRESOLVED"
      data.frame(profile_id=z$profile_id[k],left_scenario=z$scenario[k],
        right_scenario=z$scenario[k+1L],a_left=z$a[k],a_right=z$a[k+1L],
        left_week=wl$week,right_week=wr$week,
        finite_week_change=wr$week-wl$week,log_width=log(z$a[k+1L]/z$a[k]),
        same_branch=same,signed_gap_left=gap_left,signed_gap_right=gap_right,
        window_selection=boundary,abrupt_candidate=abs(wr$week-wl$week)>=8 &&
          log(z$a[k+1L]/z$a[k])<=.05,
        dominance_exchange_verified=exchange,classification=class)
    }))
  }))
  list(status="ANNUAL_COUNTERFACTUAL_SAVED",historical=replay$annual,
    joint=joint,pairs=pair,
    interpretation="Finite displacement across sampled parameter paths; no abruptness inferred")
}

p22_claims <- function(replay,mechanism,sensitivity,annual,fourier,validation) {
  m <- mechanism$intervals
  jm <- m[m$reference=="JOINT" & m$year %in% 2018:2019,]
  post <- m[m$reference %in% c("FROZEN","JOINT") & m$year>=2020,]
  post_by_reference <- vapply(c("FROZEN","JOINT"),function(name) {
    z <- post[post$reference==name,]
    sprintf("%s: %d timing changes, %d unresolved/suppressed",
      name,sum(z$classification %in% c("WITHIN_BRANCH_DISPLACEMENT",
        "WINDOW_SELECTION","DOMINANCE_EXCHANGE")),
      sum(z$classification %in% c("UNRESOLVED","NO_ELIGIBLE_PEAK")))
  },character(1))
  timing <- function(x) sum(x$classification %in%
    c("WITHIN_BRANCH_DISPLACEMENT","WINDOW_SELECTION","DOMINANCE_EXCHANGE"))
  unresolved <- function(x) sum(x$classification %in% c("UNRESOLVED","NO_ELIGIBLE_PEAK"))
  jm_verdict <- if(timing(jm)==0 || unresolved(jm)>0) "INCONCLUSIVE" else if(
    any(jm$classification=="DOMINANCE_EXCHANGE")) "MIXED" else "SUPPORTED"
  replay_difference <- mechanism$historical_classifier_difference
  if(nrow(replay_difference) && any(!replay_difference$matches))
    jm_verdict <- "INCONCLUSIVE"
  competitive <- m[grepl("^JOINT_COMPETITIVE_",m$reference) &
    m$year %in% 2018:2019,]
  if(nrow(competitive) && jm_verdict=="SUPPORTED" &&
     (any(competitive$classification=="DOMINANCE_EXCHANGE") ||
      unresolved(competitive)>0)) jm_verdict <- "MIXED"
  s <- sensitivity$estimates
  js <- s[s$reference=="JOINT" & s$year %in% 2018:2019,]
  ps <- s[s$reference %in% c("FROZEN","JOINT") & s$year>=2020,]
  post_sensitivity <- vapply(c("FROZEN","JOINT"),function(name) {
    z <- ps[ps$reference==name,]
    sprintf("%s: %d/%d paired eligible",name,sum(z$paired_eligible),nrow(z))
  },character(1))
  judge_s <- function(x) {
    eligible <- x[x$paired_eligible,]
    if(nrow(eligible)<4L || length(unique(eligible$year))<2L) return("INCONCLUSIVE")
    good <- eligible$sign_ordering & eligible$magnitude_ordering
    if(all(good)) "SUPPORTED" else if(any(good)) "MIXED" else "NOT_SUPPORTED"
  }
  sensitivity_verdict <- judge_s(js)
  competitive_s <- s[grepl("^JOINT_COMPETITIVE_",s$reference) &
    s$year %in% 2018:2019,]
  if(nrow(competitive_s) && sensitivity_verdict=="SUPPORTED") {
    competitive_verdicts <- vapply(split(competitive_s,
      competitive_s$reference),judge_s,character(1))
    if(any(competitive_verdicts=="INCONCLUSIVE"))
      sensitivity_verdict <- "INCONCLUSIVE" else
    if(any(competitive_verdicts!="SUPPORTED")) sensitivity_verdict <- "MIXED"
  }
  a <- annual$joint
  displacement <- vapply(split(a,a$profile_id),function(z) {
    z <- z[order(z$a),]
    max(z$peak_week)-min(z$peak_week)>=8L
  },logical(1))
  annual_verdict <- if(all(displacement)) "SUPPORTED" else "MIXED"
  if(annual_verdict=="SUPPORTED" && any(annual$pairs$classification %in%
    c("UNRESOLVED","DOMINANCE_EXCHANGE"))) annual_verdict <- "MIXED"
  v <- validation$aggregate
  at <- function(family,lambda=.1) v$held_log_score[v$family==family & v$lambda==lambda]
  a_score <- at("A");as_score <- at("AS")
  if(length(a_score)!=1L || length(as_score)!=1L ||
     !all(is.finite(c(a_score,as_score)))) fourier_verdict <- "INCONCLUSIVE" else {
    diff <- as_score-a_score
    fourier_verdict <- if(abs(diff)<=.02*60) "MIXED" else if(diff>0)
      "SUPPORTED" else "NOT_SUPPORTED"
    sm <- v[v$family %in% c("A","AS"),]
    by_lambda <- vapply(unique(sm$lambda),function(lambda) {
      z <- sm[sm$lambda==lambda,]
      left <- z$held_log_score[z$family=="A"]
      right <- z$held_log_score[z$family=="AS"]
      if(length(left)!=1L || length(right)!=1L) NA_real_ else right-left
    },numeric(1))
    if(anyNA(by_lambda)) fourier_verdict <- "INCONCLUSIVE" else
    if(fourier_verdict=="SUPPORTED" && any(by_lambda<=.02*60))
      fourier_verdict <- "MIXED"
  }
  data.frame(claim=c("mechanism_classification","sensitivity_sign_and_magnitude",
    "annual_only_finite_displacement","A_versus_AS_fit_advantage"),
    historical_replay=c("Table 1: 14/12 displacements, 0 exchanges",
      "Table 2: eight eligible pairs, both signed sensitivities negative; beta larger",
      "Profiles 004/009: nine/eighteen week finite displacement",
      "Table 3: AS higher unpenalized historical log likelihood"),
    same_window_JOINT=c(sprintf("%d timing changes; %d unresolved/suppressed",
      timing(jm),unresolved(jm)),
      sprintf("%d/%d paired eligible",sum(js$paired_eligible),nrow(js)),
      paste(names(displacement),displacement,collapse=";"),
      "Full-period Fourier comparison; see held-out scores"),
    post_period=c(paste(post_by_reference,collapse="; "),
      paste(post_sensitivity,collapse="; "),
      "Synthetic year 4 only; no post-period annual-only estimand",
      "2020-2023 held-out blocks included in full-period comparison"),
    eligibility_coverage=c(sprintf("%d JOINT intervals, %d unresolved/suppressed",
      nrow(jm),unresolved(jm)),
      sprintf("%d JOINT year-step pairs",sum(js$paired_eligible)),
      sprintf("%d saved paths transferred",sum(displacement)),
      sprintf("%d scored validation folds",sum(validation$fold_scores$status=="SCORED"))),
    verdict=c(jm_verdict,sensitivity_verdict,annual_verdict,fourier_verdict),
    post_verdict=c(if(timing(post)==0 || unresolved(post)>0) "INCONCLUSIVE" else if(any(post$classification=="DOMINANCE_EXCHANGE")) "MIXED" else "SUPPORTED",
      judge_s(ps),"NOT_APPLICABLE",fourier_verdict))
}

p22_figures <- function(replay,mechanism,sensitivity,annual,out) {
  dir.create(file.path(out,"figures"),showWarnings=FALSE)
  grDevices::pdf(file.path(out,"figures","figure1_mechanism.pdf"),width=8,height=5)
  tab <- xtabs(count~year+classification,mechanism$counts[
    mechanism$counts$reference=="JOINT" & mechanism$counts$year %in% 2018:2019,])
  barplot(t(tab),beside=TRUE,legend.text=TRUE,args.legend=list(x="topright",cex=.65),
    main="JOINT transmission path: audited interval classes",xlab="Evaluation year",
    ylab="Intervals")
  grDevices::dev.off()
  grDevices::pdf(file.path(out,"figures","figure2_sensitivity.pdf"),width=8,height=5)
  x <- sensitivity$estimates
  x <- x[x$reference=="JOINT" & x$year %in% 2018:2019 & x$paired_eligible,]
  if(nrow(x)) {
    plot(x$h,abs(x$E_S),log="xy",pch=16,col="steelblue",
      xlab="Symmetric log step h",ylab="Absolute weeks per log unit",
      main="JOINT paired peak timing sensitivity")
    points(x$h,abs(x$E_beta),pch=17,col="firebrick")
    legend("topright",c("S0","beta"),pch=c(16,17),
      col=c("steelblue","firebrick"),bty="n")
  } else plot.new()
  grDevices::dev.off()
  grDevices::pdf(file.path(out,"figures","figure3_annual.pdf"),width=8,height=5)
  x <- annual$joint
  plot(range(x$a),range(x$peak_week),type="n",xlab="Transmission multiplier a",
    ylab="Simulation year 4 peak week",main="Annual-only JOINT nuisance counterfactual")
  for(i in seq_along(unique(x$profile_id))) {
    z <- x[x$profile_id==unique(x$profile_id)[i],]
    z <- z[order(z$a),]
    lines(z$a,z$peak_week,type="b",pch=i,col=i)
  }
  legend("topright",unique(x$profile_id),pch=seq_along(unique(x$profile_id)),
    col=seq_along(unique(x$profile_id)),bty="n")
  grDevices::dev.off()
}
