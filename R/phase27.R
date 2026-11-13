# Phase 27 adapters retain the saved Phase 26 objective and the Phase 24/25 peak rules.
p27_write <- function(x,out,name) {
  utils::write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
  invisible(x)
}
p27_bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z)>0L,x)
  if(length(x)) do.call(rbind,x) else data.frame()
}
p27_hashes <- function(paths) {
  if(any(!file.exists(paths))) stop("Missing input: ",paste(paths[!file.exists(paths)],collapse=", "))
  h <- tools::md5sum(paths)
  if(anyNA(h)) stop("Input hash failed")
  data.frame(path=names(h),md5=unname(h),stringsAsFactors=FALSE)
}
p27_parameter_hash <- function(x) {
  p <- tempfile("phase27-vector-")
  on.exit(unlink(p),add=TRUE)
  saveRDS(unname(x),p,version=2)
  unname(tools::md5sum(p))
}
p27_saved <- function(path,stage) {
  x <- read_phase2_result(path)
  if(!identical(x$metadata$phase,26L) || !identical(x$metadata$stage,stage) ||
     !identical(x$result$config_version,"phase26-1")) stop("Wrong M2 artifact: ",path)
  h <- tools::md5sum(x$metadata$input_paths)
  if(anyNA(h) || !identical(unname(h),unname(x$metadata$input_hashes)))
    stop("M2 input dependency changed: ",path)
  x$result
}
p27_reference <- function(inputs,source_cfg,c26,c27) {
  fit <- p27_saved(c27$fit,"m2-fit")
  ev <- p27_saved(c27$evaluation,"m2-evaluate")
  if(fit$model!="M2" || ev$model!="M2" || !identical(ev$status,"EVALUATED"))
    stop("The saved endpoint is not evaluated M2")
  setup <- p22_problem(inputs)
  now <- p26_score(fit$endpoint,inputs,setup,c26,"M2",fit$penalty)
  old <- ev$score
  fields <- c("S","I","R","incidence","expected_cases")
  checks <- data.frame(check=c("parameter_vector","trajectory","expected_counts",
    "historical_level","historical_shape","future_q","historical_loglik",
    "future_loglik","seasonal_penalty","future_penalty","level_penalty",
    "shape_penalty","objective","weekly_effective_beta","likelihood_mask",
    "peak_2015","peak_2017","cog_2015","cog_2017",
    "prepared_week_53"),error=NA_real_)
  put <- function(key,value) checks$error[checks$check==key] <<- value
  put("parameter_vector",max(abs(fit$endpoint-
    head(fit$parameter_table$value,length(fit$endpoint)))))
  put("trajectory",max(abs(as.matrix(now$states[,fields])-as.matrix(old$states[,fields]))))
  put("expected_counts",max(abs(now$weekly$expected-old$weekly$expected)))
  put("historical_level",max(abs(now$weekly$log_level-old$weekly$log_level)))
  put("historical_shape",max(abs(now$weekly$log_shape-old$weekly$log_shape)))
  put("future_q",max(abs(now$weekly$future_q-old$weekly$future_q)))
  for(k in c("historical_loglik","future_loglik","seasonal_penalty","future_penalty",
             "level_penalty","shape_penalty","objective")) put(k,abs(now[[k]]-old[[k]]))
  put("weekly_effective_beta",max(abs(now$weekly$effective_beta-old$weekly$effective_beta)))
  put("likelihood_mask",sum(now$mask!=old$mask)+
    sum(now$mask!=c(FALSE,rep(TRUE,length(now$mask)-1L))))
  metrics <- p26_season_metrics(now,inputs,c26,"M2")
  for(y in c(2015L,2017L)) {
    a <- metrics[metrics$season==y,];b <- ev$metrics[ev$metrics$season==y,]
    put(paste0("peak_",y),abs(a$fitted_peak-b$fitted_peak))
    put(paste0("cog_",y),abs(a$fitted_cog-b$fitted_cog))
  }
  clock <- p26_clock(inputs)
  put("prepared_week_53",as.numeric(any(clock$year==2015L & clock$week==53L) ||
    sum(clock$year==2020L & clock$week==53L)!=1L ||
    length(p27_window(clock,2020L,season=TRUE))!=53L))
  checks$tolerance <- c26$replay_tolerance
  checks$passed <- is.finite(checks$error) & checks$error<=checks$tolerance
  if(!all(checks$passed)) stop("M2 exact-endpoint replay failed: ",
    paste(checks$check[!checks$passed],collapse=", "))
  list(fit=fit,evaluation=ev,score=now,checks=checks,metrics=metrics,
       transition=p26_transition(metrics,c26))
}
p27_replay <- function(ref,inputs,c26,c27,out) {
  p27_write(ref$checks,out,"replay_checks")
  p27_write(ref$score$weekly,out,"weekly_replay")
  p27_write(ref$metrics,out,"season_metrics")
  p27_write(ref$transition,out,"transition_metrics")
  p27_write(data.frame(parameter=head(ref$fit$parameter_table$parameter,length(ref$fit$endpoint)),
    value=ref$fit$endpoint),out,"parameter_vector")
  list(status="REPLAY_PASSED",reference_status=c27$reference_status,
    endpoint_hash=unname(tools::md5sum(c27$fit)),checks=ref$checks,
    transition=ref$transition)
}
p27_window <- function(clock,y,season=FALSE) {
  if(season) {
    ix <- which((clock$year==y & clock$week>=27L) |
      (clock$year==y+1L & clock$week<=26L))
    week_date <- function(year,week) {
      jan4 <- as.Date(sprintf("%04d-01-04",year))
      jan4-(as.integer(format(jan4,"%u"))-1L)+(week-1L)*7L
    }
    dates <- as.Date(vapply(seq_along(ix),function(i)
      as.character(week_date(clock$year[ix[i]],clock$week[ix[i]])),""))
    gap <- as.integer(diff(dates))
    skipped_week53 <- gap==14L & head(clock$week[ix],-1L)==52L &
      tail(clock$week[ix],-1L)==1L
    if(!length(ix) || !identical(ix,seq.int(min(ix),max(ix))) ||
       clock$week[ix[1L]]!=27L || clock$week[tail(ix,1L)]!=26L ||
       any(gap!=7L & !skipped_week53) ||
       length(ix)!=52L+sum(clock$week[ix]==53L))
      stop("Incomplete season window: ",y)
    return(ix)
  }
  if(y<=2019L) return(p24_window(clock,y))
  ix <- which(clock$year==y)
  if(length(ix)<52L || !identical(ix,seq.int(min(ix),max(ix))))
    stop("Incomplete future calendar window: ",y)
  ix
}
p27_profile <- function(score,inputs,u=0,v=0,origin="2013",year=NULL,z=0) {
  if(origin=="2013") {
    if(z!=0) stop("Initial-I intervention requires a matched-season origin")
    return(p26_variant(score,inputs,u,v)$expected)
  }
  p27_state_profile(score,inputs,year,u,z,v)$expected
}
p27_restart_state <- function(saved,u=0,z=0,donor_state=NULL) {
  if(!identical(names(saved),c("S","I","R")) || any(!is.finite(saved)))
    stop("Saved state must contain finite S, I and R")
  mass <- sum(saved);state <- saved
  if(is.null(donor_state)) {
    state["S"] <- saved["S"]*exp(u)
    state["I"] <- saved["I"]*exp(z)
  } else {
    if(!identical(names(donor_state),c("S","I"))) stop("Donor state must name S and I")
    state[c("S","I")] <- donor_state
  }
  state["R"] <- mass-state["S"]-state["I"]
  if(any(!is.finite(state)) || any(state<=0)) stop("Infeasible matched state")
  state
}
p27_state_profile <- function(score,inputs,year,u=0,z=0,v=0,
                              donor_state=NULL) {
  if(is.null(year)) stop("Matched intervention requires a year")
  clock <- p26_clock(inputs);ix <- p27_window(clock,year,season=TRUE);at <- ix[1L]-1L
  if(at<1L) stop("No saved pre-window state")
  p <- score$parameters;weekly <- score$weekly;pop <- inputs$replay$data$fixed_parameters$pop
  fixed <- list(pop=pop,mu=ref_or(p$mu,inputs$replay$data$fixed_parameters$mu),
    gamma=ref_or(p$gamma,inputs$replay$data$fixed_parameters$gamma),
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  state <- unlist(score$states[at,c("S","I","R")],use.names=TRUE)
  mass <- sum(state);state <- p27_restart_state(state,u,z,donor_state)
  end <- nrow(weekly)
  effective <- weekly$effective_beta[at:end]
  effective[seq.int(2L,length(ix)+1L)] <-
    effective[seq.int(2L,length(ix)+1L)]*exp(v)
  beta <- vapply(1:52,function(w) weekly$seasonal_beta[which(weekly$seasonal_week==w)[1L]],numeric(1))
  sim <- simulate_weekly_sirs(beta,state,fixed,weekly$seasonal_week[at:end],
    rep(1,length(effective)),initial_mass=mass,beta_by_transition=effective)
  list(expected=c(weekly$expected[seq_len(at)],sim$expected_cases[-1L]),
    states=sim,restart=at,initial_state=state,saved_mass=mass,
    effective_beta=effective,window=ix)
}
p27_peak_record <- function(values,clock,year,id,padding=13L,season=FALSE) {
  ix <- p27_window(clock,year,season);start <- max(1L,ix[1L]-1L)
  end <- min(length(values),tail(ix,1L)+padding)
  curve <- values[start:end]
  z <- list(states=data.frame(expected_cases=curve))
  p <- p24_peaks(z,id,"M2",year,"BETA_ONLY",0,0,0L)
  if(nrow(p)) {
    p$position <- start+p$peak_id-1L-ix[1L]+1L
    p$in_window <- p$position>=1L & p$position<=length(ix)
    p$support <- p$support & (start+p$peak_id-1L)>1L &
      (start+p$peak_id-1L)<length(values)
    p$genuine <- p$genuine & p$support
  }
  w <- values[ix];top <- which(w==max(w))
  grid <- data.frame(scenario=id,target=as.character(year),reference="M2",path="BETA_ONLY",
    u=0,v=0,selected_week=top[1L],selected_height=w[top[1L]],
    window_length=length(ix),
    peak_count=if(nrow(p)) sum(p$genuine & p$in_window) else 0L,
    selected_branch=NA_character_,match_status="UNTRACKED",ties=length(top),
    right_support=end>tail(ix,1L),stringsAsFactors=FALSE)
  list(grid=grid,peaks=p,curve=data.frame(scenario=id,season=year,
    position=seq.int(start,end)-ix[1L]+1L,expected=curve,
    in_window=seq.int(start,end) %in% ix))
}
p27_height_gaps <- function(g,p,c25) {
  p27_bind(lapply(seq_len(nrow(g)),function(k) {
    a <- g[k,];q <- p[p$scenario==a$scenario & p$genuine & p$in_window,]
    early <- q[q$position<=c25$early_max_week,,drop=FALSE]
    late <- q[q$position>=c25$late_min_week,,drop=FALSE]
    if(nrow(early)!=1L || nrow(late)!=1L) return(NULL)
    data.frame(season=as.integer(a$target),scenario=a$scenario,path=a$path,
      u=a$u,v=a$v,early_branch=early$branch,
      late_branch=late$branch,early_week=early$position,
      late_week=late$position,early_height=early$height,
      late_height=late$height,signed_gap=early$height-late$height,
      pair_count=nrow(q),trackable=all(c(early$match_status,late$match_status)!="AMBIGUOUS"))
  }))
}
p27_classify <- function(g,p,gaps) {
  if(nrow(g)<2L) return(data.frame())
  rows <- vector("list",nrow(g)-1L)
  for(k in seq_len(nrow(g)-1L)) {
    a <- g[k,];b <- g[k+1L,]
    ga <- gaps[gaps$scenario==a$scenario,,drop=FALSE]
    gb <- gaps[gaps$scenario==b$scenario,,drop=FALSE]
    boundary <- a$selected_week %in% c(1L,a$window_length) ||
      b$selected_week %in% c(1L,b$window_length) ||
      !a$right_support || !b$right_support
    ambiguous <- a$ties!=1L || b$ties!=1L ||
      a$match_status=="WINDOW_OR_AMBIGUOUS" ||
      b$match_status=="WINDOW_OR_AMBIGUOUS" ||
      any(p$match_status[p$scenario %in% c(a$scenario,b$scenario)]=="AMBIGUOUS")
    exchange <- nrow(ga)==1L && nrow(gb)==1L &&
      isTRUE(ga$trackable) && isTRUE(gb$trackable) &&
      identical(ga$early_branch,gb$early_branch) &&
      identical(ga$late_branch,gb$late_branch) &&
      ga$signed_gap*gb$signed_gap<0 &&
      !is.na(a$selected_branch) && !is.na(b$selected_branch) &&
      !identical(a$selected_branch,b$selected_branch)
    cls <- if(boundary) "WINDOW_SELECTION" else if(ambiguous) "UNRESOLVED" else
      if(exchange) "DOMINANCE_EXCHANGE" else if(a$peak_count!=b$peak_count)
        "SAMPLED_PEAK_LOSS" else if(!is.na(a$selected_branch) &&
          identical(a$selected_branch,b$selected_branch) && a$selected_week!=b$selected_week)
        "CONTINUOUS_DISPLACEMENT" else if(!is.na(a$selected_branch) &&
          identical(a$selected_branch,b$selected_branch))
        "STABLE_SAME_BRANCH" else "UNRESOLVED"
    rows[[k]] <- data.frame(season=as.integer(a$target),path=a$path,u=a$u,
      v_left=a$v,v_right=b$v,width=b$v-a$v,
      left=a$scenario,right=b$scenario,left_week=a$selected_week,right_week=b$selected_week,
      left_branch=a$selected_branch,right_branch=b$selected_branch,
      left_count=a$peak_count,right_count=b$peak_count,
      early_minus_late_left=if(nrow(ga)) ga$signed_gap else NA_real_,
      early_minus_late_right=if(nrow(gb)) gb$signed_gap else NA_real_,
      classification=cls,stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}
p27_mechanism <- function(ref,inputs,c27,c25,out) {
  clock <- p26_clock(inputs);base <- p22_mechanism_grid(list(mechanism_grid=c27$mechanism_grid))
  targeted <- p25_verify(c25,phase22_config(),phase24_config())$targets
  beta_target <- targeted[targeted$path=="BETA_ONLY",,drop=FALSE]
  joint <- targeted[targeted$path=="JOINT",,drop=FALSE]
  if(nrow(joint)!=1L || nrow(beta_target)!=1L)
    stop("Phase 25 targeted path identity changed")
  paths <- c(lapply(2013:2023,function(y)
    list(year=y,path="BETA_ONLY",u=0,origin="2013",grid=base)),
    list(list(year=2017L,path="MATCHED_BETA_ONLY",u=0,origin="matched",
        grid=sort(unique(c(base,-base,beta_target$x_left,beta_target$x_right)))),
      list(year=2017L,path="MATCHED_JOINT",u=joint$u,origin="matched",
        grid=sort(unique(c(base,-base,joint$x_left,joint$x_right))))))
  all <- list();n <- 0L
  for(policy in c27$boundary_policies) {
    score <- if(policy=="raw") ref$score else p26_score(ref$fit$endpoint,inputs,
      p22_problem(inputs),phase26_config(),"M2",ref$fit$penalty,blend=TRUE)
    for(spec in paths) {
      y <- spec$year;vals <- spec$grid;cache <- list();round <- 0L
      repeat {
        missing <- vals[!vapply(vals,function(v) any(vapply(cache,function(r)
          abs(r$grid$v-v)<1e-14,logical(1))),logical(1))]
        for(v in missing) {
          n <- n+1L;id <- sprintf("%s-%s-%d-%06d",policy,spec$path,y,n)
          values <- p27_profile(score,inputs,u=spec$u,v=v,
            origin=spec$origin,year=y)
          r <- p27_peak_record(values,clock,y,id,c27$padding_weeks,
            season=spec$origin=="matched")
          r$grid$v <- v;r$grid$u <- spec$u;r$grid$path <- spec$path
          if(nrow(r$peaks)) {
            r$peaks$v <- v;r$peaks$u <- spec$u;r$peaks$path <- spec$path
          }
          r$curve$v <- v;r$curve$u <- spec$u;r$curve$path <- spec$path
          cache[[id]] <- r
        }
        tracked <- p25_track(cache,c25)
        gaps <- p27_height_gaps(tracked$grid,tracked$peaks,c25)
        edges <- p27_classify(tracked$grid,tracked$peaks,gaps)
        trigger <- edges$classification %in% c("DOMINANCE_EXCHANGE","UNRESOLVED",
          "SAMPLED_PEAK_LOSS","WINDOW_SELECTION") |
          edges$left_week!=edges$right_week | edges$left_count!=edges$right_count |
          (is.finite(edges$early_minus_late_left) &
           is.finite(edges$early_minus_late_right) &
           edges$early_minus_late_left*edges$early_minus_late_right<=0)
        candidate <- (edges$v_left+edges$v_right)/2
        candidate <- candidate[trigger & (edges$v_right-edges$v_left)>c27$refinement_width]
        candidate <- setdiff(candidate,vals)
        round <- round+1L
        if(!length(candidate) || round>c27$refinement_rounds ||
           length(vals)>=c27$refinement_cap) break
        vals <- sort(c(vals,head(candidate,c27$refinement_cap-length(vals))))
      }
      all[[paste(policy,spec$path,y)]] <- list(grid=tracked$grid,peaks=tracked$peaks,
        curves=p27_bind(lapply(cache,`[[`,"curve")),edges=edges,
        gaps=gaps,policy=policy)
    }
  }
  table <- function(key) p27_bind(lapply(all,function(x)
    if(nrow(x[[key]])) transform(x[[key]],policy=x$policy) else data.frame()))
  g <- table("grid");p <- table("peaks");e <- table("edges")
  curves <- table("curves");gaps <- table("gaps")
  p27_write(g,out,"peak_path_grid");p27_write(p,out,"peak_branches")
  p27_write(e,out,"mechanism_intervals");p27_write(curves,out,"peak_curves")
  p27_write(gaps,out,"paired_height_gaps")
  counts <- aggregate(rep(1,nrow(e)),e[c("policy","season","path","classification")],sum)
  names(counts)[5L] <- "intervals";p27_write(counts,out,"interval_counts")
  plotted <- utils::read.csv(file.path(out,"peak_path_grid.csv"))
  grDevices::pdf(file.path(out,"peak_paths.pdf"),width=9,height=6)
  for(y in c(2015L,2017L,2021L)) {
    q <- plotted[plotted$policy=="raw" & plotted$path=="BETA_ONLY" & plotted$target==y,]
    plot(q$v,q$selected_week,type="b",main=paste("M2 season",y),xlab="Log beta multiplier",ylab="Selected week")
  }
  grDevices::dev.off()
  plotted <- if(nrow(gaps)) utils::read.csv(file.path(out,"paired_height_gaps.csv")) else
    data.frame(policy=character(),path=character(),season=integer(),
      v=numeric(),signed_gap=numeric())
  grDevices::pdf(file.path(out,"competing_peak_gaps.pdf"),width=8,height=6)
  for(y in c(2015L,2017L)) {
    q <- plotted[plotted$policy=="raw" & plotted$path=="BETA_ONLY" & plotted$season==y,]
    if(nrow(q)) {
      plot(q$v,q$signed_gap,type="b",xlab="Log beta multiplier",
        ylab="Early minus late height",main=paste("M2 season",y))
      abline(h=0,lty=2)
    } else plot.new()
  }
  grDevices::dev.off()
  list(status="MECHANISM_SAVED",grid=g,peaks=p,intervals=e,
    gaps=gaps,counts=counts,targeted_paths=targeted)
}
p27_peak_for_sensitivity <- function(values,clock,y,id,season=FALSE) {
  r <- p27_peak_record(values,clock,y,id,season=season)
  p <- r$peaks;g <- r$grid
  if(!nrow(p) || g$ties!=1L || g$selected_week %in% c(1L,g$window_length))
    return(list(status="WINDOW_OR_TIE",grid=g,peaks=p,time=NA_real_))
  win <- p[p$in_window & p$genuine & p$position==g$selected_week,]
  if(nrow(win)!=1L) return(list(status="NO_ELIGIBLE_PEAK",grid=g,peaks=p,time=NA_real_))
  at <- which(r$curve$position==win$position)
  q <- p10_interpolate(r$curve$expected,at,at,p10_defaults()$identity)
  list(status=if(q$status=="OK") "ELIGIBLE" else q$status,
    grid=g,peaks=p,time=if(q$status=="OK") win$position+q$position-at else NA_real_,
    peak=win$peak_id)
}
p27_pair_eligibility <- function(baseline_status,side_status) {
  ok <- identical(baseline_status,"ELIGIBLE") &&
    length(side_status)==4L && all(side_status=="ELIGIBLE")
  list(eligible=ok,reason=if(ok) "" else
    paste(unique(c(if(length(side_status)!=4L) "INCOMPLETE_PAIR",
      if(baseline_status!="ELIGIBLE") baseline_status,
      side_status[side_status!="ELIGIBLE"])),collapse=";"))
}
p27_joint_score <- function(inputs,source_cfg) {
  r <- p22_references(inputs,source_cfg)$JOINT
  z <- p22_score(r$coords,r$post,inputs,r$setup)
  clock <- p26_clock(inputs);n <- nrow(inputs$replay$data$time_markers)
  seasonal <- z$beta[clock$seasonal_week]
  z$weekly <- data.frame(expected=z$states$expected_cases,
    seasonal_week=clock$seasonal_week,seasonal_beta=seasonal,
    effective_beta=seasonal*c(rep(1,n),z$q),in_likelihood=z$mask,
    observed=z$observed)
  z
}
p27_sensitivity <- function(ref,inputs,c27,source_cfg,out) {
  clock <- p26_clock(inputs);rows <- list();variants <- list()
  joint <- p27_joint_score(inputs,source_cfg)
  for(key in c("M2_raw","M2_blend","JOINT_raw")) {
    reference <- sub("_.*$","",key);policy <- sub("^.*_","",key)
    score <- if(reference=="JOINT") joint else if(policy=="raw") ref$score else
      p26_score(ref$fit$endpoint,inputs,p22_problem(inputs),phase26_config(),
        "M2",ref$fit$penalty,blend=TRUE)
    for(design in c("2013_ORIGIN","MATCHED_SEASON_ORIGIN")) {
      years <- if(design=="2013_ORIGIN") 2013:2023 else 2013:2022
      for(y in years) {
        origin <- if(design=="2013_ORIGIN") "2013" else "matched"
        baseline <- p27_peak_for_sensitivity(score$weekly$expected,clock,y,
          paste(policy,design,y,"base",sep="/"),season=origin=="matched")
        for(h in c27$sensitivity_steps) {
          sides <- list()
          for(axis in c("S","beta")) for(sign in c(-1,1)) {
            key <- paste(axis,sign,sep="/")
            r <- tryCatch({
              x <- p27_profile(score,inputs,if(axis=="S") sign*h else 0,
                if(axis=="beta") sign*h else 0,origin,y)
              p27_peak_for_sensitivity(x,clock,y,paste(policy,design,y,h,key,sep="/"),
                season=origin=="matched")
            },error=function(e) list(status=paste0("INFEASIBLE: ",conditionMessage(e)),
              time=NA_real_,grid=data.frame(selected_week=NA_integer_),peaks=data.frame()))
            if(r$status=="ELIGIBLE" && baseline$status=="ELIGIBLE") {
              f <- p10_assignment(baseline$peaks,r$peaks,
                p10_defaults()$identity)
              rev <- p10_assignment(baseline$peaks,r$peaks,
                p10_defaults()$identity,TRUE)
              m <- f[f$status=="MATCHED" & f$left_peak==baseline$peak,]
              if(nrow(m)!=1L || m$right_peak!=r$peak ||
                 !p10_identity_consistent(f,rev)) r$status <- "BRANCH_MISMATCH"
            }
            sides[[key]] <- r
            variants[[length(variants)+1L]] <- data.frame(reference=reference,
              policy=policy,design=design,
              year=y,h=h,axis=axis,sign=sign,status=r$status,
              selected_week=r$grid$selected_week[1L],quadratic_time=r$time)
          }
          status <- vapply(sides,`[[`,character(1),"status")
          eligibility <- p27_pair_eligibility(baseline$status,status)
          paired <- eligibility$eligible
          slope <- function(axis,field) {
            a <- sides[[paste0(axis,"/1")]][[field]]
            b <- sides[[paste0(axis,"/-1")]][[field]]
            if(paired && is.finite(a) && is.finite(b)) (a-b)/(2*h) else NA_real_
          }
          es <- slope("S","time");eb <- slope("beta","time")
          week <- function(axis) if(paired)
            (sides[[paste0(axis,"/1")]]$grid$selected_week-
             sides[[paste0(axis,"/-1")]]$grid$selected_week)/(2*h) else NA_real_
          cog <- function(axis) {
            if(!paired) return(NA_real_)
            calc <- function(sign) {
              x <- p27_profile(score,inputs,if(axis=="S") sign*h else 0,
                if(axis=="beta") sign*h else 0,origin,y)[p27_window(clock,y,origin=="matched")]
              sum(seq_along(x)*x)/sum(x)
            }
            (calc(1)-calc(-1))/(2*h)
          }
          rows[[length(rows)+1L]] <- data.frame(reference=reference,
            policy=policy,design=design,year=y,h=h,
            baseline_status=baseline$status,paired_eligible=paired,
            exclusion=eligibility$reason,
            branch_id=if(paired) paste0("baseline_peak_",baseline$peak) else NA_character_,
            E_S=es,E_beta=eb,sign_S=sign(es),sign_beta=sign(eb),
            absolute_ratio=if(is.finite(es) && abs(es)>1e-8) abs(eb/es) else NA_real_,
            discrete_S=week("S"),discrete_beta=week("beta"),
            cog_S=cog("S"),cog_beta=cog("beta"))
        }
      }
    }
  }
  x <- p27_bind(rows);v <- p27_bind(variants)
  p27_write(x,out,"paired_sensitivity");p27_write(v,out,"sensitivity_variants")
  j <- x[x$reference=="JOINT" & x$policy=="raw",
    c("design","year","h","E_S","E_beta","paired_eligible")]
  names(j)[4:6] <- c("JOINT_E_S","JOINT_E_beta","JOINT_eligible")
  compare <- merge(x[x$reference=="M2" & x$policy=="raw",],j,
    by=c("design","year","h"),all.x=TRUE)
  p27_write(compare,out,"joint_comparison")
  plotted <- utils::read.csv(file.path(out,"paired_sensitivity.csv"))
  plotted <- plotted[plotted$reference=="M2" & plotted$policy=="raw" &
    plotted$paired_eligible,]
  grDevices::pdf(file.path(out,"paired_sensitivity.pdf"),width=8,height=5)
  if(nrow(plotted)) {
    plot(plotted$E_S,plotted$E_beta,xlab="E_S (weeks per log change)",
      ylab="E_beta (weeks per log change)",main="Paired eligible M2 derivatives")
    abline(a=0,b=1,lty=2)
  } else plot.new()
  grDevices::dev.off()
  list(status="SENSITIVITY_SAVED",estimates=x,variants=v,comparison=compare)
}
p27_annual <- function(ref,inputs,source_cfg,c27,out) {
  profiles <- utils::read.csv(source_cfg$annual_profiles)
  profiles <- profiles[profiles$profile_id %in% c27$annual_profiles,]
  scenarios <- utils::read.csv(source_cfg$annual_scenarios)
  scenarios <- scenarios[scenarios$profile_id %in% profiles$profile_id &
    scenarios$axis %in% c("baseline","beta") & scenarios$u==0 &
    scenarios$status=="VALID",]
  p <- ref$score$parameters;N <- inputs$replay$data$fixed_parameters$pop
  fixed <- list(pop=N,mu=ref_or(p$mu,inputs$replay$data$fixed_parameters$mu),
    gamma=ref_or(p$gamma,inputs$replay$data$fixed_parameters$gamma),
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  rows <- list();wins <- list();peaks <- list();policy <- p10_defaults()
  for(i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i,];pr <- profiles[profiles$profile_id==sc$profile_id,]
    beta <- sc$a*pr$b*exp(pr$amplitude*cos(2*pi*((1:52)-pr$phase_week)/52))
    maxima <- sum(beta>c(tail(beta,1L),head(beta,-1L)) &
                  beta>c(tail(beta,-1L),beta[1L]))
    if(maxima!=1L) stop("Annual effective beta does not have one sampled maximum")
    n <- 1L+52L*phase16_config()$horizon_years
    week <- ((seq_len(n)-2L)%%52L)+1L
    z <- simulate_weekly_sirs(beta,make_initial_state(p$S0,p$I0,N),fixed,week,rep(1,n))
    z$year <- c(0L,rep(seq_len(phase16_config()$horizon_years),each=52L));z$week <- week
    z$initialization_only <- seq_len(n)==1L
    record <- list(id=sc$scenario_id,a=sc$a,trajectory=z)
    pk <- p10_detect(record,4L,policy$peaks,policy$identity)
    win <- p10_window(record,4L,pk,policy$peaks)
    peaks[[sc$scenario_id]] <- pk;wins[[sc$scenario_id]] <- win
    rows[[i]] <- data.frame(scenario=sc$scenario_id,profile_id=sc$profile_id,
      a=sc$a,week=win$week,status=win$status,effective_beta_maxima=maxima)
    p27_write(data.frame(scenario=sc$scenario_id,week=seq_len(n),
      incidence=z$incidence,effective_beta=beta[week],q=1),
      out,paste0("annual_trajectory_",sc$scenario_id))
  }
  x <- p27_bind(rows);pairs <- list()
  for(pr in unique(x$profile_id)) {
    q <- x[x$profile_id==pr,];q <- q[order(q$a),]
    for(k in seq_len(nrow(q)-1L)) {
      a <- peaks[[q$scenario[k]]];b <- peaks[[q$scenario[k+1L]]]
      f <- p10_assignment(a,b,policy$identity)
      rev <- p10_assignment(a,b,policy$identity,TRUE)
      m <- f[f$status=="MATCHED",]
      wa <- wins[[q$scenario[k]]];wb <- wins[[q$scenario[k+1L]]]
      same <- nrow(m)>0L && p10_identity_consistent(f,rev) &&
        any(m$left_peak==wa$peak_id & m$right_peak==wb$peak_id)
      pairs[[length(pairs)+1L]] <- data.frame(profile_id=pr,
        left=q$scenario[k],right=q$scenario[k+1L],
        week_change=wb$week-wa$week,same_branch=same,
        classification=if(wa$status=="WINDOW_SELECTION" || wb$status=="WINDOW_SELECTION")
          "WINDOW_SELECTION" else if(same) "FINITE_WITHIN_BRANCH_DISPLACEMENT" else "UNRESOLVED")
    }
  }
  pair <- p27_bind(pairs);p27_write(x,out,"annual_scenarios")
  p27_write(pair,out,"annual_pairs")
  targets <- p16t_spec()
  target_rows <- do.call(rbind,lapply(seq_len(nrow(targets)),function(i) {
    t <- targets[i,];l <- x[x$scenario==t$left_scenario,]
    r <- x[x$scenario==t$right_scenario,]
    if(nrow(l)!=1L || nrow(r)!=1L) stop("Saved annual target scenario missing")
    chain <- x[x$profile_id==t$profile_id &
      x$a>=min(l$a,r$a)-1e-12 & x$a<=max(l$a,r$a)+1e-12,]
    chain <- chain[order(chain$a),]
    links <- paste(head(chain$scenario,-1L),tail(chain$scenario,-1L),sep="/")
    edge <- pair[match(links,paste(pair$left,pair$right,sep="/")),]
    selected <- length(links)>0L && !anyNA(edge$same_branch) &&
      all(edge$same_branch) && all(edge$classification=="FINITE_WITHIN_BRANCH_DISPLACEMENT")
    data.frame(target_id=t$target_id,profile_id=t$profile_id,
      left=t$left_scenario,right=t$right_scenario,
      left_week=l$week,right_week=r$week,week_change=r$week-l$week,
      saved_absolute_weeks=c(primary=9L,secondary_low=8L,
        secondary_high=10L)[t$target_id],
      same_branch=selected,
      transferred=selected &&
        l$status=="GENUINE_LOCAL_PEAK" && r$status=="GENUINE_LOCAL_PEAK" &&
        abs(r$week-l$week)>=c(primary=9L,secondary_low=8L,
          secondary_high=10L)[t$target_id])
  }))
  p27_write(target_rows,out,"annual_targets")
  transfer <- p27_bind(lapply(split(x,x$profile_id),function(z) {
    z <- z[order(z$a),]
    q <- pair[pair$profile_id==z$profile_id[1L],]
    data.frame(profile_id=z$profile_id[1L],left_week=z$week[1L],
      right_week=tail(z$week,1L),finite_displacement=tail(z$week,1L)-z$week[1L],
      saved_target_weeks=if(z$profile_id[1L]=="profile-004") 9L else 18L,
      transfers=abs(tail(z$week,1L)-z$week[1L])>=
        (if(z$profile_id[1L]=="profile-004") 9L else 18L) &&
        all(q$same_branch) && all(q$classification=="FINITE_WITHIN_BRANCH_DISPLACEMENT"),
      all_pairs_same_branch=all(q$same_branch),
      all_pairs_interior=all(q$classification!="WINDOW_SELECTION"))
  }))
  p27_write(transfer,out,"annual_transfer")
  plotted <- utils::read.csv(file.path(out,"annual_scenarios.csv"))
  grDevices::pdf(file.path(out,"annual_paths.pdf"),width=8,height=5)
  plot(range(plotted$a),range(plotted$week),type="n",xlab="Transmission multiplier",
    ylab="Selected week",main="Strict annual-only effective transmission")
  for(pr in unique(plotted$profile_id)) {
    q <- plotted[plotted$profile_id==pr,]
    q <- q[order(q$a),];lines(q$a,q$week,type="b")
  }
  grDevices::dev.off()
  list(status="ANNUAL_SAVED",scenarios=x,pairs=pair,transfer=transfer,
    targets=target_rows)
}
p27_fourier_start <- function(ref,inputs,c26,family) {
  setup <- p22_problem(inputs,family)
  beta <- vapply(1:52,function(w) ref$score$weekly$seasonal_beta[
    which(ref$score$weekly$seasonal_week==w)[1L]],numeric(1))
  harmonics <- setup$problem$basis$harmonics
  X <- do.call(cbind,lapply(seq_len(harmonics),function(j)
    cbind(cos(2*pi*j*(0:51)/52),sin(2*pi*j*(0:51)/52))))
  params <- ref$score$parameters
  params$beta0 <- mean(beta)
  params$free_coordinates <- as.vector(qr.solve(X,log(beta)))
  parts <- p26_split(ref$fit$endpoint,inputs,p22_problem(inputs),c26,"M2")
  c(encode_parameters(params,setup$problem$map),parts$post,
    unname(head(parts$level,-1L)),as.vector(parts$shape))
}
p27_embed_a <- function(theta,inputs) {
  a <- p22_problem(inputs,"A")$problem$map$order
  b <- p22_problem(inputs,"AS")$problem$map$order
  head <- setNames(rep(0,length(b)),b)
  head[a] <- theta[seq_along(a)]
  c(unname(head),theta[-seq_along(a)])
}
p27_independent_theta <- function(inputs,c26,family) {
  setup <- p22_problem(inputs,family)
  design <- setup$cfg$initial
  design$omega <- inputs$ref$parameters$omega
  design$delta <- inputs$ref$parameters$delta
  design$free_coordinates <- rep(0,setup$problem$basis$K-1L)
  coords <- encode_parameters(design,setup$problem$map)
  coords <- pmax(setup$problem$map$lower+1e-5,
    pmin(setup$problem$map$upper-1e-5,coords))
  c(coords,rep(0,inputs$basis$dimension),
    rep(0,length(c26$historical_seasons)-1L),
    rep(0,2L*length(c26$historical_seasons)))
}
p27_fourier_score <- function(fit,inputs,c26,family,blend) {
  p26_score(fit$endpoint,inputs,p22_problem(inputs,family),c26,"M2",
    fit$penalty,blend=blend)
}
p27_maxima <- function(x) {
  n <- length(x)
  sum(x>c(x[n],head(x,-1L)) & x>c(tail(x,-1L),x[1L]))
}
p27_fit_fourier <- function(ref,inputs,c26,c27,out) {
  rows <- list();maxima <- list();fits <- list()
  progress <- getOption("phase27.progress",function(...) invisible(NULL))
  for(policy in c27$boundary_policies) {
    blend <- policy=="blend"
    for(family in c("A","AS")) {
      dir <- file.path(out,paste(policy,family,sep="_"));dir.create(dir)
      progress("FIT_START",paste(policy,family))
      start <- if(family=="A") p27_fourier_start(ref,inputs,c26,"A") else
        p27_embed_a(fits[[paste(policy,"A",sep="_")]]$endpoint,inputs)
      fit <- p26_fit(inputs,NULL,c26,"M2",ref$fit$penalty,dir,
        setup_override=p22_problem(inputs,family),baseline_override=start,
        blend=blend)
      score <- p27_fourier_score(fit,inputs,c26,family,blend)
      replay_error <- abs(score$objective-fit$score$objective)
      descent_min <- if(all(is.na(fit$numerical$minimum_delta))) NA_real_ else
        min(fit$numerical$minimum_delta,na.rm=TRUE)
      gate <- data.frame(policy=policy,family=family,replay_error=replay_error,
        replay_tolerance=c26$replay_tolerance,
        minimum_delta=descent_min,descent_tolerance=c26$descent_tolerance,
        replay_passed=replay_error<=c26$replay_tolerance,
        descent_passed=is.finite(descent_min) &&
          descent_min>=-c26$descent_tolerance)
      p27_write(gate,dir,"endpoint_gate")
      if(!gate$replay_passed || !gate$descent_passed)
        stop(sprintf("Fourier endpoint gate failed: %s %s; replay error %.8g (limit %.8g), minimum local delta %.8g (limit >= -%.8g); see %s",
          policy,family,replay_error,c26$replay_tolerance,descent_min,
          c26$descent_tolerance,file.path(dir,"endpoint_gate.csv")))
      if(family=="AS") {
        a <- fits[[paste(policy,"A",sep="_")]]
        embedded <- p26_score(p27_embed_a(a$endpoint,inputs),inputs,
          p22_problem(inputs,"AS"),c26,"M2",ref$fit$penalty,blend=blend)
        if(abs(embedded$objective-a$score$objective)>c26$replay_tolerance)
          stop("Nested A-to-AS embedding failed")
      }
      key <- paste(policy,family,sep="_");fits[[key]] <- fit
      saveRDS(fit,file.path(dir,"fit.rds"))
      progress("FIT_COMPLETED",paste(policy,family,
        "objective",format(score$objective,digits=10)))
      rows[[key]] <- data.frame(policy=policy,family=family,
        historical_loglik=score$historical_loglik,future_loglik=score$future_loglik,
        full_loglik=score$loglik,objective=score$objective,
        seasonal_penalty=score$seasonal_penalty,future_penalty=score$future_penalty,
        level_penalty=score$level_penalty,shape_penalty=score$shape_penalty,
        parameters=length(fit$endpoint),
        shared_maxima=p22_fourier_maxima(score$parameters,p22_problem(inputs,family)$problem$basis))
      clock <- p26_clock(inputs)
      for(y in 2013:2023) {
        ix <- if(y<=2019L) p27_window(clock,y) else which(clock$year==y)
        maxima[[length(maxima)+1L]] <- data.frame(policy=policy,family=family,year=y,
          effective_beta_maxima=p27_maxima(score$weekly$effective_beta[ix]))
      }
      p27_write(score$weekly,out,paste0("weekly_",key))
    }
  }
  x <- p27_bind(rows);m <- p27_bind(maxima)
  p27_write(x,out,"fourier_fits");p27_write(m,out,"effective_beta_maxima")
  plotted <- utils::read.csv(file.path(out,"fourier_fits.csv"))
  grDevices::pdf(file.path(out,"fourier_objectives.pdf"),width=8,height=5)
  barplot(plotted$objective,names.arg=paste(plotted$policy,plotted$family,sep="/"),
    ylab="Penalized objective",main="A and AS shared seasonal factor")
  grDevices::dev.off()
  list(status="FOURIER_SAVED",fits=x,maxima=m,fit_objects=fits)
}
p27_fold_mask <- function(n,held,N) {
  if(any(held<1L) || any(held>N-n) || anyDuplicated(held)) stop("Invalid held-out rows")
  mask <- c(FALSE,rep(TRUE,n-1L),rep(TRUE,N-n))
  mask[n+held] <- FALSE
  mask
}
p27_validation <- function(ref,inputs,c26,c27,source_cfg,out) {
  clock <- p26_clock(inputs);N <- nrow(clock)
  n <- nrow(inputs$replay$data$time_markers)
  rows <- list();weekly <- list()
  progress <- getOption("phase27.progress",function(...) invisible(NULL))
  for(policy in c27$boundary_policies) for(k in seq_along(source_cfg$validation_blocks)) {
    held <- source_cfg$validation_blocks[[k]]
    mask <- p27_fold_mask(n,held,N)
    if(any(mask[n+held])) stop("Fold mask leaks held-out rows")
    trained <- list()
    for(family in c("A","AS")) {
      dir <- file.path(out,paste(policy,k,family,sep="_"));dir.create(dir)
      progress("FOLD_FIT_START",paste(policy,"fold",k,family))
      extra <- if(family=="AS") list(p27_embed_a(trained$A$endpoint,inputs)) else list()
      fit <- p26_fit(inputs,NULL,c26,"M2",ref$fit$penalty,dir,mask=mask,
        holdout=TRUE,setup_override=p22_problem(inputs,family),
        baseline_override=p27_independent_theta(inputs,c26,family),
        extra_starts=extra,blend=policy=="blend")
      score <- p26_score(fit$endpoint,inputs,p22_problem(inputs,family),
        c26,"M2",ref$fit$penalty,mask,blend=policy=="blend")
      replay_error <- abs(score$objective-fit$score$objective)
      descent_min <- min(fit$numerical$minimum_delta,na.rm=TRUE)
      gate <- data.frame(policy=policy,fold=k,family=family,
        replay_error=replay_error,replay_tolerance=c26$replay_tolerance,
        minimum_delta=descent_min,descent_tolerance=c26$descent_tolerance,
        replay_passed=replay_error<=c26$replay_tolerance,
        descent_passed=descent_min>=-c26$descent_tolerance)
      p27_write(gate,dir,"endpoint_gate")
      if(!gate$replay_passed || !gate$descent_passed)
        stop(sprintf("Fold endpoint gate failed: %s fold %d %s; replay error %.8g, minimum local delta %.8g; see %s",
          policy,k,family,replay_error,descent_min,
          file.path(dir,"endpoint_gate.csv")))
      trained[[family]] <- fit
      idx <- n+held
      per <- score$per_week[idx]
      rows[[length(rows)+1L]] <- data.frame(policy=policy,fold=k,family=family,
        held_rows=length(idx),held_log_score=sum(per),
        held_mae=mean(abs(score$observed[idx]-score$weekly$expected[idx])),
        training_objective=score$objective)
      weekly[[length(weekly)+1L]] <- data.frame(policy=policy,fold=k,family=family,
        row=idx,observed=score$observed[idx],expected=score$weekly$expected[idx],
        held_log_score=per)
      saveRDS(fit,file.path(dir,"fit.rds"))
      progress("FOLD_FIT_COMPLETED",paste(policy,"fold",k,family,
        "objective",format(score$objective,digits=10)))
    }
  }
  x <- p27_bind(rows);w <- p27_bind(weekly)
  pooled <- p27_bind(lapply(split(x,list(x$policy,x$family),drop=TRUE),function(z)
    data.frame(policy=z$policy[1L],family=z$family[1L],folds=nrow(z),
      held_rows=sum(z$held_rows),pooled_log_score=sum(z$held_log_score),
      pooled_mae=sum(z$held_mae*z$held_rows)/sum(z$held_rows))))
  p27_write(x,out,"fold_scores");p27_write(w,out,"fold_weekly")
  p27_write(pooled,out,"pooled_scores")
  plotted <- utils::read.csv(file.path(out,"fold_scores.csv"))
  grDevices::pdf(file.path(out,"fold_scores.pdf"),width=8,height=5)
  barplot(plotted$held_log_score,
    names.arg=paste(plotted$policy,plotted$fold,plotted$family,sep="/"),
    las=2,cex.names=.55,ylab="Held-out NB2 log score",
    main="Retrospective interpolation folds")
  grDevices::dev.off()
  list(status="VALIDATION_SAVED",folds=x,weekly=w,pooled=pooled)
}
p27_report <- function(ref,c27,source_cfg,paths,out) {
  stages <- lapply(names(paths),function(k) {
    x <- read_phase2_result(paths[[k]])
    if(!identical(x$metadata$phase,27L) ||
       !identical(x$metadata$stage,k) ||
       !identical(x$metadata$reference_status,"USER_SELECTED_M2") ||
       !identical(x$metadata$endpoint_parameter_md5,
         p27_parameter_hash(ref$fit$endpoint)))
      stop("Phase 27 report input identity mismatch: ",k)
    x$result
  })
  names(stages) <- names(paths)
  if(!all(vapply(stages,function(z) grepl("_SAVED$",z$status),logical(1))))
    stop("Phase 27 empirical stages are incomplete")
  m <- stages$mechanism$intervals;s <- stages$sensitivity$estimates
  a <- stages$annual$transfer;f <- stages$fourier$fits;v <- stages$validation$pooled
  raw <- m[m$policy=="raw" & m$path=="BETA_ONLY",]
  m2 <- s[s$reference=="M2" & s$policy=="raw",]
  blend <- s[s$reference=="M2" & s$policy=="blend",]
  target <- m2[m2$year %in% 2015:2017,]
  paired <- target[target$paired_eligible,]
  beta_dominates <- nrow(paired)>0L && all(is.finite(paired$E_S)) &&
    all(is.finite(paired$E_beta)) && all(paired$E_S<0) &&
    all(paired$E_beta<0) && all(abs(paired$E_beta)>abs(paired$E_S))
  compare_s <- merge(target,blend,by=c("design","year","h"),
    suffixes=c("_raw","_blend"))
  sensitivity_boundary <- any(compare_s$paired_eligible_raw !=
    compare_s$paired_eligible_blend) ||
    any(sign(abs(compare_s$E_beta_raw)-abs(compare_s$E_S_raw)) !=
      sign(abs(compare_s$E_beta_blend)-abs(compare_s$E_S_blend)),na.rm=TRUE)
  mech_boundary <- sum(m$policy=="raw" & m$classification=="DOMINANCE_EXCHANGE") !=
    sum(m$policy=="blend" & m$classification=="DOMINANCE_EXCHANGE")
  fit_order <- function(policy) {
    q <- f[f$policy==policy,]
    if(nrow(q)!=2L || !setequal(q$family,c("A","AS"))) return(NA)
    q$full_loglik[q$family=="AS"]>q$full_loglik[q$family=="A"]
  }
  cv_order <- function(policy) {
    q <- v[v$policy==policy,]
    if(nrow(q)!=2L || !setequal(q$family,c("A","AS"))) return(NA)
    q$pooled_log_score[q$family=="AS"]>q$pooled_log_score[q$family=="A"]
  }
  fourier_boundary <- !identical(fit_order("raw"),fit_order("blend")) ||
    !identical(cv_order("raw"),cv_order("blend"))
  key <- c("peak movement mechanism","S versus beta timing sensitivity",
    "annual examples and A versus AS fit")
  verdict <- c(if(any(raw$classification=="DOMINANCE_EXCHANGE") ||
      any(m$policy=="raw" & m$path=="MATCHED_JOINT" &
        m$classification=="DOMINANCE_EXCHANGE")) "REVISED_UNDER_M2" else
      if(any(raw$classification=="UNRESOLVED") || mech_boundary)
        "INCONCLUSIVE_UNDER_M2" else "SUPPORTED_UNDER_M2",
    if(!nrow(paired) || !all(2015:2017 %in% paired$year) || sensitivity_boundary)
      "INCONCLUSIVE_UNDER_M2" else if(beta_dominates)
        "SUPPORTED_UNDER_M2" else "REVISED_UNDER_M2",
    if(anyNA(c(fit_order("raw"),cv_order("raw"))) || fourier_boundary)
      "INCONCLUSIVE_UNDER_M2" else if(isTRUE(fit_order("raw")) &&
        isTRUE(cv_order("raw")) && all(a$transfers)) "SUPPORTED_UNDER_M2" else
      "REVISED_UNDER_M2")
  table <- data.frame(claim=key,verdict=verdict,
    eligible_years=c(paste(sort(unique(raw$season)),collapse=","),
      paste(sort(unique(m2$year[m2$paired_eligible])),collapse=","),
      "synthetic year 4; observed 2013-2023"),
    numerator=c(sum(raw$classification=="DOMINANCE_EXCHANGE"),
      sum(m2$paired_eligible & is.finite(m2$E_S) & is.finite(m2$E_beta) &
        m2$E_S<0 & m2$E_beta<0 & abs(m2$E_beta)>abs(m2$E_S)),
      sum(a$transfers)),
    denominator=c(nrow(raw),sum(m2$paired_eligible),nrow(a)),
    count_definition=c("verified exchanges / raw BETA_ONLY intervals",
      "eligible M2 raw year-steps with both negative and larger beta response / eligible M2 raw year-steps",
      "transferred annual profiles / two saved profiles"),
    boundary_sensitive=c(mech_boundary,sensitivity_boundary,fourier_boundary),
    source_csv=c(file.path(dirname(paths$mechanism),"mechanism_intervals.csv"),
      file.path(dirname(paths$sensitivity),"paired_sensitivity.csv"),
      file.path(dirname(paths$annual),"annual_transfer.csv")),
    stringsAsFactors=FALSE)
  p27_write(table,out,"claim_verdicts")
  old_fourier <- read_phase2_result(c27$joint_fourier)$result$summary
  old_validation <- read_phase2_result(c27$joint_validation)$result$fold_scores
  old_mechanism <- utils::read.csv(file.path(dirname(c27$phase25),
    "transition_brackets.csv"))
  p27_write(old_fourier,out,"phase22_fourier_comparator")
  p27_write(old_validation,out,"phase22_validation_comparator")
  p27_write(old_mechanism,out,"phase25_transition_comparator")
  evidence <- data.frame(claim=key,
    original_claim=c("2018/19 grid: no verified dominance exchange",
      "beta timing response larger in eight eligible 2018/19 steps",
      "nine- and eighteen-week annual-only examples; AS historical fit advantage"),
    saved_JOINT_result=c("Phase 25 JOINT targets A and B each verified a local dominance exchange; A also retained an unresolved bracket",
      "Phase 27 replays JOINT under identical origins and windows",
      "Phase 22 Fourier uses a different correction structure"),
    M2_result=verdict,
    source_csv=c(file.path(dirname(paths$mechanism),"mechanism_intervals.csv"),
      file.path(dirname(paths$sensitivity),"paired_sensitivity.csv"),
      file.path(dirname(paths$fourier),"fourier_fits.csv")),
    figure=c(file.path(dirname(paths$mechanism),"competing_peak_gaps.pdf"),
      file.path(dirname(paths$sensitivity),"paired_sensitivity.pdf"),
      file.path(dirname(paths$annual),"annual_paths.pdf")),
    claim_type=c("fitted-model conditional result","fitted-model conditional result",
      "simulation existence result and fitted-model conditional result"),
    replacement_sentence=c(
      sprintf("Under the user-selected M2 fit, %d of %d evaluated full-history transmission-grid intervals showed a verified sampled-peak dominance exchange; interval counts are not mechanism probabilities. The fitted 2015/16-to-2017/18 movement was -11 weeks versus -13 observed.",
        table$numerator[1L],table$denominator[1L]),
      sprintf("Under the user-selected M2 fit, %d of %d eligible year-step comparisons had negative timing derivatives for both S and complete beta, with the beta response larger in magnitude; comparisons use the same branch, origin, and log step.",
        table$numerator[2L],table$denominator[2L]),
      sprintf("Under M2 nuisance coordinates, %d of %d strict annual-only synthetic profiles retained their saved finite displacement. The observed-data A versus AS comparison changes the shared seasonal factor within the M2 correction family and is not a strict annual-only effective-transmission test.",
        table$numerator[3L],table$denominator[3L])))
  p27_write(evidence,out,"manuscript_evidence")
  tex <- if(file.exists("main.tex")) "main.tex" else {
    files <- list.files(".","main\\.tex$",recursive=TRUE)
    if(length(files)) files[1L] else NA_character_
  }
  hits <- if(!is.na(tex)) {
    lines <- readLines(tex,warn=FALSE)
    grep("no .*dominance exchange|no .*exchange.*verified",lines,ignore.case=TRUE)
  } else integer()
  changes <- data.frame(source=c(if(length(hits)) rep(tex,length(hits)),
    "final_1.pdf"),location=c(if(length(hits)) paste0("line ",hits),
    "Figure 1 caption"),
    current=c(if(length(hits)) lines[hits],
      "No verified dominance exchange was recorded."),
    replacement="The original grid audit found no verified exchange; the targeted Phase 25 audit verified local JOINT exchanges. The user-selected M2 mechanism results are reported separately.")
  p27_write(changes,out,"manuscript_changes")
  writeLines(c("% Phase 27 manuscript replacement candidates",
    "% Final wording requires interpretation of the empirical Phase 27 artifacts.",
    paste0("% Source: ",changes$source,"; ",changes$location),
    "The original grid audit found no verified dominance exchange, while the targeted Phase 25 audit verified local exchanges under JOINT. The user-selected M2 results are reported separately."),
    file.path(out,"replacement_snippets.tex"))
  report <- c("# Phase 27: user-selected M2 re-audit","",
    "Reference: USER_SELECTED_M2. Phase 26 selected no adequate candidate.",
    "M2 fits a -11-week 2015/16-to-2017/18 peak shift versus -13 weeks observed;",
    "its 2015/16 endpoint peak misses by 3 weeks and the center-of-gravity change by 1.19 weeks.",
    "These empirical audits do not identify a causal explanation of the actual Japanese shift.","",
    "## Fit to observed cases","",capture.output(print(f,row.names=FALSE)),"",
    "## Counterfactual peak mechanisms within M2","",
    "Grid intervals are not independent epidemics or mechanism probabilities.",
    capture.output(print(table[1L,],row.names=FALSE)),"",
    "## Timing sensitivity within M2","",
    capture.output(print(table[2L,],row.names=FALSE)),"",
    "## Strict annual-only existence examples","",
    capture.output(print(a,row.names=FALSE)),"",
    "## Seasonal-factor Fourier comparison","",
    "A versus AS changes the shared seasonal factor conditional on M2 historical shape corrections and future q.",
    "It is not a strict annual-only effective-transmission test.",
    capture.output(print(v,row.names=FALSE)),"",
    "## Causal identification","",
    "The actual Japanese shift remains causally unidentified by these fitted-model interventions.","",
    "## Manuscript changes","",
    "See manuscript_evidence.csv and replacement_snippets.tex. main.tex was not edited.","",
    "## Reproduction","",
    "Run the commands in docs/phase27_protocol.md in order, using fresh output directories.")
  writeLines(report,file.path(out,"phase27_report.md"))
  list(status="REPORT_SAVED",verdicts=table,evidence=evidence)
}
