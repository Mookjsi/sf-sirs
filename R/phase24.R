# Phase 24: saved-model audit and matched season-start interventions.
p24_read_stage <- function(path,stage,cfg) {
  x <- read_phase2_result(path)
  if(!identical(x$metadata$stage,stage) || !identical(x$metadata$phase,24L))
    stop("Wrong Phase 24 stage: ",path)
  h <- tools::md5sum(x$metadata$input_paths)
  if(anyNA(h) || !identical(unname(h),unname(x$metadata$input_hashes)))
    stop("Stage dependency changed: ",path)
  if(!identical(x$result$config_version,cfg$version)) stop("Phase 24 configuration version changed")
  x$result
}
p24_date <- function(y,w) {
  jan4 <- as.Date(sprintf("%04d-01-04",y))
  as.character(jan4-(as.integer(format(jan4,"%u"))-1L)+(w-1L)*7L)
}
p24_clock <- function(inputs) {
  h <- inputs$replay$data$time_markers
  data.frame(row=seq_len(nrow(h)+nrow(inputs$calendar)),
    year=c(h$year,inputs$calendar$year),
    week=c(h$week,inputs$calendar$source_week),
    seasonal_week=c(inputs$replay$data$seasonal_week,inputs$calendar$seasonal_week))
}
p24_window <- function(clock,y) {
  ix <- which((clock$year==y & clock$week>=27L) |
              (clock$year==y+1L & clock$week<=26L))
  if(length(ix)!=52L || !identical(ix,seq.int(min(ix),max(ix))) ||
     clock$year[ix[1L]]!=y || clock$week[ix[1L]]!=27L ||
     clock$year[ix[52L]]!=y+1L || clock$week[ix[52L]]!=26L ||
     any(clock$week[ix]>52L)) stop("Season calendar is not 52 contiguous prepared weeks: ",y)
  ix
}
p24_refs <- function(inputs,source_cfg) {
  refs <- p22_references(inputs,source_cfg)[c("HIST","JOINT")]
  fit <- p22_saved_joint(inputs,source_cfg)$fit$primary
  rival <- setdiff(fit$competitive,fit$best)
  if(length(rival)) {
    a <- fit$attempts[[rival[1L]]]
    d <- length(refs$JOINT$setup$problem$map$order)
    if(!is.null(a$coordinates) && length(a$coordinates)==d+inputs$basis$dimension)
      refs$JOINT_COMPETING <- list(coords=a$coordinates[seq_len(d)],
        post=tail(a$coordinates,inputs$basis$dimension),setup=refs$JOINT$setup,
        label=paste("saved competitive endpoint",rival[1L]))
  }
  refs
}
p24_scores <- function(inputs,source_cfg,tol) {
  refs <- p24_refs(inputs,source_cfg)
  out <- list();audit <- list()
  for(nm in names(refs)) {
    z <- tryCatch(p22_score(refs[[nm]]$coords,refs[[nm]]$post,inputs,refs[[nm]]$setup),error=function(e)e)
    if(inherits(z,"error")) {
      if(nm=="JOINT_COMPETING") next
      stop("Reference replay failed: ",nm,": ",conditionMessage(z))
    }
    if(nm=="HIST") {
      old <- inputs$ref$evaluation$states
      err <- max(abs(as.matrix(z$states[seq_len(nrow(old)),c("S","I","R","incidence","expected_cases")])-
        as.matrix(old[,c("S","I","R","incidence","expected_cases")])))
      if(!is.finite(err) || err>tol) stop("Accepted HIST replay differs: ",err)
    } else {
      fit <- p22_saved_joint(inputs,source_cfg)$fit$primary
      id <- if(nm=="JOINT") fit$best else setdiff(fit$competitive,fit$best)[1L]
      old <- fit$attempts[[id]]$score$states
      err <- max(abs(as.matrix(z$states[,c("S","I","R","incidence","expected_cases")])-
        as.matrix(old[,c("S","I","R","incidence","expected_cases")])))
      if(!is.finite(err) || err>tol) {
        if(nm=="JOINT_COMPETING") next
        stop("JOINT saved replay differs: ",err)
      }
    }
    out[[nm]] <- z
    audit[[nm]] <- data.frame(reference=nm,full_replay_max_absolute_error=err,
      tolerance=tol,passed=TRUE,source=refs[[nm]]$label)
  }
  list(scores=out,replay=do.call(rbind,audit))
}
p24_stats <- function(v,clock,ix,reference,season,series) {
  w <- seq_along(ix);a <- v[ix];top <- which(a==max(a));pk <- top[1L]
  data.frame(reference=reference,season=season,series=series,
    peak_season_week=pk,peak_year=clock$year[ix[pk]],peak_calendar_week=clock$week[ix[pk]],
    peak_date=p24_date(clock$year[ix[pk]],clock$week[ix[pk]]),peak_height=a[pk],
    tied_season_weeks=paste(top,collapse=";"),ties=length(top),
    center_of_gravity=if(sum(a)>0) sum(w*a)/sum(a) else NA_real_,
    center_of_gravity_date=if(sum(a)>0) p24_date(clock$year[ix[round(sum(w*a)/sum(a))]],
       clock$week[ix[round(sum(w*a)/sum(a))]]) else NA_character_,
    total=sum(a),stringsAsFactors=FALSE)
}
p24_validate <- function(cfg,source_cfg) {
  inputs <- p22_inputs(source_cfg);clock <- p24_clock(inputs)
  prior <- lapply(cfg$phase23,read_phase2_result)
  if(!identical(vapply(prior,function(z) z$metadata$stage,character(1)),
                c("mechanism","sensitivity","matched","report")))
    stop("Phase 23 saved-stage identity mismatch")
  prior_summary <- prior[[1L]]$result$summary
  prior_summary <- prior_summary[prior_summary$reference %in% c("HIST","JOINT") &
    prior_summary$year %in% cfg$seasons,,drop=FALSE]
  if(any(clock$year==2015L & clock$week==53L) ||
     !identical(as.integer(inputs$replay$data$time_markers$week[which(inputs$replay$data$time_markers$year==2015L)]),1:52))
    stop("Accepted 2015 W53 merge policy changed")
  for(y in cfg$seasons) p24_window(clock,y)
  weekly <- utils::read.csv(cfg$weekly)
  if(!all(c("reference","row","year","week","observed","expected") %in% names(weekly))) stop("Weekly schema changed")
  result <- p24_scores(inputs,source_cfg,cfg$replay_tolerance)
  for(nm in names(result$scores)) {
    a <- weekly[weekly$reference==if(nm=="JOINT_COMPETING") "JOINT" else nm,]
    if(nrow(a)!=nrow(clock) || any(a$year!=clock$year | a$week!=clock$week) ||
       (nm!="JOINT_COMPETING" &&
        max(abs(a$expected-result$scores[[nm]]$states$expected_cases))>cfg$replay_tolerance))
      stop("Phase 22 weekly calendar or fit mismatch: ",nm)
  }
  list(status="INPUTS_VALIDATED",config_version=cfg$version,replay=result$replay,
    calendar=clock,source_hashes=tools::md5sum(c(source_cfg$reference,source_cfg$replay,
      source_cfg$joint_fit,cfg$weekly,cfg$phase23)),
    phase23_beta_reference=prior_summary,
    week53_policy="2015-W53 merged into 2015-W52 before fitting")
}
p24_fixed <- function(score,inputs) {
  p <- score$parameters;d <- inputs$replay$data$fixed_parameters
  list(pop=d$pop,mu=ref_or(p$mu,d$mu),gamma=ref_or(p$gamma,d$gamma),
       omega=ref_or(p$omega,inputs$ref$parameters$omega),
       delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
}
p24_replay_segment <- function(score,inputs,ix,padding,tol) {
  z <- score$states;start <- ix[1L]-1L
  if(start<1L || inputs$replay$data$time_markers$week[start]!=26L) stop("Start must follow week 26")
  end <- min(nrow(z),ix[52L]+padding)
  state <- unlist(z[start,c("S","I","R")],use.names=TRUE)
  weeks <- z$week[start:end]
  q <- c(rep(1,nrow(inputs$replay$data$time_markers)),score$q)
  fixed <- p24_fixed(score,inputs)
  r <- simulate_weekly_sirs(score$beta,state,fixed,weeks,q[start:end],initial_mass=sum(state))
  cols <- c("S","I","R","incidence","expected_cases")
  # Row one is initialization; compare transitions only.
  err <- max(abs(as.matrix(r[-1L,cols])-as.matrix(z[(start+1L):end,cols])))
  if(!is.finite(err) || err>tol) stop("Week-26 replay alignment failed: ",err)
  list(state=state,error=err,end=end)
}
p24_baseline <- function(cfg,source_cfg,validate_path) {
  valid <- p24_read_stage(validate_path,"validate-inputs",cfg)
  inputs <- p22_inputs(source_cfg);clock <- p24_clock(inputs)
  if(!identical(clock,valid$calendar)) stop("Calendar changed")
  s <- p24_scores(inputs,source_cfg,cfg$replay_tolerance)$scores
  weekly <- utils::read.csv(cfg$weekly)
  metrics <- starts <- checks <- curves <- annual <- list()
  for(nm in names(s)) for(y in cfg$seasons) {
    ix <- p24_window(clock,y);z <- s[[nm]];r <- p24_replay_segment(z,inputs,ix,cfg$padding_weeks,cfg$replay_tolerance)
    starts[[length(starts)+1L]] <- data.frame(reference=nm,season=y,start_row=ix[1L]-1L,
      start_year=y,start_week=26L,S=r$state["S"],I=r$state["I"],R=r$state["R"],
      mass=sum(r$state),row.names=NULL)
    checks[[length(checks)+1L]] <- data.frame(reference=nm,season=y,
      transition_start=ix[1L],replayed_through=r$end,max_absolute_error=r$error,passed=TRUE)
    obs <- weekly$observed[weekly$reference=="HIST"]
    metrics[[length(metrics)+1L]] <- p24_stats(obs,clock,ix,nm,y,"OBSERVED")
    metrics[[length(metrics)+1L]] <- p24_stats(z$states$expected_cases,clock,ix,nm,y,"FITTED")
    curves[[length(curves)+1L]] <- data.frame(reference=nm,season=y,season_week=1:52,
      year=clock$year[ix],calendar_week=clock$week[ix],
      date=mapply(p24_date,clock$year[ix],clock$week[ix]),
      observed=obs[ix],expected=z$states$expected_cases[ix])
  }
  for(nm in names(s)) for(y in cfg$seasons) {
    ix <- which(clock$year==y)
    obs <- weekly$observed[weekly$reference=="HIST"]
    for(kind in c("OBSERVED","FITTED")) {
      a <- if(kind=="OBSERVED") obs[ix] else s[[nm]]$states$expected_cases[ix]
      top <- which(a==max(a))
      annual[[length(annual)+1L]] <- data.frame(reference=nm,calendar_year=y,
        series=kind,peak_calendar_week=clock$week[ix[top[1L]]],
        peak_date=p24_date(y,clock$week[ix[top[1L]]]),
        ties=paste(clock$week[ix[top]],collapse=";"),peak_height=max(a))
    }
  }
  metrics <- do.call(rbind,metrics);starts <- do.call(rbind,starts)
  changes <- do.call(rbind,lapply(split(metrics,list(metrics$reference,metrics$series),drop=TRUE),function(a) {
    a <- a[order(a$season),];data.frame(reference=a$reference[1],series=a$series[1],
      peak_change=a$peak_season_week[3]-a$peak_season_week[1],
      cog_change=a$center_of_gravity[3]-a$center_of_gravity[1])
  }))
  gate <- do.call(rbind,lapply(names(s),function(nm) {
    a <- metrics[metrics$reference==nm & metrics$series=="OBSERVED",]
    b <- metrics[metrics$reference==nm & metrics$series=="FITTED",]
    a <- a[order(a$season),];b <- b[order(b$season),]
    pe <- (b$peak_season_week[3]-b$peak_season_week[1])-(a$peak_season_week[3]-a$peak_season_week[1])
    ce <- (b$center_of_gravity[3]-b$center_of_gravity[1])-(a$center_of_gravity[3]-a$center_of_gravity[1])
    se <- max(abs(b$peak_season_week-a$peak_season_week))
    data.frame(reference=nm,gate="OBSERVED_TRANSITION_REPRODUCED",peak_change_error=pe,
      cog_change_error=ce,max_season_peak_error=se,
      passed=abs(pe)<=cfg$observed_gate$max_peak_change_error_weeks &&
        abs(ce)<=cfg$observed_gate$max_cog_change_error_weeks &&
        se<=cfg$observed_gate$max_each_season_peak_error_weeks)
  }))
  list(status="BASELINE_SAVED",config_version=cfg$version,metrics=metrics,changes=changes,
    gate=gate,starts=starts,replay=do.call(rbind,checks),curves=do.call(rbind,curves),
    calendar_peaks=do.call(rbind,annual),
    source_hashes=valid$source_hashes,phase23_beta_reference=valid$phase23_beta_reference)
}
p24_scenario <- function(score,inputs,clock,y,start,u=0,v=0,swap=NULL,padding=13L) {
  ix <- p24_window(clock,y);at <- ix[1L]-1L
  state <- unlist(score$states[at,c("S","I","R")],use.names=TRUE)
  mass <- sum(state)
  if(is.null(swap)) state["S"] <- state["S"]*exp(u) else {
    if(!is.na(swap["S"])) state["S"] <- swap["S"]
    if(!is.na(swap["I"])) state["I"] <- swap["I"]
  }
  state["R"] <- mass-state["S"]-state["I"]
  if(any(!is.finite(state)) || any(state<0)) return(list(status="INFEASIBLE_STATE"))
  end <- min(nrow(score$states),ix[52L]+padding)
  q <- c(rep(1,nrow(inputs$replay$data$time_markers)),score$q)
  mult <- q[at:end];mult[2:53] <- mult[2:53]*exp(v)
  z <- simulate_weekly_sirs(score$beta,state,p24_fixed(score,inputs),
    score$states$week[at:end],mult,initial_mass=mass)
  a <- z$expected_cases[2:53];top <- which(a==max(a));pk <- top[1L]
  list(status="VALID",states=z,start=state,selected_position=pk,
    selected_height=a[pk],ties=length(top),
    cog=if(sum(a)>0) sum((1:52)*a)/sum(a) else NA_real_,
    total=sum(a),last_supported=nrow(z)-1L)
}
p24_peaks <- function(z,id,reference,y,path,u,v,level) {
  a <- z$states$expected_cases;n <- length(a)
  policy <- list(tie_abs=1e-8,tie_rel=1e-10,small_height=0)
  local <- p4_weekly_peaks(a,seq_len(n),policy)$local
  if(!nrow(local)) return(data.frame())
  at <- as.integer(local$first)
  data.frame(scenario=id,reference=reference,season=y,path=path,u=u,v=v,level=level,
    peak_id=at,position=at-1L,height=local$height,
    left_height=ifelse(at>1L,a[pmax(1L,at-1L)],NA_real_),
    right_height=ifelse(local$last<n,a[pmin(n,local$last+1L)],NA_real_),
    plateau=local$plateau,support=!local$boundary,
    genuine=local$isolated,in_window=at>=2L & at<=53L,
    branch=NA_character_,stringsAsFactors=FALSE)
}
p24_build_grid <- function(cfg) {
  g <- cfg$log_grid
  if(length(g)%%2L!=1L || !any(abs(g)<1e-14) ||
     max(abs(g))>log(1.21) || max(abs(g+rev(g)))>1e-12) stop("Invalid preregistered grid")
  rbind(data.frame(path="S_ONLY",u=g,v=0,level=0L),
        data.frame(path="BETA_ONLY",u=0,v=g,level=0L),
        transform(expand.grid(u=g,v=g),path="JOINT",level=0L))
}
p24_refine_grid <- function(grid,selected,cfg) {
  added <- list()
  for(path in c("S_ONLY","BETA_ONLY")) {
    p <- grid[grid$path==path,]; x <- if(path=="S_ONLY") p$u else p$v
    for(level in seq_len(cfg$refine_levels)) {
      ord <- order(x);p <- p[ord,];x <- x[ord]
      z <- selected[match(paste(path,round(p$u,12),round(p$v,12)),
                          paste(selected$path,round(selected$u,12),round(selected$v,12))),]
      change <- which(diff(z$selected_position)!=0L | z$ties[-nrow(z)]>1L | z$ties[-1L]>1L)
      if(!length(change)) break
      mids <- (x[change]+x[change+1L])/2
      mids <- mids[!vapply(mids,function(t) any(abs(x-t)<1e-12),logical(1))]
      mids <- head(mids,cfg$max_refined_per_path-sum(p$level>0L))
      if(!length(mids)) break
      new <- data.frame(path=path,u=if(path=="S_ONLY") mids else 0,
        v=if(path=="BETA_ONLY") mids else 0,level=level)
      p <- rbind(p,new);x <- if(path=="S_ONLY") p$u else p$v
      added[[length(added)+1L]] <- new
    }
  }
  if(length(added)) do.call(rbind,added) else grid[FALSE,]
}
p24_interventions <- function(cfg,source_cfg,baseline_path) {
  base <- p24_read_stage(baseline_path,"baseline",cfg)
  if(!all(base$replay$passed)) stop("Blocking replay failed")
  inputs <- p22_inputs(source_cfg);clock <- p24_clock(inputs)
  scores <- p24_scores(inputs,source_cfg,cfg$replay_tolerance)$scores
  grid <- p24_build_grid(cfg)
  scenarios <- list();peaks <- list();curves <- list();counter <- 0L
  run <- function(nm,y,p,u,v,level,swap=NULL) {
    counter <<- counter+1L;id <- sprintf("scenario-%06d",counter)
    z <- p24_scenario(scores[[nm]],inputs,clock,y,NULL,u,v,swap,cfg$padding_weeks)
    row <- data.frame(scenario=id,reference=nm,season=y,path=p,u=u,v=v,level=level,
      status=z$status,selected_position=if(z$status=="VALID") z$selected_position else NA_integer_,
      selected_height=if(z$status=="VALID") z$selected_height else NA_real_,
      ties=if(z$status=="VALID") z$ties else NA_integer_,
      cog=if(z$status=="VALID") z$cog else NA_real_,
      total=if(z$status=="VALID") z$total else NA_real_,stringsAsFactors=FALSE)
    scenarios[[length(scenarios)+1L]] <<- row
    if(z$status=="VALID") {
      peaks[[length(peaks)+1L]] <<- p24_peaks(z,id,nm,y,p,u,v,level)
      curves[[length(curves)+1L]] <<- data.frame(scenario=id,week=1:52,expected=z$states$expected_cases[2:53])
    }
    invisible(row)
  }
  for(nm in names(scores)) for(y in cfg$seasons) {
    for(i in seq_len(nrow(grid))) run(nm,y,grid$path[i],grid$u[i],grid$v[i],grid$level[i])
    this <- do.call(rbind,scenarios);this <- this[this$reference==nm & this$season==y,]
    new <- p24_refine_grid(grid,this,cfg)
    if(nrow(new)) for(i in seq_len(nrow(new))) run(nm,y,new$path[i],new$u[i],new$v[i],new$level[i])
  }
  swaps <- list()
  for(nm in names(scores)) for(y in c(2015L,2017L)) {
    other <- if(y==2015L) 2017L else 2015L
    b <- base$starts[base$starts$reference==nm & base$starts$season==other,]
    for(kind in c("S","I","SI")) {
      sw <- c(S=if(grepl("S",kind)) b$S else NA_real_,
              I=if(grepl("I",kind)) b$I else NA_real_)
      z <- p24_scenario(scores[[nm]],inputs,clock,y,NULL,swap=sw,padding=cfg$padding_weeks)
      swaps[[length(swaps)+1L]] <- data.frame(reference=nm,baseline_season=y,
        donor_season=other,swap=kind,status=z$status,
        peak_season_week=if(z$status=="VALID") z$selected_position else NA_integer_,
        center_of_gravity=if(z$status=="VALID") z$cog else NA_real_,
        S=if(z$status=="VALID") z$start["S"] else NA_real_,
        I=if(z$status=="VALID") z$start["I"] else NA_real_,row.names=NULL)
    }
  }
  bind <- function(x) {x <- Filter(function(a)nrow(a)>0L,x);if(length(x)) do.call(rbind,x) else data.frame()}
  list(status="INTERVENTIONS_SAVED",config_version=cfg$version,
    grid=bind(scenarios),peaks=bind(peaks),curves=bind(curves),swaps=bind(swaps),
    coarse_grid=grid,refined_points=bind(lapply(scenarios,function(z)z[z$level>0L,])),
    starts=base$starts,baseline_gate=base$gate,baseline_metrics=base$metrics,
    baseline_changes=base$changes,baseline_curves=base$curves,
    calendar_peaks=base$calendar_peaks,
    replay=base$replay,source_hashes=base$source_hashes,
    phase23_beta_reference=base$phase23_beta_reference)
}
p24_link_group <- function(peaks,grid,cfg,prefix,sort_grid=TRUE) {
  if(sort_grid) {
    axis <- if(all(grid$path=="S_ONLY")) grid$u else grid$v
    grid <- grid[order(axis),]
  }
  next_id <- 0L;previous <- NULL
  for(k in seq_len(nrow(grid))) {
    id <- grid$scenario[k];ix <- which(peaks$scenario==id)
    if(!length(ix)) {previous <- NULL;next;}
    current <- peaks[ix,,drop=FALSE]
    labels <- current$branch
    if(!is.null(previous) && nrow(previous)) {
      f <- p10_assignment(previous,current,cfg$identity)
      r <- p10_assignment(previous,current,cfg$identity,reverse=TRUE)
      if(p10_identity_consistent(f,r)) {
        m <- f[f$status=="MATCHED",,drop=FALSE]
        for(j in seq_len(nrow(m))) {
          left <- match(m$left_peak[j],previous$peak_id)
          right <- match(m$right_peak[j],current$peak_id)
          if(!is.na(left) && !is.na(right) && is.na(labels[right])) labels[right] <- previous$branch[left]
        }
      }
    }
    for(j in which(is.na(labels))) {next_id <- next_id+1L;labels[j] <- paste0(prefix,"/b",next_id)}
    peaks$branch[ix] <- labels;current$branch <- labels;previous <- current
  }
  peaks
}
p24_identity <- function(peaks,grid,cfg) {
  if(!nrow(peaks)) return(peaks)
  for(nm in unique(grid$reference)) for(y in unique(grid$season)) {
    for(path in c("S_ONLY","BETA_ONLY")) {
      g <- grid[grid$reference==nm & grid$season==y & grid$path==path & grid$status=="VALID",]
      if(nrow(g)) peaks <- p24_link_group(peaks,g,cfg,paste(nm,y,path,sep="/"))
    }
    j <- grid[grid$reference==nm & grid$season==y & grid$path=="JOINT" & grid$status=="VALID",]
    zero <- j[abs(j$v)<1e-12,];zero <- zero[order(zero$u),]
    if(nrow(zero)) peaks <- p24_link_group(peaks,zero,cfg,paste(nm,y,"JOINT",sep="/"),FALSE)
    for(u in sort(unique(j$u))) for(sign in c(-1,1)) {
      g <- j[abs(j$u-u)<1e-12 & (if(sign<0) j$v<=0 else j$v>=0),]
      g <- g[order(abs(g$v)),]
      if(nrow(g)>1L) peaks <- p24_link_group(peaks,g,cfg,paste(nm,y,"JOINT",round(u,7),sign,sep="/"),FALSE)
    }
  }
  peaks
}
p24_selected <- function(grid,peaks) {
  grid$selected_branch <- NA_character_;grid$selection_status <- "UNRESOLVED"
  for(i in seq_len(nrow(grid))) {
    if(grid$status[i]!="VALID") {grid$selection_status[i] <- grid$status[i];next}
    if(grid$ties[i]>1L) {grid$selection_status[i] <- "TIE";next}
    p <- peaks[peaks$scenario==grid$scenario[i] & peaks$position==grid$selected_position[i] & peaks$genuine,]
    if(grid$selected_position[i] %in% c(1L,52L))
      grid$selection_status[i] <- "WINDOW_SELECTION"
    else if(nrow(p)==1L) {
      grid$selected_branch[i] <- p$branch
      grid$selection_status[i] <- "LOCAL_PEAK"
    } else grid$selection_status[i] <- "NO_ELIGIBLE_PEAK"
  }
  grid
}
p24_gaps <- function(grid,peaks) {
  rows <- list()
  for(i in seq_len(nrow(grid))) {
    z <- grid[i,];p <- peaks[peaks$scenario==z$scenario & peaks$in_window & peaks$genuine,]
    p <- p[order(-p$height),]
    rows[[i]] <- data.frame(scenario=z$scenario,reference=z$reference,season=z$season,
      path=z$path,u=z$u,v=z$v,level=z$level,
      branch1=if(nrow(p)) p$branch[1] else NA_character_,
      height1=if(nrow(p)) p$height[1] else NA_real_,
      branch2=if(nrow(p)>1L) p$branch[2] else NA_character_,
      height2=if(nrow(p)>1L) p$height[2] else NA_real_,
      signed_gap=if(nrow(p)>1L) p$height[1]-p$height[2] else NA_real_)
  }
  do.call(rbind,rows)
}
p24_edges <- function(grid,peaks) {
  rows <- list()
  groups <- split(grid,list(grid$reference,grid$season,grid$path,
                            ifelse(grid$path=="JOINT",round(grid$u,10),0)),drop=TRUE)
  for(g in groups) {
    x <- if(g$path[1]=="S_ONLY") g$u else g$v;g <- g[order(x),]
    if(nrow(g)<2L) next
    for(k in seq_len(nrow(g)-1L)) {
      a <- g[k,];b <- g[k+1L,]
      cls <- "UNRESOLVED";gap_left <- gap_right <- NA_real_
      if(any(c(a$selection_status,b$selection_status)=="WINDOW_SELECTION")) cls <- "WINDOW_SELECTION"
      else if(any(c(a$selection_status,b$selection_status)=="NO_ELIGIBLE_PEAK")) cls <- "NO_ELIGIBLE_PEAK"
      else if(any(c(a$selection_status,b$selection_status)=="TIE")) cls <- "TIE"
      else if(all(c(a$selection_status,b$selection_status)=="LOCAL_PEAK")) {
        if(a$selected_branch==b$selected_branch) cls <- if(a$selected_position==b$selected_position) "STABLE_SAME_BRANCH" else "DISPLACEMENT"
        else {
          pa <- peaks[peaks$scenario==a$scenario & peaks$in_window & peaks$genuine,]
          pb <- peaks[peaks$scenario==b$scenario & peaks$in_window & peaks$genuine,]
          aa <- pa$height[match(a$selected_branch,pa$branch)]
          ab <- pa$height[match(b$selected_branch,pa$branch)]
          ba <- pb$height[match(a$selected_branch,pb$branch)]
          bb <- pb$height[match(b$selected_branch,pb$branch)]
          if(all(is.finite(c(aa,ab,ba,bb)))) {
            gap_left <- aa-ab;gap_right <- ba-bb
            if(gap_left*gap_right<0) cls <- "DOMINANCE_EXCHANGE"
          }
        }
      }
      rows[[length(rows)+1L]] <- data.frame(reference=a$reference,season=a$season,
        path=a$path,left=a$scenario,right=b$scenario,
        x_left=if(a$path=="S_ONLY") a$u else a$v,
        x_right=if(b$path=="S_ONLY") b$u else b$v,
        classification=cls,left_branch=a$selected_branch,right_branch=b$selected_branch,
        signed_gap_left=gap_left,signed_gap_right=gap_right,
        peak_week_change=b$selected_position-a$selected_position,
        cog_change=b$cog-a$cog)
    }
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p24_interaction <- function(grid) {
  j <- grid[grid$path=="JOINT",];rows <- list()
  for(i in seq_len(nrow(j))) {
    z <- j[i,];a <- j[j$reference==z$reference & j$season==z$season & abs(j$u-z$u)<1e-12 & abs(j$v)<1e-12,]
    b <- j[j$reference==z$reference & j$season==z$season & abs(j$u)<1e-12 & abs(j$v-z$v)<1e-12,]
    o <- j[j$reference==z$reference & j$season==z$season & abs(j$u)<1e-12 & abs(j$v)<1e-12,]
    if(nrow(a)!=1L || nrow(b)!=1L || nrow(o)!=1L) stop("Incomplete interaction grid")
    four <- rbind(z,a,b,o)
    compatible <- all(four$selection_status=="LOCAL_PEAK") &&
      length(unique(four$selected_branch))==1L
    rows[[i]] <- data.frame(reference=z$reference,season=z$season,u=z$u,v=z$v,
      status=if(compatible) "SAME_ELIGIBLE_BRANCH" else "CATEGORICAL_BRANCH_OR_WINDOW",
      peak_interaction=if(compatible) z$selected_position-a$selected_position-b$selected_position+o$selected_position else NA_real_,
      cog_interaction=if(all(is.finite(four$cog))) z$cog-a$cog-b$cog+o$cog else NA_real_)
  }
  do.call(rbind,rows)
}
p24_adjudicate <- function(cfg,interventions_path) {
  x <- p24_read_stage(interventions_path,"interventions",cfg)
  peaks <- p24_identity(x$peaks,x$grid,cfg)
  grid <- p24_selected(x$grid,peaks)
  edges <- p24_edges(grid,peaks)
  coarse_edges <- p24_edges(grid[grid$level==0L,,drop=FALSE],peaks)
  gaps <- p24_gaps(grid,peaks)
  interaction <- p24_interaction(grid)
  mechanisms <- do.call(rbind,lapply(split(edges,list(edges$reference,edges$season,edges$path),drop=TRUE),function(z)
    data.frame(reference=z$reference[1],season=z$season[1],path=z$path[1],
      coarse_intervals=sum(coarse_edges$reference==z$reference[1] &
        coarse_edges$season==z$season[1] & coarse_edges$path==z$path[1]),
      refined_intervals=nrow(z),displacement=sum(z$classification=="DISPLACEMENT"),
      dominance_exchange=sum(z$classification=="DOMINANCE_EXCHANGE"),
      window_selection=sum(z$classification=="WINDOW_SELECTION"),
      unresolved=sum(z$classification %in% c("UNRESOLVED","TIE","NO_ELIGIBLE_PEAK")))))
  # The center-of-gravity interaction is a 52-week window statistic even when peak identity differs.
  x$grid <- grid;x$peaks <- peaks;x$edges <- edges;x$coarse_edges <- coarse_edges;x$gaps <- gaps
  x$interaction <- interaction;x$mechanisms <- mechanisms
  x$status <- "ADJUDICATED";x
}
p24_verdicts <- function(x,cfg) {
  rows <- list()
  for(nm in unique(x$baseline_gate$reference)) {
    gate <- x$baseline_gate[x$baseline_gate$reference==nm,]
    obs <- x$baseline_changes[x$baseline_changes$reference==nm & x$baseline_changes$series=="OBSERVED",]
    fit <- x$baseline_changes[x$baseline_changes$reference==nm & x$baseline_changes$series=="FITTED",]
    rows[[length(rows)+1L]] <- data.frame(reference=nm,claim="A_OBSERVED_TRANSITION",
      verdict=if(gate$passed) "REPRODUCED_WITHIN_PREDECLARED_TOLERANCES" else "NOT_REPRODUCED",
      detail=sprintf("peak change observed %.1f, fitted %.1f; CoG change observed %.2f, fitted %.2f; errors %.1f and %.2f weeks",
        obs$peak_change,fit$peak_change,obs$cog_change,fit$cog_change,
        gate$peak_change_error,gate$cog_change_error))
    for(path in c("S_ONLY","BETA_ONLY","JOINT")) {
      e <- x$mechanisms[x$mechanisms$reference==nm & x$mechanisms$path==path,]
      rows[[length(rows)+1L]] <- data.frame(reference=nm,claim=paste0("B_",path),
        verdict=if(any(e$dominance_exchange>0)) "DOMINANCE_EXCHANGE_PRESENT" else if(any(e$window_selection>0)) "WINDOW_SELECTION_PRESENT" else if(any(e$displacement>0)) "DISPLACEMENT_PRESENT" else "UNRESOLVED_OR_STABLE",
        detail=sprintf("Across three seasons: displacement %d, dominance exchange %d, window selection %d, unresolved %d; interval counts include refinement",
          sum(e$displacement),sum(e$dominance_exchange),sum(e$window_selection),sum(e$unresolved)))
    }
    base <- x$baseline_metrics[x$baseline_metrics$reference==nm & x$baseline_metrics$season==2015L & x$baseline_metrics$series=="FITTED",]
    sw <- x$swaps[x$swaps$reference==nm & x$swaps$baseline_season==2015L,]
    si <- sw[sw$swap=="SI",]
    movement <- si$peak_season_week-base$peak_season_week
    cog_movement <- si$center_of_gravity-base$center_of_gravity
    enough <- length(movement)==1L && is.finite(movement) && is.finite(cog_movement) &&
      abs(movement-obs$peak_change)<=cfg$observed_gate$max_peak_change_error_weeks &&
      abs(cog_movement-obs$cog_change)<=cfg$observed_gate$max_cog_change_error_weeks
    rows[[length(rows)+1L]] <- data.frame(reference=nm,claim="C_START_STATE_SWAP",
      verdict=if(enough) "SUFFICIENT_MAGNITUDE_IN_CONTROLLED_SWAP" else "INSUFFICIENT_OR_UNRESOLVED",
      detail=sprintf("2015/16 baseline to fitted 2017/18 S/I swap: peak movement %.1f versus observed %.1f weeks; CoG movement %.2f versus observed %.2f weeks",
        movement,obs$peak_change,cog_movement,obs$cog_change))
  }
  do.call(rbind,rows)
}
p24_plot <- function(file,expr,width=10,height=6) {
  grDevices::pdf(file,width=width,height=height)
  tryCatch(force(expr),finally=grDevices::dev.off())
}
p24_report <- function(cfg,adjudicate_path,out) {
  x <- p24_read_stage(adjudicate_path,"adjudicate",cfg)
  if(!identical(x$status,"ADJUDICATED") || !all(x$replay$passed)) stop("Adjudicated replay required")
  verdicts <- p24_verdicts(x,cfg)
  write <- function(z,name) utils::write.csv(z,file.path(out,paste0(name,".csv")),row.names=FALSE)
  write(x$baseline_metrics,"observed_vs_fitted_seasons")
  write(x$baseline_curves,"observed_vs_fitted_weekly")
  write(x$calendar_peaks,"calendar_year_peaks")
  write(x$starts,"season_start_states")
  write(x$grid,"intervention_grid")
  write(x$peaks,"all_peak_branches")
  write(x$gaps,"competitor_height_gaps")
  write(x$edges[,c("reference","season","path","left","right","left_branch","right_branch",
                  "signed_gap_left","signed_gap_right","classification")],"paired_branch_gaps")
  write(x$mechanisms,"mechanism_by_path")
  write(x$swaps,"state_swaps")
  write(x$interaction,"interaction")
  write(verdicts,"claim_verdicts")
  write(x$replay,"replay_checks")
  write(x$baseline_gate,"observed_transition_gate")
  write(x$phase23_beta_reference,"phase23_beta_reference")
  write(x$edges,"interval_classifications")
  write(x$coarse_edges,"coarse_interval_classifications")
  write(x$baseline_changes,"timing_changes")
  write(data.frame(path=names(x$source_hashes),md5=unname(x$source_hashes)),"provenance")
  p24_plot(file.path(out,"observed_vs_fitted_seasons.pdf"),{
    op <- graphics::par(mfrow=c(3,1),mar=c(4,4,2,1))
    for(y in cfg$seasons) {
      z <- x$baseline_curves[x$baseline_curves$reference=="HIST" & x$baseline_curves$season==y,]
      j <- x$baseline_curves[x$baseline_curves$reference=="JOINT" & x$baseline_curves$season==y,]
      m <- x$baseline_metrics[x$baseline_metrics$reference=="HIST" & x$baseline_metrics$season==y,]
      q <- x$baseline_metrics[x$baseline_metrics$reference=="JOINT" & x$baseline_metrics$season==y & x$baseline_metrics$series=="FITTED",]
      graphics::plot(z$season_week,z$observed,type="l",col="black",xlab="Season week (week 27 to week 26)",ylab="Cases",main=paste(y,y+1,sep="/"),ylim=range(c(z$observed,z$expected,j$expected)))
      graphics::lines(z$season_week,z$expected,col="steelblue")
      graphics::lines(j$season_week,j$expected,col="firebrick")
      graphics::abline(v=m$peak_season_week,col=c("black","steelblue"),lty=2)
      graphics::abline(v=m$center_of_gravity,col=c("black","steelblue"),lty=3)
      graphics::abline(v=q$peak_season_week,col="firebrick",lty=2)
      graphics::abline(v=q$center_of_gravity,col="firebrick",lty=3)
    }
    graphics::par(op)
  })
  p24_plot(file.path(out,"S_only_vs_beta_only_peak_paths.pdf"),{
    op <- graphics::par(mfrow=c(2,1))
    for(path in c("S_ONLY","BETA_ONLY")) {
      z <- x$grid[x$grid$reference=="HIST" & x$grid$season==2015L & x$grid$path==path,]
      a <- if(path=="S_ONLY") z$u else z$v
      graphics::plot(a,z$selected_position,pch=16,xlab=paste(path,"log multiplier"),ylab="Selected season week",main=paste(path,"2015/16"),ylim=c(0,65))
      graphics::abline(h=c(1,52),lty=3)
      p <- x$peaks[x$peaks$reference=="HIST" & x$peaks$season==2015L & x$peaks$path==path,]
      for(id in unique(p$branch)) {
        b <- p[p$branch==id,];b <- b[order(if(path=="S_ONLY") b$u else b$v),]
        graphics::lines(if(path=="S_ONLY") b$u else b$v,b$position,col="steelblue")
        if(nrow(b)) graphics::text(if(path=="S_ONLY") tail(b$u,1L) else tail(b$v,1L),
          tail(b$position,1L),labels=sub("^.*/","",id),cex=.55,pos=4)
      }
    }
    graphics::par(op)
  })
  p24_plot(file.path(out,"competitor_heights_and_gaps.pdf"),{
    z <- x$gaps[x$gaps$reference=="HIST" & x$gaps$season==2015L & x$gaps$path=="S_ONLY",]
    z <- z[order(z$u),];op <- graphics::par(mfrow=c(2,1))
    graphics::plot(z$u,z$height1,type="b",xlab="Log S multiplier",ylab="Peak height",main="Leading eligible branches")
    graphics::lines(z$u,z$height2,type="b",col="steelblue")
    e <- x$edges[x$edges$reference=="HIST" & x$edges$season==2015L & x$edges$path=="S_ONLY" &
                 is.finite(x$edges$signed_gap_left),]
    if(nrow(e)) {
      graphics::plot(range(c(e$x_left,e$x_right)),range(c(e$signed_gap_left,e$signed_gap_right,0)),
        type="n",xlab="Log S multiplier",ylab="Left winner minus right winner height")
      graphics::segments(e$x_left,e$signed_gap_left,e$x_right,e$signed_gap_right,
        col=ifelse(e$classification=="DOMINANCE_EXCHANGE","firebrick","grey50"))
      graphics::abline(h=0,lty=2)
    } else graphics::plot.new()
    graphics::par(op)
  })
  p24_plot(file.path(out,"timing_change_comparison.pdf"),{
    z <- verdicts[verdicts$reference %in% c("HIST","JOINT") & verdicts$claim=="A_OBSERVED_TRANSITION",]
    refs <- z$reference;vals <- matrix(NA_real_,nrow=length(refs),ncol=3,
      dimnames=list(refs,c("Observed","Fitted","S/I swap")))
    for(i in seq_along(refs)) {
      a <- x$baseline_changes[x$baseline_changes$reference==refs[i],]
      vals[i,1:2] <- a$peak_change[match(c("OBSERVED","FITTED"),a$series)]
      b <- x$baseline_metrics[x$baseline_metrics$reference==refs[i] & x$baseline_metrics$season==2015L & x$baseline_metrics$series=="FITTED",]
      s <- x$swaps[x$swaps$reference==refs[i] & x$swaps$baseline_season==2015L & x$swaps$swap=="SI",]
      vals[i,3] <- s$peak_season_week-b$peak_season_week
    }
    graphics::barplot(t(vals),beside=TRUE,legend.text=TRUE,ylab="2015/16 to 2017/18 peak change (weeks)")
  })
  lines <- c("# Phase 24 saved-model audit","",
    paste("Configuration:",cfg$version,"; deterministic seed:",cfg$seed),
    paste("Replay tolerance:",format(cfg$replay_tolerance),"absolute counts/states."),
    "2015-W53 was merged into 2015-W52 in the accepted prepared historical data.",
    "Each primary season has 52 prepared weeks, from week 27 through week 26.","",
    "## Claim verdicts",capture.output(print(verdicts,row.names=FALSE)),"",
    "The beta-only path is hypothetical: fitted seasonal beta is shared across 2013-2019.",
    "Counterfactual mechanisms do not identify the cause of the observed transition.",
    "Grid intervals and optimizer endpoints are neither independent events nor confidence intervals.","",
    "## Source hashes",capture.output(print(x$source_hashes)))
  writeLines(lines,file.path(out,"phase24_report.md"))
  list(status="REPORT_BUILT",config_version=cfg$version,verdicts=verdicts,
    adjudicate_hash=unname(tools::md5sum(adjudicate_path)))
}
