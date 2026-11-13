# Deterministic counterfactuals from an explicitly accepted Phase 7C lambda=1 endpoint.
p8_defaults <- function() list(source_lambda=1,replay_tolerance=1e-7,
  time=list(horizon_rows=NULL,years=NULL,include_incomplete=FALSE),
  interventions=list(s_values=NULL,s_multipliers=c(.8,1,1.2),a_values=c(.8,1,1.2),
    infected_policy="fixed_baseline_fraction",joint=FALSE,max_scenarios=100L),
  annual=list(enabled=TRUE,mean_beta0=NULL,amplitude=.3,peak_week=1L),
  peaks=list(periodic=FALSE,tie_abs=1e-8,tie_rel=1e-6,small_height=1e-6,
    min_prominence=0,min_relative_prominence=.01,min_displacement=8,
    match_distance=12,max_log_height_change=log(2),height_weight=2,identity_margin=.5,
    max_interval=list(s=.01,a=.25)),
  refinement=list(enabled=FALSE,levels=2L,max_points=12L))
p8_critical_files <- function(root) normalizePath(file.path(root,"R",c("sirs_core.R","nb2.R","cyclic_spline.R",
  "parameter_transforms.R","penalized_objective.R","data_adapter.R")),mustWork=TRUE)
p8_source <- function(source,data,s,root=".") {
  if(!identical(source$type,"phase7c_tangent_refit") || !isTRUE(source$gate$passed) ||
    !isTRUE(source$config$independent) || !identical(s$source_lambda,1)) stop("Accepted independent Phase 7C lambda=1 source required")
  if(!identical(source$config,source$source_config)) stop("Saved scientific configuration mismatch")
  r <- source$records[["1"]]
  if(is.null(r) || r$lambda!=1 || !isTRUE(r$accepted) || !isTRUE(r$optimizer_terminated) ||
    !isTRUE(r$independent_replay) || !identical(r$candidate_kind,"OPTIMIZER_ENDPOINT") ||
    !identical(r$status,"IMPROVED_OPTIMIZER_ENDPOINT_ACCEPTED") ||
    !identical(r$retained_coordinates,r$candidate_coordinates)) stop("Reject unaccepted/abnormal lambda=1 diagnostic candidate; no fallback")
  paths <- unique(c(p8_critical_files(root),normalizePath(p3_inputs(source$config),mustWork=TRUE)))
  old <- source$identity[paths]
  if(anyNA(old) || !identical(unname(old),unname(tools::md5sum(paths)))) stop("Scientific code or input identity mismatch")
  cfg <- source$config;cfg$spline$lambda <- 1;p <- p6_problem(cfg,data,"SPL")
  if(p$basis$K!=8L || !isTRUE(p$map$joint_mode)) stop("K=8 joint spline source required")
  q <- r$retained_coordinates;check_optimizer_coordinates(q,p$map)
  e <- evaluate_candidate(q,p);old <- r$retained_evaluation
  tol <- s$replay_tolerance
  if(!e$valid || !isTRUE(old$valid) || !p6r_close(c(e$loglik,e$penalty,e$objective),c(old$loglik,old$penalty,old$objective),tol) ||
    !p6r_close(unlist(e$parameters),unlist(old$parameters),tol) || !p6r_close(e$beta_by_week,old$beta_by_week,tol) ||
    !p6r_close(as.numeric(as.matrix(e$states)),as.numeric(as.matrix(old$states)),tol)) stop("Accepted retained endpoint baseline replay mismatch")
  list(parameters=e$parameters,beta=e$beta_by_week,evaluation=e,problem=p,data=data,coordinates=q,
    source_hash=source$artifact_hash,source_identity=source$identity,source_config=source$config,
    source_ancestry=source$source_hash,replay="ACCEPTED_RETAINED_TRAJECTORY_MATCH",
    warnings=c("Acceptance is not an identifiability or global-optimum guarantee",
      paste("Source numerical_convergence_supported:",r$numerical_convergence_supported)))
}
p8_settings <- function(s,baseline) {
  if(!identical(s$source_lambda,1)) stop("Phase 8 baseline is lambda=1 only")
  scalar_number(s$replay_tolerance,"replay tolerance",0,strict_lower=TRUE)
  p <- baseline$parameters;n <- length(baseline$data$cases)
  h <- if(is.null(s$time$horizon_rows)) n else s$time$horizon_rows
  integer_scalar(h,"horizon rows",2L);if(h>n) stop("No population/calendar extrapolation beyond validated input support")
  values <- if(is.null(s$interventions$s_values)) p$S0*s$interventions$s_multipliers else s$interventions$s_values
  finite_numeric(values,"s grid");finite_numeric(s$interventions$a_values,"a grid")
  if(any(values<=0 | values>=1-p$I0) || any(diff(values)<=0) || !any(values==p$S0)) stop("Increasing interior s grid must contain exact saved S0")
  a <- s$interventions$a_values
  if(any(a<=0) || any(diff(a)<=0) || !any(a==1)) stop("Increasing positive a grid must contain exact 1")
  if(!identical(s$interventions$infected_policy,"fixed_baseline_fraction")) stop("Only explicitly fixed baseline I0 policy implemented")
  for(flag in list(s$time$include_incomplete,s$interventions$joint,s$annual$enabled,s$peaks$periodic,s$refinement$enabled))
    if(!is.logical(flag) || length(flag)!=1L || is.na(flag)) stop("Explicit logical experiment switches required")
  integer_scalar(s$interventions$max_scenarios,"scenario cap",1L)
  integer_scalar(s$refinement$levels,"refinement levels",0L);integer_scalar(s$refinement$max_points,"refinement point cap",0L)
  for(name in c("tie_abs","tie_rel","small_height","min_prominence","min_relative_prominence","min_displacement",
    "match_distance","max_log_height_change","height_weight","identity_margin")) scalar_number(s$peaks[[name]],name,0)
  if(s$peaks$match_distance>26 || s$peaks$min_displacement>26 || s$peaks$min_relative_prominence>1) stop("Invalid circular/prominence threshold")
  for(axis in c("s","a")) scalar_number(s$peaks$max_interval[[axis]],paste(axis,"max interval"),0,strict_lower=TRUE)
  markers <- p3_markers(baseline$data)[seq_len(h),,drop=FALSE]
  years <- sort(unique(markers$year))
  if(!is.null(s$time$years)) {
    finite_numeric(s$time$years,"analysis years")
    if(any(!s$time$years %in% years) || anyDuplicated(s$time$years)) stop("Analysis years outside horizon")
    years <- sort(s$time$years)
  }
  complete <- vapply(years,function(y) identical(as.integer(markers$week[markers$year==y & seq_len(h)>1]),1:52),logical(1))
  if(!any(complete | s$time$include_incomplete)) stop("No complete annual incidence windows; choose a longer supported horizon")
  list(horizon=h,s_values=values,a_values=a,years=years,complete=complete,settings=s)
}
p8_profile <- function(baseline,s,kind="fitted") {
  if(kind=="fitted") return(list(beta=baseline$beta,label="FITTED_INPUT_FIXED_SHAPE",annual_validation=NULL))
  a <- s$annual;scalar_number(a$amplitude,"annual amplitude",0,strict_lower=TRUE)
  integer_scalar(a$peak_week,"annual peak week",1L);if(a$peak_week>52) stop("Annual phase must be an integer week 1..52")
  mean <- if(is.null(a$mean_beta0)) baseline$parameters$beta0 else a$mean_beta0
  scalar_number(mean,"hypothetical annual mean scale",0,strict_lower=TRUE)
  basis <- new_fourier_basis(1L);phase <- 2*pi*(a$peak_week-1)/52
  coefficients <- a$amplitude*c(cos(phase),sin(phase))
  profile <- fourier_profile(coefficients,mean,basis)
  fit <- list(evaluation=list(parameters=list(beta0=mean,free_coordinates=coefficients)),map=list(fourier=TRUE),basis=basis)
  peaks <- p5_beta_peaks(fit,baseline$source_config$peak_detection)
  ties <- which(abs(profile$beta_by_week-max(profile$beta_by_week))<=1e-10*max(profile$beta_by_week))
  if(any(profile$beta_by_week<=0) || peaks$count!=1L || length(ties)!=1L || ties!=a$peak_week)
    stop("Annual single continuous/discrete maximum not numerically verified; no amplitude/phase tuning")
  list(beta=profile$beta_by_week,label="HYPOTHETICAL_ANNUAL_NOT_FITTED",
    annual_validation=data.frame(mean_beta0=mean,b_raw=mean*exp(-profile$log_normalizer),amplitude=a$amplitude,
      peak_week=a$peak_week,continuous_maxima=peaks$count,weekly_maxima=length(ties),
      definition="A>0 cosine log forcing: one continuous maximum analytically; integer phase ensures unique weekly maximum"))
}
p8_weekly <- function(baseline,check,profile,susceptible,a) {
  p <- baseline$parameters;N <- baseline$data$fixed_parameters$pop
  scalar_number(susceptible,"counterfactual s",0,1-p$I0,strict_lower=TRUE)
  if(susceptible>=1-p$I0) stop("Counterfactual R0 must be positive")
  scalar_number(a,"counterfactual a",0,strict_lower=TRUE)
  h <- check$horizon;weeks <- baseline$data$seasonal_week[seq_len(h)]
  fixed <- list(pop=N,mu=p$mu,gamma=p$gamma,omega=p$omega,delta=p$delta,rho=p$rho)
  beta <- a*profile$beta
  z <- simulate_weekly_sirs(beta,make_initial_state(susceptible,p$I0,N),fixed,weeks,baseline$data$npieff[seq_len(h)])
  z$year <- baseline$data$time_markers$year[seq_len(h)];z$initialization_only <- seq_len(h)==1L
  z$beta <- beta[z$week];z$npieff <- baseline$data$npieff[seq_len(h)]
  mass <- z$S+z$I+z$R;z$state_mass <- mass;z$population_denominator <- N
  z$susceptible_fraction_actual_mass <- z$S/mass;z$susceptible_fraction_denominator <- z$S/N
  z$delta_I <- c(NA_real_,diff(z$I));z$infectious_removals <- z$incidence-z$delta_I
  z$recoveries <- z$infectious_removals-z$deaths_I
  z$susceptible_exits <- z$incidence+z$deaths_S;z$susceptible_exits[1] <- NA_real_
  z$recovered_exits <- c(NA_real_,head(z$R,-1))+z$recoveries-z$R
  z$waning_inflow <- z$recovered_exits-z$deaths_R
  z$force_of_infection <- z$beta*(c(NA_real_,head(z$I,-1))+p$omega)/N*z$npieff
  z$effective_transmission_exposure <- z$beta*c(NA_real_,head(z$S,-1))/N*z$npieff
  z$incidence_change <- c(NA_real_,diff(z$incidence));z$incidence_change[2] <- NA_real_
  z$incidence_share_in_year <- NA_real_
  for(y in unique(z$year)) {
    ids <- which(z$year==y & !z$initialization_only);total <- sum(z$incidence[ids])
    if(total>0) z$incidence_share_in_year[ids] <- z$incidence[ids]/total
  }
  z
}
p8_discrete_peaks <- function(y,weeks,policy,periodic=policy$periodic) {
  finite_numeric(y,"incidence",length(weeks));finite_numeric(weeks,"weeks",length(y))
  if(any(y<0) || any(diff(weeks)!=1)) stop("Consecutive nonnegative weekly incidence required")
  circular <- isTRUE(periodic) && identical(as.integer(weeks),1:52)
  tol <- p4_tolerance(y,policy);height <- max(y);global <- weeks[height-y<=tol]
  empty <- data.frame(peak_id=integer(),week=numeric(),x=numeric(),height=numeric(),plateau=logical(),boundary=logical(),
    prominence=numeric(),relative_prominence=numeric(),width_weeks=numeric(),retained=logical())
  if(max(y)-min(y)<=tol) return(list(local=empty,global_week=global[1],global_ties=global,
    cumulative=sum(y),flat=TRUE,small=height<=policy$small_height,periodic=circular,tolerance=tol,seam_gap=abs(y[1]-tail(y,1))))
  raw <- if(circular) p4_weekly_peaks(rep(y,3L),seq_len(3L*length(y)),policy)$local else p4_weekly_peaks(y,seq_along(y),policy)$local
  m <- length(y);rows <- list()
  if(circular) raw <- raw[raw$last>=m+1 & raw$first<=2*m,,drop=FALSE]
  used <- numeric()
  for(i in seq_len(nrow(raw))) {
    ids <- seq.int(raw$first[i],raw$last[i]);ids <- unique((ids-1L)%%m+1L)
    at <- ids[which.min(weeks[ids])];week <- weeks[at]
    if(week %in% used) next
    used <- c(used,week);value <- max(y[ids])
    rows[[length(rows)+1L]] <- data.frame(peak_id=length(rows)+1L,week=week,x=week-1,height=value,
      plateau=length(ids)>1L,boundary=any(ids %in% c(1L,m)),at=at)
  }
  local <- if(length(rows)) do.call(rbind,rows) else NULL
  if(is.null(local)) local <- empty else {
    for(i in seq_len(nrow(local))) {
      at <- local$at[i];higher <- local$at[local$height>local$height[i]+tol]
      if(circular) {
        left <- if(length(higher)) min((at-higher)%%m) else m
        right <- if(length(higher)) min((higher-at)%%m) else m
        lo <- min(y[((at-(0:left)-1)%%m)+1]);hi <- min(y[((at+(0:right)-1)%%m)+1])
      } else {
        left <- if(any(higher<at)) max(higher[higher<at]) else 1L
        right <- if(any(higher>at)) min(higher[higher>at]) else m
        lo <- min(y[left:at]);hi <- min(y[at:right])
      }
      prom <- max(0,local$height[i]-max(lo,hi));level <- local$height[i]-prom/2
      width <- NA_real_
      if(prom>tol) {
        walk <- function(direction) {
          for(k in seq_len(m)) {
            j <- at+direction*k
            if(!circular && (j<1 || j>m)) return(NA_integer_)
            j <- (j-1L)%%m+1L
            if(y[j]<level) return(k)
          }
          NA_integer_
        }
        l <- walk(-1L);r <- walk(1L);if(!is.na(l) && !is.na(r)) width <- min(m,l+r-1L)
      }
      local$prominence[i] <- prom;local$relative_prominence[i] <- prom/local$height[i];local$width_weeks[i] <- width
    }
    local$retained <- local$prominence>=policy$min_prominence & local$relative_prominence>=policy$min_relative_prominence
    local$at <- NULL
  }
  list(local=local,global_week=global[1],global_ties=global,cumulative=sum(y),flat=FALSE,
    small=height<=policy$small_height,periodic=circular,tolerance=tol,seam_gap=abs(y[1]-tail(y,1)))
}
p8_match <- function(left,right,policy) {
  L <- left$local[left$local$retained,,drop=FALSE];R <- right$local[right$local$retained,,drop=FALSE]
  distance <- abs(outer(L$x,R$x,`-`));distance <- pmin(distance,52-distance)
  ratio <- abs(outer(log(L$height),log(R$height),`-`))
  allowed <- ratio<=policy$max_log_height_change
  cost <- distance+policy$height_weight*ratio
  result <- tryCatch(p5_circular_peak_match(L$x,R$x,policy$match_distance,edge_allowed=allowed,edge_cost=cost),error=function(e) e)
  if(inherits(result,"error")) return(list(status="MATCH_UNAVAILABLE",reason=conditionMessage(result),table=data.frame(),L=L,R=R,ambiguous=TRUE))
  result$left_peak <- L$peak_id[match(result$reference_x,L$x)];result$right_peak <- R$peak_id[match(result$target_x,R$x)]
  result$log_height_change <- ratio[cbind(match(result$reference_x,L$x),match(result$target_x,R$x))]
  ambiguous <- FALSE
  for(i in seq_len(nrow(L))) {
    vals <- sort(cost[i,allowed[i,] & distance[i,]<=policy$match_distance])
    if(length(vals)>1 && vals[2]-vals[1]<=policy$identity_margin) ambiguous <- TRUE
  }
  for(j in seq_len(nrow(R))) {
    vals <- sort(cost[allowed[,j] & distance[,j]<=policy$match_distance,j])
    if(length(vals)>1 && vals[2]-vals[1]<=policy$identity_margin) ambiguous <- TRUE
  }
  # A global assignment tie can exist even when every row/column has a unique cheapest edge.
  paired <- which(result$status=="MATCHED")
  total_cost <- function(t) {
    at <- which(t$status=="MATCHED")
    if(!length(at)) return(0)
    sum(cost[cbind(match(t$reference_x[at],L$x),match(t$target_x[at],R$x))])
  }
  best_cost <- total_cost(result)
  for(k in paired) {
    mask <- allowed
    mask[match(result$reference_x[k],L$x),match(result$target_x[k],R$x)] <- FALSE
    alt <- p5_circular_peak_match(L$x,R$x,policy$match_distance,edge_allowed=mask,edge_cost=cost)
    if(sum(alt$status=="MATCHED")==length(paired) && total_cost(alt)-best_cost<=policy$identity_margin)
      ambiguous <- TRUE
  }
  list(status="MATCHED_WITH_UNMATCHED_RECORDED",table=result,L=L,R=R,ambiguous=ambiguous,reason="")
}
p8_classify <- function(left,right,interval,axis,policy) {
  match <- p8_match(left,right,policy);L <- match$L;R <- match$R;t <- match$table
  distance <- min(abs(left$global_week-right$global_week),52-abs(left$global_week-right$global_week))
  Dleft <- Dright <- NA_real_;mechanism <- "AMBIGUOUS_UNRESOLVED";reason <- "Peak identity unresolved"
  good <- if(nrow(t)) t[t$status=="MATCHED",,drop=FALSE] else t
  unique_globals <- length(left$global_ties)==1L && length(right$global_ties)==1L
  lwin <- L$peak_id[L$week==left$global_week];rwin <- R$peak_id[R$week==right$global_week]
  switching <- displacement <- FALSE
  if(length(lwin)==1L && length(rwin)==1L && nrow(good)) {
    j <- good$right_peak[good$left_peak==lwin];i <- good$left_peak[good$right_peak==rwin]
    if(length(i)==1L && length(j)==1L) {
      displacement <- j==rwin
      if(!displacement) {
        Dleft <- L$height[L$peak_id==lwin]-L$height[L$peak_id==i]
        Dright <- R$height[R$peak_id==j]-R$height[R$peak_id==rwin]
        switching <- Dleft>0 && Dright<0
      }
    }
  }
  unmatched <- if(nrow(t)) sum(t$status!="MATCHED") else nrow(L)+nrow(R)
  boundary <- any(L$boundary | L$plateau) || any(R$boundary | R$plateau)
  if(left$flat || right$flat || left$small || right$small) reason <- "Flat/tiny incidence: timing unstable"
  else if(interval>policy$max_interval[[axis]]) reason <- "Parameter interval exceeds prespecified maximum"
  else if(!unique_globals) reason <- "Global tie or near-tie; preserve tied weeks"
  else if(boundary) reason <- "Annual boundary/plateau can alter peak topology"
  else if(match$ambiguous || match$status=="MATCH_UNAVAILABLE") reason <- "Matching ambiguity or unavailable matching"
  else if(nrow(left$local)!=nrow(right$local) && nrow(L)!=nrow(R) && unmatched>0) {
    mechanism <- "PEAK_CREATION_ANNIHILATION";reason <- "Discrete interior local-maximum count change; continuous bifurcation not certified"
  } else if(unmatched>0) reason <- "Unmatched/threshold-limited branches; do not force displacement"
  else if(switching && distance>=policy$min_displacement) {
    mechanism <- "PEAK_SWITCHING";reason <- "Two persistent branches reverse height ordering and global winner"
  } else if(displacement && distance>=policy$min_displacement) {
    mechanism <- "PEAK_DISPLACEMENT";reason <- "Same matched global branch moves on a finite grid; not a discontinuity theorem"
  } else {mechanism <- "NO_ABRUPT_CANDIDATE";reason <- "Below displacement threshold or no global branch exchange"}
  list(row=data.frame(mechanism=mechanism,reason=reason,global_circular_displacement=distance,parameter_interval=interval,
    D_left=Dleft,D_right=Dright,persistent_branch_exchange=switching,raw_peaks_left=nrow(left$local),raw_peaks_right=nrow(right$local),
    retained_peaks_left=nrow(L),retained_peaks_right=nrow(R),unmatched=unmatched,
    matching_ambiguous=match$ambiguous,unique_global_maxima=unique_globals,
    evidence="Finite-grid weekly candidate; not a mathematical discontinuity or causal claim"),matching=match)
}
p8_evaluate <- function(baseline,check,profile,susceptible,a,id,axis,stage="coarse",level=0L,parent="") {
  tryCatch({
    z <- p8_weekly(baseline,check,profile,susceptible,a);annual <- list();rows <- list()
    for(year in check$years) {
      at <- which(z$year==year & !z$initialization_only);complete <- identical(as.integer(z$week[at]),1:52)
      included <- complete || check$settings$time$include_incomplete
      peaks <- if(included) p8_discrete_peaks(z$incidence[at],z$week[at],check$settings$peaks) else NULL
      if(!is.null(peaks)) {peaks$complete <- complete;annual[[as.character(year)]] <- peaks}
      rows[[length(rows)+1L]] <- data.frame(year=year,n_weeks=length(at),complete=complete,
        status=if(included) "ANALYZED" else "INCOMPLETE_EXCLUDED",primary_outcome="incident_infections",
        global_week=if(included) peaks$global_week else NA_real_,
        global_ties=if(included) paste(peaks$global_ties,collapse=";") else "",
        global_height=if(included) max(z$incidence[at]) else NA_real_,
        raw_peak_count=if(included) nrow(peaks$local) else NA_integer_,
        retained_peak_count=if(included) sum(peaks$local$retained) else NA_integer_,
        cumulative_incidence=if(included) peaks$cumulative else NA_real_,
        periodic_analysis=if(included) peaks$periodic else FALSE,
        flat=if(included) peaks$flat else NA,seam_gap=if(included) peaks$seam_gap else NA_real_)
    }
    z$analysis_window_included <- z$year %in% as.integer(names(annual)) & !z$initialization_only
    list(id=id,profile=profile$kind,axis=axis,stage=stage,level=level,parent=parent,s=susceptible,a=a,
      status="VALID",reason="",trajectory=z,annual=annual,annual_summary=p3_bind_rows(rows))
  },error=function(e) list(id=id,profile=profile$kind,axis=axis,stage=stage,level=level,parent=parent,s=susceptible,a=a,
    status="FAILED",reason=conditionMessage(e),annual=list(),trajectory=NULL))
}
p8_edge <- function(left,right,axis,year,policy,stage="coarse",level=0L) {
  theta <- function(r) if(axis=="s") r$s else r$a
  L <- left$annual[[as.character(year)]];R <- right$annual[[as.character(year)]]
  if(left$status!="VALID" || right$status!="VALID" || is.null(L) || is.null(R)) {
    row <- data.frame(mechanism="AMBIGUOUS_UNRESOLVED",reason="Failed scenario or excluded/incomplete window",
      global_circular_displacement=NA_real_,parameter_interval=theta(right)-theta(left))
    match <- list(table=data.frame(),status="UNAVAILABLE",ambiguous=TRUE)
  } else {
    value <- p8_classify(L,R,theta(right)-theta(left),axis,policy);row <- value$row;match <- value$matching
    if(!isTRUE(L$complete) || !isTRUE(R$complete)) {row$mechanism <- "AMBIGUOUS_UNRESOLVED";row$reason <- "Incomplete annual window"}
  }
  row$left <- left$id;row$right <- right$id;row$profile <- left$profile;row$axis <- axis;row$year <- year
  row$theta_left <- theta(left);row$theta_right <- theta(right);row$fixed_other <- if(axis=="s") left$a else left$s
  row$stage <- stage;row$level <- level
  list(row=row,matching=match,left=left$id,right=right$id,axis=axis,year=year)
}
p8_series <- function(records) {
  ids <- Filter(function(r) r$axis %in% c("s","a","joint"),records);groups <- list()
  for(r in ids) for(axis in if(r$axis=="joint") c("s","a") else r$axis) {
    key <- paste(r$profile,axis,format(if(axis=="s") r$a else r$s,digits=17),sep=":")
    groups[[key]] <- c(groups[[key]],r$id)
  }
  groups
}
p8_track <- function(records,groups,years,policy,stage="coarse") {
  edges <- branches <- list()
  for(key in names(groups)) for(year in years) {
    ids <- groups[[key]];axis <- strsplit(key,":",fixed=TRUE)[[1]][2]
    values <- vapply(records[ids],function(r) if(axis=="s") r$s else r$a,numeric(1))
    ids <- ids[order(values)];values <- sort(values);ids <- ids[!duplicated(values)]
    branch <- integer();next_id <- 0L
    for(i in seq_along(ids)) {
      r <- records[[ids[i]]];peaks <- r$annual[[as.character(year)]]
      local <- if(is.null(peaks)) NULL else peaks$local[peaks$local$retained,,drop=FALSE]
      current <- if(is.null(local)) integer() else setNames(rep(NA_integer_,nrow(local)),local$peak_id)
      if(i>1L) {
        e <- p8_edge(records[[ids[i-1L]]],r,axis,year,policy,stage);edges[[length(edges)+1L]] <- e
        t <- e$matching$table
        if(nrow(t) && !isTRUE(e$matching$ambiguous)) for(j in which(t$status=="MATCHED"))
          current[as.character(t$right_peak[j])] <- branch[as.character(t$left_peak[j])]
      }
      for(j in which(is.na(current))) {next_id <- next_id+1L;current[j] <- next_id}
      if(length(current)) for(j in seq_along(current)) branches[[length(branches)+1L]] <- data.frame(group=key,stage=stage,
        year=year,scenario=r$id,theta=if(axis=="s") r$s else r$a,peak_id=local$peak_id[j],
        branch=paste(key,year,current[j],sep="/"),week=local$week[j],height=local$height[j],prominence=local$prominence[j])
      branch <- current
    }
  }
  list(edges=edges,branches=p3_bind_rows(branches))
}
p8_run <- function(source,data,s=p8_defaults(),root=".") {
  baseline <- p8_source(source,data,s,root);check <- p8_settings(s,baseline)
  profiles <- list(fitted=p8_profile(baseline,s));if(s$annual$enabled) profiles$annual <- p8_profile(baseline,s,"annual")
  for(kind in names(profiles)) profiles[[kind]]$kind <- kind
  count <- length(profiles)*(1L+length(check$s_values)+length(check$a_values)+
    if(s$interventions$joint) length(check$s_values)*length(check$a_values) else 0L)
  if(count>s$interventions$max_scenarios) stop("Declared coarse grid exceeds scenario cap")
  records <- list();p <- baseline$parameters
  add <- function(profile,sval,a,axis,stage="coarse",level=0L,parent="") {
    id <- sprintf("scenario-%04d",length(records)+1L)
    r <- p8_evaluate(baseline,check,profiles[[profile]],sval,a,id,axis,stage,level,parent)
    records[[id]] <<- r;id
  }
  for(profile in names(profiles)) {
    add(profile,p$S0,1,"baseline")
    for(value in check$s_values) add(profile,value,1,"s")
    for(value in check$a_values) add(profile,p$S0,value,"a")
    if(s$interventions$joint) for(sval in check$s_values) for(a in check$a_values) add(profile,sval,a,"joint")
  }
  fitted <- Filter(function(r) r$profile=="fitted" && r$axis=="baseline",records)[[1]]
  reference <- baseline$evaluation$states[seq_len(check$horizon),,drop=FALSE]
  if(fitted$status!="VALID" || !p6r_close(as.numeric(as.matrix(fitted$trajectory[,names(reference)])),
    as.numeric(as.matrix(reference)),s$replay_tolerance)) stop("Baseline intervention invariance/replay failed")
  groups <- p8_series(records);coarse <- p8_track(records,groups,check$years,s$peaks)
  refined <- list();queue <- Filter(function(e) isTRUE(e$row$persistent_branch_exchange) && !isTRUE(e$matching$ambiguous),coarse$edges)
  added <- 0L;refinement_status <- if(s$refinement$enabled) "NO_DECLARED_SWITCH_CANDIDATES" else "NOT_REQUESTED"
  if(s$refinement$enabled) for(level in seq_len(s$refinement$levels)) {
    next_queue <- list()
    for(e in queue) {
      if(added>=s$refinement$max_points || length(records)>=s$interventions$max_scenarios) {refinement_status <- "POINT_CAP_REACHED";break}
      left <- records[[e$left]];right <- records[[e$right]]
      theta <- (e$row$theta_left+e$row$theta_right)/2
      if(theta %in% c(e$row$theta_left,e$row$theta_right)) next
      id <- add(left$profile,if(e$axis=="s") theta else left$s,if(e$axis=="a") theta else left$a,
        e$axis,"refined",level,paste(e$left,e$right,sep="/"));added <- added+1L
      refinement_status <- "MIDPOINT_DIAGNOSTICS_STORED"
      for(pair in list(c(e$left,id),c(id,e$right))) {
        child <- p8_edge(records[[pair[1]]],records[[pair[2]]],e$axis,e$year,s$peaks,"refined",level)
        refined[[length(refined)+1L]] <- child
        a <- records[[pair[1]]]$annual[[as.character(e$year)]];b <- records[[pair[2]]]$annual[[as.character(e$year)]]
        tied <- !is.null(a) && !is.null(b) && (length(a$global_ties)>1 || length(b$global_ties)>1) &&
          !isTRUE(child$matching$ambiguous) && nrow(child$matching$table)>=2
        if(child$row$mechanism=="PEAK_SWITCHING" || tied) next_queue[[length(next_queue)+1L]] <- child
      }
    }
    queue <- next_queue;if(!length(queue)) break
  }
  augmented <- if(added) p8_track(records,p8_series(records),check$years,s$peaks,"augmented_grid")$branches else NULL
  list(type="phase8_shift_mechanisms",status="COUNTERFACTUAL_SCENARIOS_STORED_NOT_CAUSAL_INFERENCE",config=source$config,
    diagnostic_settings=s,source_hash=source$artifact_hash,source_identity=source$identity,source_config=source$config,
    source_ancestry=baseline$source_ancestry,baseline_parameters=p,baseline_coordinates=baseline$coordinates,
    baseline_replay=baseline$replay,profiles=profiles,time_contract=check,records=records,
    coarse_edges=coarse$edges,refined_edges=refined,branches=coarse$branches,augmented_branches=augmented,
    refinement_status=refinement_status,source_warnings=baseline$warnings,
    outcome="Production StoI incidence; initialization row excluded; no stochastic observations or refitting")
}
p8_report <- function(result,out) {
  if(!identical(result$type,"phase8_shift_mechanisms")) stop("Saved Phase 8 result required")
  csv <- function(t,n) utils::write.csv(t,file.path(out,paste0(n,".csv")),row.names=FALSE)
  manifest <- weekly <- annual <- local <- matches <- list()
  for(r in result$records) {
    manifest[[length(manifest)+1L]] <- data.frame(scenario=r$id,profile=r$profile,axis=r$axis,s=r$s,a=r$a,
      stage=r$stage,level=r$level,parent=r$parent,status=r$status,reason=r$reason,
      role=if(r$profile=="annual") "HYPOTHETICAL_ANNUAL_COUNTERFACTUAL" else if(r$axis=="baseline") "FITTED_BASELINE_REPLAY" else "FITTED_SHAPE_COUNTERFACTUAL")
    if(r$status!="VALID") next
    t <- r$trajectory;t$scenario <- r$id;t$profile <- r$profile;t$s <- r$s;t$a <- r$a;weekly[[r$id]] <- t
    t <- r$annual_summary;t$scenario <- r$id;t$profile <- r$profile;t$axis <- r$axis;t$s <- r$s;t$a <- r$a;t$stage <- r$stage;annual[[r$id]] <- t
    for(year in names(r$annual)) {
      t <- r$annual[[year]]$local;t$scenario <- rep(r$id,nrow(t));t$year <- rep(as.integer(year),nrow(t));local[[paste(r$id,year)]] <- t
    }
  }
  edges <- c(result$coarse_edges,result$refined_edges)
  transitions <- p3_bind_rows(lapply(edges,`[[`,"row"))
  for(i in seq_along(edges)) {
    t <- edges[[i]]$matching$table
    if(!nrow(t)) next
    t$left <- edges[[i]]$left;t$right <- edges[[i]]$right;t$year <- edges[[i]]$year;t$stage <- edges[[i]]$row$stage
    matches[[i]] <- t
  }
  m <- p3_bind_rows(manifest);w <- p3_bind_rows(weekly);y <- p3_bind_rows(annual)
  csv(m,"scenario_manifest");csv(w,"weekly_trajectories");csv(y,"annual_peak_summary");csv(p3_bind_rows(local),"local_peaks")
  csv(result$branches,"peak_branches");csv(result$augmented_branches,"augmented_peak_branches")
  csv(transitions,"shift_transitions");csv(transitions,"mechanism_classification");csv(p3_bind_rows(matches),"peak_matching")
  annual_summary <- list()
  if("annual" %in% names(result$profiles)) for(year in result$time_contract$years) {
    t <- transitions[transitions$profile=="annual" & transitions$year==year,,drop=FALSE]
    annual_summary[[length(annual_summary)+1L]] <- data.frame(year=year,
      continuous_transmission_maxima=result$profiles$annual$annual_validation$continuous_maxima,
      discrete_transmission_maxima=result$profiles$annual$annual_validation$weekly_maxima,
      switching_candidates=sum(t$mechanism=="PEAK_SWITCHING"),displacement_candidates=sum(t$mechanism=="PEAK_DISPLACEMENT"),
      creation_annihilation_candidates=sum(t$mechanism=="PEAK_CREATION_ANNIHILATION"),
      status=if(any(t$mechanism %in% c("PEAK_SWITCHING","PEAK_DISPLACEMENT","PEAK_CREATION_ANNIHILATION") &
        t$global_circular_displacement>=result$diagnostic_settings$peaks$min_displacement,na.rm=TRUE))
        "CANDIDATE_COUNTEREXAMPLE_REQUIRES_VALIDATION" else "NO_COUNTEREXAMPLE_DETECTED_IN_DECLARED_GRID",
      necessity_conclusion="NOT_ESTABLISHED_BY_AUTOMATIC_CLASSIFICATION")
  }
  csv(p3_bind_rows(annual_summary),"annual_counterexample_summary")
  if("annual" %in% names(result$profiles)) csv(result$profiles$annual$annual_validation,"annual_forcing_validation")
  csv(data.frame(parameter=names(result$baseline_parameters),value=vapply(result$baseline_parameters,function(x) paste(x,collapse=";"),character(1))),"fixed_source_parameters")
  grDevices::pdf(file.path(out,"shift_mechanisms.pdf"),width=11,height=8)
  tryCatch({
    for(profile in names(result$profiles)) {
      base <- Filter(function(r) r$profile==profile && r$axis=="baseline" && r$status=="VALID",result$records)
      if(length(base)) {
        z <- base[[1]]$trajectory
        graphics::par(mfrow=c(2,1))
        graphics::plot(1:52,result$profiles[[profile]]$beta,type="l",xlab="Seasonal week",ylab="beta",main=paste(profile,result$profiles[[profile]]$label))
        graphics::plot(z$row,z$incidence,type="l",xlab="Saved-support row (initial row excluded from peaks)",ylab="Incident infections",main="Deterministic baseline/scenario; no observed data claim")
        graphics::par(mfrow=c(1,1))
      }
      for(year in result$time_contract$years) for(axis in c("s","a")) {
        t <- y[y$profile==profile & y$year==year & y$axis==axis & y$stage=="coarse" & y$status=="ANALYZED",,drop=FALSE]
        if(!nrow(t)) next
        t <- t[order(t[[axis]]),]
        graphics::plot(t[[axis]],t$global_week,type="b",ylim=c(1,52),xlab=axis,ylab="Global incidence peak week",
          main=paste(profile,year,"finite-grid timing; ambiguity in tables"))
        b <- result$branches
        b <- b[b$year==year & grepl(paste0("^",profile,":",axis,":"),b$group),,drop=FALSE]
        if(nrow(b)) {
          graphics::plot(range(b$theta),range(b$height),type="n",xlab=axis,ylab="Local incidence peak heights",main=paste(profile,year,"operational branches"))
          ids <- unique(b$branch)
          for(i in seq_along(ids)) {v <- b[b$branch==ids[i],];graphics::lines(v$theta,v$height,type="b",col=i)}
        }
      }
      t <- transitions[transitions$profile==profile,,drop=FALSE]
      if(nrow(t)) {
        left <- m[match(t$left,m$scenario),];right <- m[match(t$right,m$scenario),]
        labels <- unique(t$mechanism);colors <- match(t$mechanism,labels)
        graphics::plot((left$s+right$s)/2,(left$a+right$a)/2,col=colors,pch=19,xlab="s midpoint",ylab="a midpoint",
          main=paste(profile,"neighbor-edge classifications (no interpolation)"))
        graphics::legend("topright",labels,col=seq_along(labels),pch=19,bty="n",cex=.65)
      }
    }
    if(!length(result$profiles)) {graphics::plot.new();graphics::title("No simulations available")}
  },finally=grDevices::dev.off())
  writeLines(c("# Phase 8 SIRS shift mechanisms","",paste("Source Phase 7C hash:",result$source_hash),
    "Baseline is ONLY the accepted retained lambda=1 endpoint. Original endpoints/configs/certification remain unchanged.",
    paste("Replay:",result$baseline_replay),paste("Horizon rows:",result$time_contract$horizon,"| years:",paste(result$time_contract$years,collapse=",")),
    result$outcome,"",
    "All interventions are deterministic counterfactuals. I0 fraction is fixed; s redistributes S/R. a changes beta scale only. Annual forcing is hypothetical, not fitted.",
    "Initial row has no infection interval and is excluded. Incomplete years are excluded by default; each analyzed year is retained separately. No burn-in or attractor-convergence assumption is used.",
    "Peak matching uses circular distance AND height continuity. Boundary/plateau/tie/threshold/matching ambiguity remains unresolved. Weekly peak-count changes are not certified continuous-time bifurcations.",
    "A finite-grid displacement is not a mathematical discontinuity. Switching requires two persistent branches with reversed heights and exchanged global winner, not global week alone.",
    paste("Refinement:",result$refinement_status,"| original coarse records preserved"),
    "Annual counterexample labels are candidates only; failure to find one does not establish necessity. No cause of Japan's historical shift or identifiability is established.",
    "The exact map gives delta I = incidence - infectious removals. This is a prevalence growth balance, NOT an incidence peak condition. beta*S/N=gamma is not asserted as an exact criterion.",
    "Both actual state mass S+I+R and fixed population denominator N are exported; births/deaths do not guarantee later mass equals N.",
    "Saved-only report: no propagation, fit, peak detection or new refinement.","",
    "## Recorded thresholds","```",capture.output(str(result$diagnostic_settings)),"```",result$source_warnings),file.path(out,"phase8_summary.md"))
  invisible(result)
}
