# Phase 12 diagnostics use the accepted weekly transition without changing it.
p12_bind <- function(rows) {
  rows <- Filter(function(x) !is.null(x) && nrow(x)>0L,rows)
  if(!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows,names),use.names=FALSE))
  do.call(rbind,lapply(rows,function(x) {
    for(n in setdiff(columns,names(x))) x[[n]] <- rep(NA,nrow(x))
    x[,columns,drop=FALSE]
  }))
}
p12_defaults <- function() {
  h <- log(1.2)
  list(source_lambda=1,replay_tolerance=1e-7,years=2018:2019,
    peak=p10_defaults()$identity,peaks=p10_defaults()$peaks,
    a=list(start=h/4,halvings=4L,min_step=h/64,absolute_tolerance=.5,
      relative_tolerance=.1,near_zero_scale=.25,noise_tolerance=.05,
      max_scenarios=21L),
    b=list(v_grid=seq(-h,h,length.out=9L),u_grid=seq(-h,h,length.out=9L),
      levels=4L,max_refinement=48L,per_year=24L,max_scenarios=80L,
      rapid_weeks=2,near_height_fraction=.02),
    c=list(grid=seq(-h,h,length.out=5L),levels=1L,max_scenarios=90L,
      additivity_tolerance_weeks=.5,numerical_uncertainty_weeks=.25),
    d=list(calendar_start_week=1L,npieff=1,initial_s_multipliers=c(.9,1,1.1),
      initial_I_multipliers=c(.8,1,1.2),a_values=c(1/1.2,1,1.2),
      checkpoints=c(25L,50L,100L,200L),periods=1:4,max_years=200L,
      consecutive_cycles=3L,state_abs=1e-5,state_rel=1e-7,
      profile_abs=1e-5,profile_rel=1e-6,low_incidence=1e-6,
      max_scenarios=27L,store_full_trajectories=FALSE))
}
p12_validate <- function(s,module,s0,I0) {
  if(!identical(s$source_lambda,1)) stop("Phase 12 requires lambda=1")
  scalar_number(s$replay_tolerance,"replay_tolerance",0,strict_lower=TRUE)
  if(!identical(as.integer(s$years),2018:2019)) stop("Phase 12 requires complete 2018 and 2019")
  if(module=="a") {
    z <- s$a;integer_scalar(z$halvings,"halvings",1L)
    for(n in c("start","min_step","absolute_tolerance","relative_tolerance","near_zero_scale","noise_tolerance")) scalar_number(z[[n]],n,0,strict_lower=TRUE)
    if(z$start/2^z$halvings<z$min_step) stop("Nested sequence below minimum step")
    integer_scalar(z$max_scenarios,"max_scenarios",1L)
    if(1+4*(z$halvings+1L)>z$max_scenarios) stop("12A scenario budget too small")
  } else if(module=="b") {
    z <- s$b
    for(n in c("v_grid","u_grid")) {
      g <- z[[n]];finite_numeric(g,n)
      if(length(g)<3L || any(diff(g)<=0) || !any(g==0)) stop("Ordered axis grid containing zero required")
    }
    if(!any(abs(z$v_grid-log(1/1.2))<1e-14)) stop("12B must include a=1/1.2 exactly")
    for(n in c("levels","max_refinement","per_year","max_scenarios")) integer_scalar(z[[n]],n,if(n=="max_scenarios") 1L else 0L)
    if(length(z$v_grid)+length(z$u_grid)-1L+z$max_refinement>z$max_scenarios) stop("12B budget below declared allocation")
    for(n in c("rapid_weeks","near_height_fraction")) scalar_number(z[[n]],n,0)
  } else if(module=="c") {
    z <- s$c;g <- z$grid;finite_numeric(g,"joint grid")
    if(length(g)<3L || any(diff(g)<=0) || !any(g==0) || !isTRUE(all.equal(g,-rev(g),tolerance=1e-12))) stop("Symmetric joint grid required")
    for(n in c("levels","max_scenarios")) integer_scalar(z[[n]],n,if(n=="levels") 0L else 1L)
    if((2^z$levels*(length(g)-1L)+1L)^2>z$max_scenarios) stop("12C budget below nested grid")
    for(n in c("additivity_tolerance_weeks","numerical_uncertainty_weeks")) scalar_number(z[[n]],n,0,strict_lower=TRUE)
  } else if(module=="d") {
    z <- s$d;integer_scalar(z$calendar_start_week,"calendar_start_week",1L)
    if(z$calendar_start_week>52L) stop("Invalid seasonal phase")
    scalar_number(z$npieff,"npieff",0)
    for(n in c("initial_s_multipliers","initial_I_multipliers","a_values")) {finite_numeric(z[[n]],n);if(any(z[[n]]<=0)) stop("Positive initial and transmission multipliers required")}
    if(!any(z$a_values==1) || !any(abs(z$a_values-1/1.2)<1e-14)) stop("12D must cover baseline and 12B lower endpoint")
    for(n in c("checkpoints","periods")) {finite_numeric(z[[n]],n);if(any(z[[n]]!=as.integer(z[[n]]) | z[[n]]<=0)) stop("Positive integer years/periods required")}
    for(n in c("max_years","consecutive_cycles","max_scenarios")) integer_scalar(z[[n]],n,1L)
    if(max(z$checkpoints)>z$max_years || length(z$a_values)*length(z$initial_s_multipliers)*length(z$initial_I_multipliers)>z$max_scenarios) stop("12D budget/horizon exceeded")
    for(n in c("state_abs","state_rel","profile_abs","profile_rel","low_incidence")) scalar_number(z[[n]],n,0,strict_lower=TRUE)
    if(!is.logical(z$store_full_trajectories) || length(z$store_full_trajectories)!=1L) stop("Explicit storage policy required")
  } else stop("Unknown Phase 12 module")
  if(module!="d") {
    u <- if(module=="a") c(-s$a$start,s$a$start) else if(module=="b") s$b$u_grid else s$c$grid
    if(any(s0*exp(u)<=0 | s0*exp(u)>=1-I0)) stop("Susceptibility range violates simplex")
  }
  invisible(s)
}
p12_context <- function(source,data,s,root,module) {
  b <- p8_source(source,data,s,root);p12_validate(s,module,b$parameters$S0,b$parameters$I0)
  cfg <- p10_defaults();cfg$time$years <- s$years
  u <- if(module=="a") c(-s$a$start,s$a$start) else if(module=="b") s$b$u_grid else if(module=="c") s$c$grid else 0
  v <- if(module=="a") c(-s$a$start,s$a$start) else if(module=="b") s$b$v_grid else if(module=="c") s$c$grid else 0
  cfg$interventions$s_values <- sort(unique(c(b$parameters$S0,b$parameters$S0*exp(u))))
  cfg$interventions$a_values <- sort(unique(c(1,exp(v))))
  cfg$interventions$max_scenarios <- if(module=="a") s$a$max_scenarios else if(module=="b") s$b$max_scenarios else if(module=="c") s$c$max_scenarios else 1L
  check <- p8_settings(cfg,b)
  if(!all(check$complete)) stop("Complete reporting windows required")
  profile <- p8_profile(b,cfg);profile$kind <- "fitted"
  base <- p8_evaluate(b,check,profile,b$parameters$S0,1,"baseline","baseline")
  ref <- b$evaluation$states[seq_len(check$horizon),,drop=FALSE]
  if(base$status!="VALID" || !p6r_close(as.numeric(as.matrix(base$trajectory[,names(ref)])),as.numeric(as.matrix(ref)),s$replay_tolerance)) stop("Exact baseline replay failed")
  list(baseline=b,check=check,profile=profile,base=base)
}
p12_scenarios <- function(ctx,coordinates,max_scenarios) {
  coordinates <- unique(coordinates[,c("u","v"),drop=FALSE]);coordinates <- rbind(data.frame(u=0,v=0),coordinates)
  coordinates <- coordinates[!duplicated(p11_key(coordinates$u,coordinates$v)),,drop=FALSE]
  if(nrow(coordinates)>max_scenarios) stop("Scenario cap exceeded")
  coordinates$id <- sprintf("scenario-%04d",seq_len(nrow(coordinates)))
  records <- setNames(vector("list",nrow(coordinates)),coordinates$id)
  for(i in seq_len(nrow(coordinates))) {
    row <- coordinates[i,];records[[row$id]] <- if(row$u==0 && row$v==0) ctx$base else
      p8_evaluate(ctx$baseline,ctx$check,ctx$profile,ctx$baseline$parameters$S0*exp(row$u),exp(row$v),row$id,
        if(row$u==0) "beta" else if(row$v==0) "S" else "joint","phase12")
    records[[row$id]]$id <- row$id
    records[[row$id]]$u <- row$u;records[[row$id]]$v <- row$v
  }
  coordinates$s <- ctx$baseline$parameters$S0*exp(coordinates$u);coordinates$a <- exp(coordinates$v)
  coordinates$status <- vapply(records,`[[`,character(1),"status")
  list(plan=coordinates,records=records)
}
p12_observe <- function(records,s) {
  x <- p11_outcomes(records,s$years,s)
  list(peaks=x$peaks,windows=x$windows,annual=x$annual)
}
p12_peak <- function(peaks,window) {
  if(nrow(window)!=1L || window$status!="GENUINE_LOCAL_PEAK") return(NULL)
  p <- peaks[peaks$scenario==window$scenario & peaks$year==window$year & peaks$peak_id==window$peak_id,,drop=FALSE]
  if(nrow(p)!=1L || !isTRUE(p$genuine)) return(NULL)
  p
}
p12_pair <- function(L,R,d,idL,idR,year) {
  f <- p10_assignment(L,R,d);v <- p10_assignment(L,R,d,TRUE)
  ok <- isTRUE(p10_identity_consistent(f,v))
  for(t in list(f,v)) if(nrow(t)) {
    t$left <- idL;t$right <- idR;t$year <- year;t$direction_consistent <- ok
    if(identical(t$direction[1],"forward")) forward <- t else reverse <- t
  }
  list(forward=if(exists("forward")) forward else data.frame(),reverse=if(exists("reverse")) reverse else data.frame(),consistent=ok)
}
p12_match_winner <- function(peaks,wl,wc,wr,d) {
  if(any(vapply(list(wl,wc,wr),function(w) nrow(w)!=1L || w$status!="GENUINE_LOCAL_PEAK",logical(1)))) return(FALSE)
  left <- peaks[peaks$scenario==wl$scenario & peaks$year==wl$year,,drop=FALSE]
  center <- peaks[peaks$scenario==wc$scenario & peaks$year==wc$year,,drop=FALSE]
  right <- peaks[peaks$scenario==wr$scenario & peaks$year==wr$year,,drop=FALSE]
  for(q in list(list(left,center,wl$peak_id,wc$peak_id),list(center,right,wc$peak_id,wr$peak_id))) {
    f <- p10_assignment(q[[1]],q[[2]],d);v <- p10_assignment(q[[1]],q[[2]],d,TRUE)
    if(!p10_identity_consistent(f,v) || !any(f$status=="MATCHED" & f$left_peak==q[[3]] & f$right_peak==q[[4]])) return(FALSE)
  }
  TRUE
}
p12_a <- function(source,data,s=p12_defaults(),root=".") {
  ctx <- p12_context(source,data,s,root,"a")
  h <- s$a$start/2^(0:s$a$halvings)
  coordinates <- data.frame(u=c(h,-h,rep(0,2*length(h))),v=c(rep(0,2*length(h)),h,-h))
  sc <- p12_scenarios(ctx,coordinates,s$a$max_scenarios);obs <- p12_observe(sc$records,s)
  rows <- list();est <- c("quadratic","cubic")
  id <- function(u,v) sc$plan$id[match(p11_key(u,v),p11_key(sc$plan$u,sc$plan$v))]
  for(year in s$years) for(axis in c("u","v")) for(step in h) {
    ids <- c(id(if(axis=="u") -step else 0,if(axis=="v") -step else 0),id(0,0),
      id(if(axis=="u") step else 0,if(axis=="v") step else 0))
    ww <- lapply(ids,function(i) obs$windows[obs$windows$scenario==i & obs$windows$year==year,,drop=FALSE])
    pp <- lapply(ww,function(w) p12_peak(obs$peaks,w))
    valid <- p12_match_winner(obs$peaks,ww[[1]],ww[[2]],ww[[3]],s$peak)
    for(method in est) {
      field <- paste0(method,"_position");flag <- paste0(method,"_status")
      T <- vapply(pp,function(p) if(!is.null(p) && p[[flag]]==if(method=="cubic") "KNOT_MAXIMUM_AUDIT" else "OK") p[[field]] else NA_real_,numeric(1))
      good <- valid && all(is.finite(T))
      x <- c(-step,0,step)
      rows[[length(rows)+1L]] <- data.frame(year=year,axis=axis,h=step,method=method,
        minus=ids[1],center=ids[2],plus=ids[3],minus_coordinate=x[1],plus_coordinate=x[3],
        weekly_minus=if(nrow(ww[[1]])) ww[[1]]$week else NA_real_,weekly_center=if(nrow(ww[[2]])) ww[[2]]$week else NA_real_,weekly_plus=if(nrow(ww[[3]])) ww[[3]]$week else NA_real_,
        ties_minus=if(nrow(ww[[1]])) ww[[1]]$ties else "",ties_center=if(nrow(ww[[2]])) ww[[2]]$ties else "",ties_plus=if(nrow(ww[[3]])) ww[[3]]$ties else "",
        branch_minus=if(nrow(ww[[1]])) ww[[1]]$peak_id else NA_integer_,branch_center=if(nrow(ww[[2]])) ww[[2]]$peak_id else NA_integer_,branch_plus=if(nrow(ww[[3]])) ww[[3]]$peak_id else NA_integer_,
        support_minus=!is.null(pp[[1]]) && pp[[1]]$support,support_center=!is.null(pp[[2]]) && pp[[2]]$support,
        support_plus=!is.null(pp[[3]]) && pp[[3]]$support,
        curvature_minus=if(!is.null(pp[[1]])) pp[[1]]$curvature else NA_real_,
        curvature_center=if(!is.null(pp[[2]])) pp[[2]]$curvature else NA_real_,
        curvature_plus=if(!is.null(pp[[3]])) pp[[3]]$curvature else NA_real_,
        moving_weekly_stencil=if(all(vapply(ww,nrow,integer(1))==1L)) length(unique(vapply(ww,`[[`,numeric(1),"position")))>1L else NA,
        Tminus=T[1],T0=T[2],Tplus=T[3],sensitivity=if(good && method=="quadratic") p10_three_point(x,T,2L) else NA_real_,
        finite_minus_baseline_week=if(all(vapply(ww[1:2],nrow,integer(1))==1L)) ww[[1]]$week-ww[[2]]$week else NA_real_,
        finite_plus_baseline_week=if(all(vapply(ww[2:3],nrow,integer(1))==1L)) ww[[3]]$week-ww[[2]]$week else NA_real_,
        interpolated_minus_baseline=if(all(is.finite(T[1:2]))) T[1]-T[2] else NA_real_,
        interpolated_plus_baseline=if(all(is.finite(T[2:3]))) T[3]-T[2] else NA_real_,
        weekly_argmax_difference=if(all(vapply(ww,nrow,integer(1))==1L)) ww[[3]]$week-ww[[1]]$week else NA_real_,
        valid=good && method=="quadratic",reason=if(good && method=="cubic") "KNOT_AUDIT_ONLY" else
          if(good) "BRANCH_CONSISTENT" else "INVALID_BRANCH_STENCIL_OR_INTERPOLATION",stringsAsFactors=FALSE)
    }
  }
  tab <- do.call(rbind,rows);tab$change_from_previous <- NA_real_;tab$stability <- "NOT_ESTIMABLE"
  for(year in s$years) for(axis in c("u","v")) for(method in "quadratic") {
    ix <- which(tab$year==year & tab$axis==axis & tab$method==method);ix <- ix[order(-tab$h[ix])]
    if(length(ix)>1L) for(k in 2:length(ix)) {
      i <- ix[k];j <- ix[k-1L];d <- abs(tab$sensitivity[i]-tab$sensitivity[j]);tab$change_from_previous[i] <- d
      scale <- max(s$a$near_zero_scale,abs(tab$sensitivity[i]),abs(tab$sensitivity[j]))
      tab$stability[i] <- if(!tab$valid[i] || !tab$valid[j]) "INVALID_BRANCH_STENCIL" else if(d<=s$a$absolute_tolerance+s$a$relative_tolerance*scale) "STABLE_WITHIN_DECLARED_TOLERANCE" else if(d<=s$a$noise_tolerance) "UNRESOLVED_NUMERICAL_FLOOR" else "STEP_DEPENDENT"
    }
  }
  comparison <- list()
  for(year in s$years) for(step in h) for(method in "quadratic") {
    a <- tab[tab$year==year & tab$h==step & tab$method==method & tab$axis=="u",]
    b <- tab[tab$year==year & tab$h==step & tab$method==method & tab$axis=="v",]
    good <- nrow(a)==1L && nrow(b)==1L && a$valid && b$valid
    comparison[[length(comparison)+1L]] <- data.frame(year=year,h=step,method=method,E_S=if(good) a$sensitivity else NA_real_,
      E_beta=if(good) b$sensitivity else NA_real_,abs_E_S=if(good) abs(a$sensitivity) else NA_real_,abs_E_beta=if(good) abs(b$sensitivity) else NA_real_,
      magnitude_difference=if(good) abs(b$sensitivity)-abs(a$sensitivity) else NA_real_,
      beta_larger=if(good) abs(b$sensitivity)>abs(a$sensitivity) else NA,valid=good)
  }
  comp <- do.call(rbind,comparison)
  stability <- list()
  for(year in s$years) {
    q <- comp[comp$year==year,,drop=FALSE]
    d <- tab[tab$year==year & tab$method=="quadratic",,drop=FALSE]
    status <- if(!any(q$valid)) "NOT_ESTIMABLE" else if(any(!q$valid)) "INVALID_BRANCH_STENCIL" else
      if(any(d$stability=="STEP_DEPENDENT")) "STEP_DEPENDENT" else
      if(any(q$magnitude_difference<=s$a$near_zero_scale)) "UNRESOLVED_NUMERICAL_FLOOR" else
      if(length(unique(q$beta_larger))>1L) "STEP_DEPENDENT" else
      if(all(d$stability[d$stability!="NOT_ESTIMABLE"]=="STABLE_WITHIN_DECLARED_TOLERANCE"))
        "STABLE_WITHIN_DECLARED_TOLERANCE" else "NOT_ESTIMABLE"
    stability[[length(stability)+1L]] <- data.frame(year=year,status=status,
      beta_larger_across_valid_steps=if(any(q$valid)) all(q$beta_larger[q$valid]) else NA,
      supported_estimator="quadratic",valid_steps=sum(q$valid))
  }
  list(type="phase12_a",status="DIAGNOSTICS_STORED",settings=s,source_hash=source$artifact_hash,
    source_ancestry=ctx$baseline$source_ancestry,baseline_replay=ctx$baseline$replay,
    plan=sc$plan,records=sc$records,peaks=obs$peaks,windows=obs$windows,annual=obs$annual,
    sensitivities=tab,comparisons=comp,stability=do.call(rbind,stability))
}
p12_b_edges <- function(plan,obs,axis,s) {
  q <- plan[if(axis=="v") plan$u==0 else plan$v==0,,drop=FALSE]
  q <- q[order(q[[axis]]),,drop=FALSE];rows <- list();matches <- list();heights <- list()
  for(year in s$years) for(k in seq_len(nrow(q)-1L)) {
    l <- q[k,];r <- q[k+1L,];L <- obs$peaks[obs$peaks$scenario==l$id & obs$peaks$year==year,,drop=FALSE]
    R <- obs$peaks[obs$peaks$scenario==r$id & obs$peaks$year==year,,drop=FALSE]
    if(anyNA(c(L$genuine,R$genuine,L$retained,R$retained)) ||
       any(!is.finite(c(L$position,R$position,L$height,R$height))))
      stop("Invalid peak evidence in Phase 12B")
    wl <- obs$windows[obs$windows$scenario==l$id & obs$windows$year==year,,drop=FALSE]
    wr <- obs$windows[obs$windows$scenario==r$id & obs$windows$year==year,,drop=FALSE]
    if(nrow(wl)>1L || nrow(wr)>1L) stop("Duplicate annual window row in Phase 12B")
    pair <- p12_pair(L,R,s$peak,l$id,r$id,year);f <- pair$forward
    for(m in list(pair$forward,pair$reverse)) if(nrow(m)) {m$axis <- axis;matches[[length(matches)+1L]] <- m}
    good <- f[which(f$status=="MATCHED"),,drop=FALSE]
    if(!pair$consistent) good <- good[FALSE,,drop=FALSE]
    same <- nrow(wl)==1L && nrow(wr)==1L &&
      isTRUE(any(good$left_peak==wl$peak_id & good$right_peak==wr$peak_id,na.rm=TRUE))
    switch <- FALSE;near <- FALSE
    if(nrow(good)>=2L) for(i in seq_len(nrow(good)-1L)) for(j in (i+1L):nrow(good)) {
      a <- L[match(good$left_peak[c(i,j)],L$peak_id),,drop=FALSE]
      b <- R[match(good$right_peak[c(i,j)],R$peak_id),,drop=FALSE]
      dl <- a$height[1]-a$height[2];dr <- b$height[1]-b$height[2]
      tol <- s$peaks$tie_abs+s$peaks$tie_rel*max(c(a$height,b$height))
      reversal <- isTRUE(dl*dr<0 && abs(dl)>tol && abs(dr)>tol)
      exchange <- reversal && !same && nrow(wl)==1L && nrow(wr)==1L &&
        isTRUE(wl$status=="GENUINE_LOCAL_PEAK") && isTRUE(wr$status=="GENUINE_LOCAL_PEAK") &&
        isTRUE(wl$peak_id %in% a$peak_id) && isTRUE(wr$peak_id %in% b$peak_id) &&
        isTRUE(a$peak_id[which.max(a$height)]==wl$peak_id) && isTRUE(b$peak_id[which.max(b$height)]==wr$peak_id)
      switch <- switch || exchange
      near <- near || isTRUE(min(abs(dl),abs(dr))/max(c(a$height,b$height))<=s$b$near_height_fraction)
      heights[[length(heights)+1L]] <- data.frame(year=year,axis=axis,left=l$id,right=r$id,
        left_peak1=a$peak_id[1],left_peak2=a$peak_id[2],right_peak1=b$peak_id[1],right_peak2=b$peak_id[2],
        left_height1=a$height[1],left_height2=a$height[2],right_height1=b$height[1],right_height2=b$height[2],
        weekly_order_reversal=reversal,global_exchange=exchange)
    }
    missing_status <- (nrow(wl)==1L && is.na(wl$status)) || (nrow(wr)==1L && is.na(wr$status))
    unresolved <- !pair$consistent || any(f$status=="AMBIGUOUS",na.rm=TRUE) ||
      nrow(wl)!=1L || nrow(wr)!=1L ||
      !isTRUE(wl$status=="GENUINE_LOCAL_PEAK") || !isTRUE(wr$status=="GENUINE_LOCAL_PEAK")
    window <- p12_b_window_support(obs,l$id,r$id,year,wl,wr,s$peak)
    movement <- if(nrow(wl)==1L && nrow(wr)==1L) wr$position-wl$position else NA_real_
    if(length(movement)!=1L || !is.numeric(movement)) stop("Invalid annual window position in Phase 12B")
    if(!is.finite(movement)) movement <- NA_real_
    movement_known <- is.finite(movement)
    raw_change <- sum(L$genuine)!=sum(R$genuine)
    filter_change <- sum(L$retained)!=sum(R$retained)
    classification <- if(switch && !unresolved) "SWITCH" else if(window && pair$consistent) "WINDOW_SELECTION" else
      if(same && !unresolved && movement_known && movement!=0) "DISPLACEMENT" else "UNRESOLVED"
    reason <- if(classification!="UNRESOLVED") NA_character_ else if(missing_status) "MISSING_PEAK_STATUS" else
      if(!movement_known) "MISSING_PEAK_MOVEMENT" else if(any(f$status=="AMBIGUOUS",na.rm=TRUE)) "AMBIGUOUS_BRANCH_MATCH" else
      "IDENTITY_OR_WINDOW_UNRESOLVED"
    priority <- if(switch) 1L else if(near) 2L else if(same && movement_known && abs(movement)>=s$b$rapid_weeks) 3L else
      if(window || unresolved) 4L else if(raw_change) 5L else 99L
    rows[[length(rows)+1L]] <- data.frame(year=year,axis=axis,left=l$id,right=r$id,
      coordinate_left=l[[axis]],coordinate_right=r[[axis]],a_left=exp(l$v),a_right=exp(r$v),
      signed_annual_week_difference=if(nrow(wl)==1L && nrow(wr)==1L) wr$week-wl$week else NA_real_,
      chronological_difference=movement,window_left=if(nrow(wl)) wl$status else "MISSING",window_right=if(nrow(wr)) wr$status else "MISSING",
      winner_left=if(nrow(wl)) wl$peak_id else NA_integer_,winner_right=if(nrow(wr)) wr$peak_id else NA_integer_,
      same_branch=same,identity_consistent=pair$consistent,raw_peak_count_left=sum(L$genuine),raw_peak_count_right=sum(R$genuine),
      retained_left=sum(L$retained),retained_right=sum(R$retained),raw_topology_change=raw_change,
      filtering_effect=filter_change && !raw_change,near_height=near,priority=priority,
      classification=classification,classification_reason=reason,stringsAsFactors=FALSE)
  }
  bind <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
  list(edges=bind(rows),matches=bind(matches),heights=bind(heights))
}
p12_b_context <- function(records,s) {
  rows <- list()
  for(r in records) if(r$status=="VALID") for(year in unique(r$trajectory$year)) {
    p <- p10_detect(r,year,s$peaks,s$peak)
    if(nrow(p)) rows[[length(rows)+1L]] <- p
  }
  p12_bind(rows)
}
p12_b_window_support <- function(obs,left,right,year,wl,wr,d) {
  if(nrow(wl)!=1L || nrow(wr)!=1L || is.null(obs$context_peaks) || !nrow(obs$context_peaks)) return(FALSE)
  if(!any(c(wl$status,wr$status)=="WINDOW_SELECTION",na.rm=TRUE)) return(FALSE)
  L <- obs$context_peaks[obs$context_peaks$scenario==left & abs(obs$context_peaks$year-year)<=1L,,drop=FALSE]
  R <- obs$context_peaks[obs$context_peaks$scenario==right & abs(obs$context_peaks$year-year)<=1L,,drop=FALSE]
  f <- p10_assignment(L,R,d);rev <- p10_assignment(L,R,d,TRUE)
  if(!isTRUE(p10_identity_consistent(f,rev))) return(FALSE)
  good <- f[which(f$status=="MATCHED"),,drop=FALSE]
  if(isTRUE(wl$status=="WINDOW_SELECTION") && isTRUE(wr$status=="GENUINE_LOCAL_PEAK")) {
    outside <- L$peak_id[which(L$year!=year & L$genuine & abs(L$position-wl$position)<=d$match_distance)]
    if(any(good$left_peak %in% outside & good$right_peak==wr$peak_id,na.rm=TRUE)) return(TRUE)
  }
  if(isTRUE(wr$status=="WINDOW_SELECTION") && isTRUE(wl$status=="GENUINE_LOCAL_PEAK")) {
    outside <- R$peak_id[which(R$year!=year & R$genuine & abs(R$position-wr$position)<=d$match_distance)]
    if(any(good$right_peak %in% outside & good$left_peak==wl$peak_id,na.rm=TRUE)) return(TRUE)
  }
  FALSE
}
p12_b_branches <- function(plan,peaks,windows,matches,years) {
  rows <- list()
  for(axis in c("u","v")) for(year in years) {
    q <- plan[if(axis=="v") plan$u==0 else plan$v==0,,drop=FALSE]
    q <- q[order(q[[axis]]),,drop=FALSE];prior <- character();next_id <- 0L
    for(k in seq_len(nrow(q))) {
      p <- peaks[peaks$scenario==q$id[k] & peaks$year==year & peaks$genuine,,drop=FALSE]
      current <- setNames(rep(NA_character_,nrow(p)),as.character(p$peak_id))
      if(k>1L && nrow(matches)) {
        m <- matches[matches$axis==axis & matches$year==year & matches$left==q$id[k-1L] &
          matches$right==q$id[k] & matches$direction=="forward" & matches$status=="MATCHED" &
          matches$direction_consistent,,drop=FALSE]
        for(j in seq_len(nrow(m))) current[as.character(m$right_peak[j])] <- prior[as.character(m$left_peak[j])]
      }
      for(j in which(is.na(current))) {next_id <- next_id+1L;current[j] <- paste(axis,year,next_id,sep="/")}
      if(nrow(p)) {
        p$axis <- axis;p$coordinate <- q[[axis]][k]
        p$branch <- unname(current[as.character(p$peak_id)])
        p$winner <- vapply(seq_len(nrow(p)),function(i) any(windows$scenario==p$scenario[i] &
          windows$year==year & windows$peak_id==p$peak_id[i] & windows$status=="GENUINE_LOCAL_PEAK"),logical(1))
        rows[[length(rows)+1L]] <- p
      }
      prior <- current
    }
  }
  out <- p12_bind(rows)
  if(nrow(out)) {
    out$accumulated_movement <- NA_real_
    for(branch in unique(out$branch)) {
      ix <- which(out$branch==branch);ix <- ix[order(out$coordinate[ix])]
      out$accumulated_movement[ix] <- out$position[ix]-out$position[ix[1]]
    }
  }
  out
}
p12_b_chronological <- function(plan,peaks,s) {
  matches <- branches <- list()
  for(axis in c("u","v")) {
    q <- plan[if(axis=="v") plan$u==0 else plan$v==0,,drop=FALSE]
    q <- q[order(q[[axis]]),,drop=FALSE]
    prior <- character();next_id <- 0L
    for(k in seq_len(nrow(q))) {
      p <- peaks[peaks$scenario==q$id[k] & peaks$genuine,,drop=FALSE]
      current <- setNames(rep(NA_character_,nrow(p)),as.character(p$peak_id))
      if(k>1L) {
        L <- peaks[peaks$scenario==q$id[k-1L],,drop=FALSE]
        R <- peaks[peaks$scenario==q$id[k],,drop=FALSE]
        pair <- p12_pair(L,R,s$peak,q$id[k-1L],q$id[k],NA_integer_)
        for(m in list(pair$forward,pair$reverse)) if(nrow(m)) {
          m$axis <- axis;m$coordinate_left <- q[[axis]][k-1L];m$coordinate_right <- q[[axis]][k]
          matches[[length(matches)+1L]] <- m
        }
        m <- pair$forward
        if(nrow(m) && pair$consistent) {
          m <- m[m$status=="MATCHED",,drop=FALSE]
          for(j in seq_len(nrow(m))) current[as.character(m$right_peak[j])] <- prior[as.character(m$left_peak[j])]
        }
      }
      for(j in which(is.na(current))) {next_id <- next_id+1L;current[j] <- paste(axis,"chronological",next_id,sep="/")}
      if(nrow(p)) {
        p$axis <- axis;p$coordinate <- q[[axis]][k]
        p$branch <- unname(current[as.character(p$peak_id)])
        branches[[length(branches)+1L]] <- p
      }
      prior <- current
    }
  }
  out <- p12_bind(branches)
  if(nrow(out)) {
    out$accumulated_movement <- NA_real_
    for(branch in unique(out$branch)) {
      ix <- which(out$branch==branch);ix <- ix[order(out$coordinate[ix])]
      out$accumulated_movement[ix] <- out$position[ix]-out$position[ix[1]]
    }
  }
  list(matches=p12_bind(matches),branches=out)
}
p12_b <- function(source,data,s=p12_defaults(),root=".") {
  ctx <- p12_context(source,data,s,root,"b")
  coords <- unique(rbind(data.frame(u=0,v=s$b$v_grid),data.frame(u=s$b$u_grid,v=0)))
  sc <- p12_scenarios(ctx,coords,s$b$max_scenarios);history <- decisions <- list();counts <- setNames(integer(length(s$years)),as.character(s$years));added <- 0L
  for(level in 0:s$b$levels) {
    obs <- p12_observe(sc$records,s);obs$context_peaks <- p12_b_context(sc$records,s)
    V <- p12_b_edges(sc$plan,obs,"v",s);U <- p12_b_edges(sc$plan,obs,"u",s)
    edges <- rbind(V$edges,U$edges);edges$level <- level;history[[length(history)+1L]] <- edges
    if(level==s$b$levels) break
    candidates <- edges[edges$priority<99L,,drop=FALSE]
    candidates <- candidates[order(candidates$priority,candidates$year,candidates$axis,candidates$coordinate_left),,drop=FALSE]
    new <- 0L
    for(i in seq_len(nrow(candidates))) {
      e <- candidates[i,];mid <- (e$coordinate_left+e$coordinate_right)/2
      u <- if(e$axis=="u") mid else 0;v <- if(e$axis=="v") mid else 0
      present <- any(abs(sc$plan$u-u)<1e-14 & abs(sc$plan$v-v)<1e-14)
      decision <- if(present) "SHARED_MIDPOINT" else if(counts[as.character(e$year)]>=s$b$per_year) "YEAR_CAP" else
        if(added>=s$b$max_refinement || nrow(sc$plan)>=s$b$max_scenarios) "SCENARIO_CAP" else "SELECTED"
      decisions[[length(decisions)+1L]] <- data.frame(year=e$year,axis=e$axis,level=level,coordinate_left=e$coordinate_left,
        coordinate_right=e$coordinate_right,priority=e$priority,midpoint=mid,decision=decision)
      if(decision!="SELECTED") next
      id <- sprintf("scenario-%04d",nrow(sc$plan)+1L)
      record <- p8_evaluate(ctx$baseline,ctx$check,ctx$profile,ctx$baseline$parameters$S0*exp(u),exp(v),id,
        if(e$axis=="u") "S" else "beta","refined",level+1L,paste(e$left,e$right,sep="/"))
      record$u <- u;record$v <- v;sc$records[[id]] <- record
      sc$plan <- rbind(sc$plan,data.frame(u=u,v=v,id=id,s=ctx$baseline$parameters$S0*exp(u),a=exp(v),status=record$status))
      counts[as.character(e$year)] <- counts[as.character(e$year)]+1L;added <- added+1L;new <- new+1L
    }
    if(!new) break
  }
  obs <- p12_observe(sc$records,s);obs$context_peaks <- p12_b_context(sc$records,s)
  V <- p12_b_edges(sc$plan,obs,"v",s);U <- p12_b_edges(sc$plan,obs,"u",s)
  terminal <- rbind(V$edges,U$edges);unresolved <- terminal[terminal$classification=="UNRESOLVED" | terminal$priority<99L,,drop=FALSE]
  if(nrow(unresolved)) unresolved$reason <- ifelse(unresolved$priority<99L,"REFINEMENT_LIMIT_OR_CAP","IDENTITY_OR_WINDOW_UNRESOLVED")
  all_matches <- p12_bind(list(V$matches,U$matches))
  branches <- p12_b_branches(sc$plan,obs$peaks,obs$windows,all_matches,s$years)
  chronological <- p12_b_chronological(sc$plan,obs$context_peaks,s)
  list(type="phase12_b",status="DIAGNOSTICS_STORED",settings=s,source_hash=source$artifact_hash,
    source_ancestry=ctx$baseline$source_ancestry,baseline_replay=ctx$baseline$replay,plan=sc$plan,
    records=sc$records,peaks=obs$peaks,context_peaks=obs$context_peaks,windows=obs$windows,annual=obs$annual,terminal=terminal,
    matches=all_matches,branches=branches,chronological_matches=chronological$matches,
    chronological_branches=chronological$branches,competing_heights=p12_bind(list(V$heights,U$heights)),
    refinement_history=do.call(rbind,history),refinement_decisions=p12_bind(decisions),unresolved=unresolved)
}
p12_c <- function(source,data,s=p12_defaults(),root=".") {
  ctx <- p12_context(source,data,s,root,"c")
  g <- s$c$grid
  for(level in seq_len(s$c$levels)) g <- sort(unique(c(g,head(g,-1L)+diff(g)/2)))
  plan <- expand.grid(u=g,v=g);plan$stage <- "joint_grid"
  sc <- p12_scenarios(ctx,plan,s$c$max_scenarios);sc$plan$stage <- "joint_grid"
  obs <- p12_observe(sc$records,s)
  if(nrow(obs$peaks)) {
    edges <- p11_edges(sc$plan,obs$peaks,s$years,s$peak)
    loops <- if(nrow(edges)) p11_loops(sc$plan,edges,obs$peaks,s$years) else data.frame()
    regions <- p11_regions(sc$plan,obs$peaks,obs$windows,edges,data.frame(),loops)
    interaction <- p11_interactions(sc$plan,obs$windows,obs$peaks,edges,loops,s$years,s$peak)
  } else {
    edges <- loops <- data.frame();regions <- list(branches=data.frame(),winners=data.frame())
    interaction <- list(branch=data.frame(),global=data.frame())
  }
  x <- interaction$branch
  if(nrow(x)) {
    if(!nrow(loops)) {x$valid <- FALSE;x$delta_interaction <- NA_real_;x$reason <- "GRID_LOOP_IDENTITY_UNAVAILABLE"}
    x$interpolation_spread_weeks <- NA_real_
    for(i in which(x$valid)) {
      ids <- sc$plan$id[match(p11_key(c(0,x$u[i],0,x$u[i]),c(0,0,x$v[i],x$v[i])),
        p11_key(sc$plan$u,sc$plan$v))]
      bids <- c(x$branch00[i],x$branch_u0[i],x$branch_0v[i],x$branch_uv[i])
      p <- lapply(seq_along(ids),function(k) obs$peaks[obs$peaks$scenario==ids[k] &
        obs$peaks$year==x$year[i] & obs$peaks$peak_id==bids[k],,drop=FALSE])
      if(all(vapply(p,nrow,integer(1))==1L)) x$interpolation_spread_weeks[i] <-
        sum(vapply(p,function(q) abs(q$quadratic_position-q$position),numeric(1)))
    }
    x$numerical_uncertainty_weeks <- pmax(s$c$numerical_uncertainty_weeks,x$interpolation_spread_weeks,na.rm=TRUE)
    x$additivity <- ifelse(!x$valid,"NOT_ESTIMABLE",ifelse(abs(x$delta_interaction)<=s$c$additivity_tolerance_weeks+
      x$numerical_uncertainty_weeks,"APPROXIMATELY_ADDITIVE_WITHIN_TOLERANCE","NONADDITIVE_NUMERICALLY"))
  }
  list(type="phase12_c",status="DIAGNOSTICS_STORED",settings=s,source_hash=source$artifact_hash,
    source_ancestry=ctx$baseline$source_ancestry,baseline_replay=ctx$baseline$replay,plan=sc$plan,
    records=sc$records,peaks=obs$peaks,windows=obs$windows,annual=obs$annual,edges=edges,
    loops=loops,branches=regions$branches,winners=regions$winners,interaction=x,
    global_contrasts=interaction$global,unresolved=if(nrow(loops)) loops[loops$status!="CONSISTENT",,drop=FALSE] else data.frame())
}
p12_d_simulate <- function(beta,fixed,initial,a,z) {
  n <- 1L+52L*z$max_years
  week <- ((z$calendar_start_week-1L+seq_len(n)-1L)%%52L)+1L
  sim <- simulate_weekly_sirs(a*beta,initial,fixed,week,rep(z$npieff,n))
  sim$model_year <- c(0L,rep(seq_len(z$max_years),each=52L))
  sim$initialization_only <- seq_len(n)==1L
  sim
}
p12_d_checkpoint <- function(sim,year,z) {
  rows <- list();profiles <- split(sim$incidence[-1L],rep(seq_len(z$max_years),each=52L))
  states <- as.matrix(sim[1L+52L*seq_len(z$max_years),c("S","I","R")])
  for(period in z$periods) {
    ids <- seq.int(max(1L,year-z$consecutive_cycles+1L),year)
    if(min(ids)-period<1L) next
    state_error <- profile_error <- numeric(length(ids));state_ok <- profile_ok <- logical(length(ids))
    for(k in seq_along(ids)) {
      y <- ids[k];prev <- y-period
      a <- states[y,];b <- states[prev,]
      state_error[k] <- max(abs(a-b));state_ok[k] <- all(abs(a-b)<=z$state_abs+z$state_rel*pmax(abs(a),abs(b)))
      a <- profiles[[y]];b <- profiles[[prev]]
      profile_error[k] <- max(abs(a-b));profile_ok[k] <- all(abs(a-b)<=z$profile_abs+z$profile_rel*pmax(abs(a),abs(b)))
    }
    rows[[length(rows)+1L]] <- data.frame(checkpoint=year,period=period,state_error=max(state_error),
      profile_error=max(profile_error),state_matches=all(state_ok),profile_matches=all(profile_ok),
      recurrence_candidate=all(state_ok)&all(profile_ok))
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p12_d <- function(source,data,s=p12_defaults(),root=".") {
  ctx <- p12_context(source,data,s,root,"d");z <- s$d;p <- ctx$baseline$parameters
  fixed <- list(pop=ctx$baseline$data$fixed_parameters$pop,mu=p$mu,gamma=p$gamma,omega=p$omega,delta=p$delta,rho=p$rho)
  grid <- expand.grid(a=z$a_values,s_multiplier=z$initial_s_multipliers,I_multiplier=z$initial_I_multipliers)
  grid$id <- sprintf("scenario-%04d",seq_len(nrow(grid)))
  grid$s <- p$S0*grid$s_multiplier;grid$I0 <- p$I0*grid$I_multiplier
  if(any(grid$s<=0 | grid$I0<=0 | grid$s+grid$I0>=1)) stop("12D initial-condition grid violates simplex")
  grid$family <- ifelse(grid$a==1,"A_INITIAL_CONDITIONS","B_PERMANENT_TRANSMISSION")
  grid$regime <- "HYPOTHETICAL_REPEATING_FITTED_BETA_CONSTANT_NPIEFF"
  grid$status <- "PLANNED"
  checkpoints <- summaries <- profiles <- trajectories <- list()
  for(i in seq_len(nrow(grid))) {
    row <- grid[i,];sim <- tryCatch(p12_d_simulate(ctx$baseline$beta,fixed,make_initial_state(row$s,row$I0,fixed$pop),row$a,z),error=function(e) e)
    if(inherits(sim,"error")) {grid$status[i] <- paste("FAILED",conditionMessage(sim));next}
    grid$status[i] <- "VALID"
    if(z$store_full_trajectories) trajectories[[row$id]] <- sim
    for(year in z$checkpoints) {
      ix <- (1L+52L*(year-1L)+1L):(1L+52L*year)
      wave <- sim$incidence[ix];ties <- which(max(wave)-wave<=s$peaks$tie_abs+s$peaks$tie_rel*max(wave))
      q <- p12_d_checkpoint(sim,year,z)
      if(nrow(q)) {q$scenario <- row$id;checkpoints[[length(checkpoints)+1L]] <- q}
      summaries[[length(summaries)+1L]] <- data.frame(scenario=row$id,checkpoint=year,
        S=sim$S[tail(ix,1)],I=sim$I[tail(ix,1)],R=sim$R[tail(ix,1)],
        state_mass=sum(sim[tail(ix,1),c("S","I","R")]),peak_week=ties[1],near_ties=paste(ties,collapse=";"),
        peak_height=max(wave),cumulative_incidence=sum(wave),low_incidence=max(wave)<=z$low_incidence)
      for(phase_year in seq.int(max(1L,year-max(z$periods)+1L),year)) {
        phase_ix <- (1L+52L*(phase_year-1L)+1L):(1L+52L*phase_year)
        profiles[[length(profiles)+1L]] <- data.frame(scenario=row$id,checkpoint=year,
          phase_year=phase_year,week=1:52,incidence=sim$incidence[phase_ix])
      }
    }
  }
  cp <- p12_bind(checkpoints);su <- p12_bind(summaries);pr <- p12_bind(profiles)
  if(nrow(su)) {
    su$recurrence_status <- "UNRESOLVED_LONG_RUN"
    for(i in seq_len(nrow(su))) {
      q <- cp[cp$scenario==su$scenario[i] & cp$checkpoint==su$checkpoint[i] & cp$recurrence_candidate,,drop=FALSE]
      su$recurrence_status[i] <- if(su$low_incidence[i]) "LOW_INCIDENCE_CANDIDATE" else if(nrow(q))
        if(min(q$period)==1L) "ANNUAL_PERIODIC_CANDIDATE" else "MULTIYEAR_PERIODIC_CANDIDATE" else "UNRESOLVED_LONG_RUN"
    }
  }
  comparisons <- list()
  for(a in z$a_values) for(year in z$checkpoints) {
    q <- grid[grid$a==a & grid$status=="VALID",,drop=FALSE]
    if(nrow(q)<2L) next
    ref <- q$id[1];r <- su[su$scenario==ref & su$checkpoint==year,,drop=FALSE]
    rp <- pr$incidence[pr$scenario==ref & pr$checkpoint==year & pr$phase_year==year]
    for(id in q$id[-1L]) {
      t <- su[su$scenario==id & su$checkpoint==year,,drop=FALSE]
      tp <- pr$incidence[pr$scenario==id & pr$checkpoint==year & pr$phase_year==year]
      sd <- max(abs(as.numeric(t[,c("S","I","R")])-as.numeric(r[,c("S","I","R")])));pd <- max(abs(tp-rp))
      close <- all(abs(as.numeric(t[,c("S","I","R")])-as.numeric(r[,c("S","I","R")]))<=
        z$state_abs+z$state_rel*pmax(as.numeric(t[,c("S","I","R")]),as.numeric(r[,c("S","I","R")]))) &&
        all(abs(tp-rp)<=z$profile_abs+z$profile_rel*pmax(tp,rp))
      comparisons[[length(comparisons)+1L]] <- data.frame(a=a,checkpoint=year,reference=ref,scenario=id,
        max_state_difference=sd,max_profile_difference=pd,timing_difference_week=t$peak_week-r$peak_week,
        numerically_indistinguishable=close,status=if(close) "DECAYED_WITHIN_TESTED_TOLERANCE" else "PERSISTING_AT_TESTED_HORIZON")
    }
  }
  transmission <- list()
  dependence <- list()
  cmp <- p12_bind(comparisons)
  for(a in z$a_values) for(year in z$checkpoints) {
    q <- grid[grid$a==a,,drop=FALSE]
    states <- su[su$scenario %in% q$id & su$checkpoint==year,,drop=FALSE]
    c0 <- cmp[cmp$a==a & cmp$checkpoint==year,,drop=FALSE]
    all_recurrent <- nrow(states)==nrow(q) && all(states$recurrence_status %in% c("ANNUAL_PERIODIC_CANDIDATE","MULTIYEAR_PERIODIC_CANDIDATE","LOW_INCIDENCE_CANDIDATE"))
    dependence[[length(dependence)+1L]] <- data.frame(a=a,checkpoint=year,
      status=if(nrow(q)<2L || !all_recurrent || nrow(c0)!=nrow(q)-1L) "UNRESOLVED_LONG_RUN" else
        if(all(c0$numerically_indistinguishable)) "INITIAL_CONDITIONS_DECAYED_WITHIN_TOLERANCE" else "INITIAL_CONDITION_DEPENDENT_CANDIDATE")
  }
  reference <- grid$id[grid$a==1 & grid$s_multiplier==1 & grid$I_multiplier==1]
  if(length(reference)==1L) for(a in z$a_values[z$a!=1]) {
    candidate <- grid$id[grid$a==a & grid$s_multiplier==1 & grid$I_multiplier==1]
    if(length(candidate)!=1L) next
    for(year in z$checkpoints) {
      x <- su[su$scenario==reference & su$checkpoint==year,,drop=FALSE]
      y <- su[su$scenario==candidate & su$checkpoint==year,,drop=FALSE]
      periodic <- nrow(x)==1L && nrow(y)==1L && x$recurrence_status %in% c("ANNUAL_PERIODIC_CANDIDATE","MULTIYEAR_PERIODIC_CANDIDATE") &&
        y$recurrence_status %in% c("ANNUAL_PERIODIC_CANDIDATE","MULTIYEAR_PERIODIC_CANDIDATE")
      transmission[[length(transmission)+1L]] <- data.frame(a=a,checkpoint=year,baseline=reference,scenario=candidate,
        baseline_status=if(nrow(x)) x$recurrence_status else "MISSING",scenario_status=if(nrow(y)) y$recurrence_status else "MISSING",
        signed_week_difference=if(nrow(x) && nrow(y)) y$peak_week-x$peak_week else NA_real_,
        circular_phase_distance=if(periodic) min(abs(y$peak_week-x$peak_week),52-abs(y$peak_week-x$peak_week)) else NA_real_,
        status=if(periodic) "PERIODIC_COMPARISON_CANDIDATE" else "AWAITING_INDEPENDENT_RECURRENCE")
    }
  }
  list(type="phase12_d",status="DIAGNOSTICS_STORED",settings=s,source_hash=source$artifact_hash,
    source_ancestry=ctx$baseline$source_ancestry,baseline_replay=ctx$baseline$replay,
    regime=list(beta="fitted 52-week vector repeated",calendar_start_week=z$calendar_start_week,
      population_denominator=fixed$pop,npieff=z$npieff,initialization="saved fractions multiplied by configured factors",
      maximum_years=z$max_years),plan=grid,checkpoint_diagnostics=cp,annual_summaries=su,
    annual_profiles=pr,initial_condition_comparisons=cmp,
    initial_condition_dependence=p12_bind(dependence),
    transmission_comparisons=p12_bind(transmission),trajectories=trajectories)
}
p12_report <- function(r,out) {
  if(!r$type %in% c("phase12_a","phase12_b","phase12_c","phase12_d")) stop("Saved Phase 12 result required")
  module <- substring(r$type,9L)
  tables <- switch(module,
    a=list(scenario_manifest=r$plan,local_peaks=r$peaks,annual_windows=r$windows,
      annual_outcomes=r$annual,sensitivity_steps=r$sensitivities,sensitivity_comparisons=r$comparisons,
      sensitivity_stability=r$stability),
    b=list(scenario_manifest=r$plan,local_peaks=r$peaks,context_peaks=r$context_peaks,annual_windows=r$windows,
      annual_outcomes=r$annual,chronological_matching=r$matches,competing_heights=r$competing_heights,
      annual_branches=r$branches,full_chronological_matches=r$chronological_matches,
      full_chronological_branches=r$chronological_branches,
      mechanism_intervals=r$terminal,refinement_history=r$refinement_history,
      refinement_decisions=r$refinement_decisions,unresolved_intervals=r$unresolved),
    c=list(scenario_manifest=r$plan,local_peaks=r$peaks,annual_windows=r$windows,
      annual_outcomes=r$annual,neighbor_matching=r$edges,grid_loops=r$loops,
      branch_regions=r$branches,winner_regions=r$winners,branch_interaction=r$interaction,
      global_contrasts=r$global_contrasts,unresolved_cells=r$unresolved),
    d=list(scenario_manifest=r$plan,checkpoint_diagnostics=r$checkpoint_diagnostics,
      annual_summaries=r$annual_summaries,annual_profiles=r$annual_profiles,
      initial_condition_comparisons=r$initial_condition_comparisons,
      initial_condition_dependence=r$initial_condition_dependence,
      transmission_comparisons=r$transmission_comparisons))
  for(name in names(tables)) utils::write.csv(tables[[name]],file.path(out,paste0(name,".csv")),row.names=FALSE)
  grDevices::pdf(file.path(out,paste0("phase12",module,"_diagnostics.pdf")),width=11,height=8)
  tryCatch({
    graphics::plot.new();graphics::title(paste("Phase 12",toupper(module),"saved diagnostic summary"))
    if(module=="a") for(year in r$settings$years) {
      q <- r$sensitivities[r$sensitivities$year==year,,drop=FALSE]
      graphics::plot(q$h,q$sensitivity,col=ifelse(q$axis=="u",2,4),pch=ifelse(q$method=="quadratic",16,1),
        xlab="Log perturbation",ylab="Weeks per unit log change",main=paste(year,"local sensitivity"))
    } else if(module=="b") for(year in r$settings$years) {
      q <- r$terminal[r$terminal$year==year & r$terminal$axis=="v",,drop=FALSE]
      if(nrow(q)) graphics::plot(q$a_left,q$signed_annual_week_difference,col=as.integer(factor(q$classification)),pch=16,
        xlab="Transmission multiplier",ylab="Signed annual-week difference",main=paste(year,"sampled mechanism intervals"))
    } else if(module=="c") for(year in r$settings$years) {
      q <- r$interaction[r$interaction$year==year,,drop=FALSE]
      if(nrow(q)) graphics::plot(q$u,q$v,col=ifelse(q$valid,4,2),pch=16,
        xlab="u",ylab="v",main=paste(year,"computed valid and excluded joint cells"))
    } else {
      q <- r$annual_summaries
      if(nrow(q)) graphics::plot(q$checkpoint,q$peak_week,col=as.integer(factor(q$recurrence_status)),pch=16,
        xlab="Model years",ylab="Peak week",main="Hypothetical long-run checkpoints")
    }
  },finally=grDevices::dev.off())
  writeLines(c(paste0("# Phase 12",toupper(module)," saved diagnostics"),"",
    paste("Source MD5:",r$source_hash),paste("Baseline replay:",r$baseline_replay),
    "Local sensitivity is conditional on the fitted endpoint and numerical peak definition.",
    "Finite timing contrasts are distinct from derivatives.",
    "Numerical recurrence is not proof of an attractor.",
    "Simulated mechanisms do not establish historical causal attribution.",
    "Saved-only reporting reads stored results and performs no simulation or new diagnostics."),
    file.path(out,paste0("phase12",module,"_summary.md")))
  invisible(r)
}
