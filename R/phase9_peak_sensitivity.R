# Numerical peak diagnostics on the accepted Phase 7C endpoint.
p9_defaults <- function() {
  s <- p8_defaults()
  s$time <- list(horizon_rows=NULL,years=2018:2019,include_incomplete=FALSE)
  s$interventions <- list(s_values=NULL,s_multipliers=1,a_values=seq(.94,1.02,by=.0025),
    infected_policy="fixed_baseline_fraction",joint=FALSE,max_scenarios=160L)
  s$annual$enabled <- FALSE
  s$peaks$periodic <- FALSE
  s$peaks$min_relative_prominence <- 0
  s$refinement <- list(enabled=TRUE,levels=4L,max_points=127L)
  s$diagnostics <- list(rapid_week_change=2,rapid_interpolated_change=1,
    flat_curvature_relative=1e-4,interpolation_disagreement=.75,
    match_distance=12,max_log_height_change=log(2),identity_margin=.5,
    interpolation_tolerance=1e-8,stability_tolerance=.5,
    descriptive_eight_week_threshold=8)
  s
}
p9_settings <- function(s,b) {
  expected <- seq(.94,1.02,by=.0025)
  if(!identical(s$time$years,2018:2019) || !is.null(s$time$horizon_rows) ||
     isTRUE(s$time$include_incomplete) || isTRUE(s$annual$enabled) || isTRUE(s$interventions$joint) ||
     !identical(s$interventions$a_values,expected) || !identical(s$interventions$s_multipliers,1) ||
     !is.null(s$interventions$s_values) || isTRUE(s$peaks$periodic))
    stop("Phase 9 requires fixed 2018-2019, a=0.94:0.0025:1.02, saved S0, nonperiodic weekly peaks")
  d <- s$diagnostics
  required <- c("rapid_week_change","rapid_interpolated_change","flat_curvature_relative",
    "interpolation_disagreement","match_distance","max_log_height_change","identity_margin",
    "interpolation_tolerance","stability_tolerance","descriptive_eight_week_threshold")
  if(!identical(sort(names(d)),sort(required))) stop("Phase 9 diagnostic tolerances must be predeclared")
  for(k in required) scalar_number(d[[k]],k,0,strict_lower=TRUE)
  if(d$descriptive_eight_week_threshold!=8) stop("Phase 8B descriptive threshold must remain eight weeks")
  if(s$interventions$max_scenarios<length(expected)) stop("Scenario cap below coarse grid")
  ch <- p8_settings(s,b)
  if(!all(ch$complete)) stop("Phase 9 requires complete 2018 and 2019 annual windows")
  ch
}
p9_interpolate <- function(y,at,d) {
  n <- length(y);out <- list(week=NA_real_,height=NA_real_,status="FAILED")
  if(at<=1L || at>=n) {out$status <- "BOUNDARY";return(list(quadratic=out,cubic=out,curvature=NA_real_))}
  q <- y[(at-1L):(at+1L)];curvature <- q[1]-2*q[2]+q[3]
  quadratic <- out;cubic <- out
  if(is.finite(curvature) && curvature < -d$interpolation_tolerance*max(1,abs(q[2]))) {
    offset <- (q[1]-q[3])/(2*curvature)
    if(is.finite(offset) && abs(offset)<=1+d$interpolation_tolerance) {
      quadratic <- list(week=at+offset,height=q[2]-(q[3]-q[1])^2/(8*curvature),status="OK")
    }
  } else quadratic$status <- "FLAT_OR_NONCONCAVE"
  if(at<=2L || at>=n-1L) cubic$status <- "BOUNDARY_STENCIL" else {
    cubic <- tryCatch({
      f <- stats::splinefun((at-2L):(at+2L),y[(at-2L):(at+2L)],method="monoH.FC")
      opt <- stats::optimize(function(x) -f(x),c(at-1,at+1),tol=1e-8)
      x <- opt$minimum;h <- -opt$objective
      if(!is.finite(h) || abs(x-at)>1+d$interpolation_tolerance || h+ d$interpolation_tolerance*max(1,y[at])<y[at])
        list(week=NA_real_,height=NA_real_,status="INVALID_EXTREMUM")
      else list(week=x,height=h,status="OK")
    },error=function(e) list(week=NA_real_,height=NA_real_,status="SPLINE_FAILURE"))
  }
  list(quadratic=quadratic,cubic=cubic,curvature=curvature)
}
p9_peaks <- function(record,year,d) {
  z <- record$trajectory;at <- which(z$year==year & !z$initialization_only)
  y <- z$incidence[at];annual <- record$annual[[as.character(year)]]
  if(!identical(as.integer(z$week[at]),1:52)) stop("Incomplete Phase 9 annual window")
  local <- annual$local
  if(!nrow(local)) return(data.frame(scenario=character(),year=integer(),a=numeric(),stage=character(),
    level=integer(),peak_id=integer(),weekly_week=numeric(),weekly_height=numeric(),global_week=numeric(),
    global_tie=logical(),retained=logical(),plateau=logical(),boundary=logical(),prominence=numeric(),
    quadratic_week=numeric(),quadratic_height=numeric(),quadratic_status=character(),cubic_week=numeric(),
    cubic_height=numeric(),cubic_status=character(),curvature=numeric(),curvature_relative=numeric(),
    interpolation_disagreement=numeric(),interpolation_uncertain=logical(),seam_case=logical()))
  rows <- lapply(seq_len(nrow(local)),function(i) {
    p <- local[i,];week <- as.integer(p$week);interp <- p9_interpolate(y,week,d)
    if(isTRUE(p$plateau)) {
      interp$quadratic <- list(week=NA_real_,height=NA_real_,status="PLATEAU")
      interp$cubic <- list(week=NA_real_,height=NA_real_,status="PLATEAU")
    }
    data.frame(scenario=record$id,year=year,a=record$a,stage=record$stage,level=record$level,
      peak_id=p$peak_id,weekly_week=p$week,weekly_height=p$height,global_week=annual$global_week,
      global_tie=length(annual$global_ties)>1L,retained=p$retained,plateau=p$plateau,boundary=p$boundary,
      prominence=p$prominence,quadratic_week=interp$quadratic$week,quadratic_height=interp$quadratic$height,
      quadratic_status=interp$quadratic$status,cubic_week=interp$cubic$week,cubic_height=interp$cubic$height,
      cubic_status=interp$cubic$status,curvature=interp$curvature,
      curvature_relative=interp$curvature/max(1,p$height),
      interpolation_disagreement=abs(interp$quadratic$week-interp$cubic$week),
      interpolation_uncertain=interp$quadratic$status!="OK" || interp$cubic$status!="OK" ||
        (is.finite(interp$quadratic$week) && is.finite(interp$cubic$week) &&
           abs(interp$quadratic$week-interp$cubic$week)>d$interpolation_disagreement),
      seam_case=week %in% c(1L,52L))
  })
  do.call(rbind,rows)
}
p9_match <- function(L,R,d) {
  L <- L[L$retained,,drop=FALSE];R <- R[R$retained,,drop=FALSE]
  candidates <- list()
  for(i in seq_len(nrow(L))) for(j in seq_len(nrow(R))) {
    dist <- abs(L$weekly_week[i]-R$weekly_week[j]);ratio <- abs(log(L$weekly_height[i]/R$weekly_height[j]))
    if(dist<=d$match_distance && ratio<=d$max_log_height_change)
      candidates[[length(candidates)+1L]] <- data.frame(i=i,j=j,cost=dist+2*ratio)
  }
  cands <- if(length(candidates)) do.call(rbind,candidates) else data.frame(i=integer(),j=integer(),cost=numeric())
  ambiguous <- rep(FALSE,nrow(cands))
  if(nrow(cands)) for(k in seq_len(nrow(cands))) {
    rival <- cands$cost[cands$i==cands$i[k] | cands$j==cands$j[k]]
    ambiguous[k] <- sum(rival<=cands$cost[k]+d$identity_margin)>1L
  }
  usedL <- usedR <- integer();rows <- list()
  if(nrow(cands)) for(k in order(cands$cost,cands$i,cands$j)) {
    i <- cands$i[k];j <- cands$j[k]
    if(i %in% usedL || j %in% usedR) next
    rows[[length(rows)+1L]] <- data.frame(left_peak=L$peak_id[i],right_peak=R$peak_id[j],
      status=if(ambiguous[k]) "AMBIGUOUS" else "MATCHED",cost=cands$cost[k])
    usedL <- c(usedL,i);usedR <- c(usedR,j)
  }
  for(i in setdiff(seq_len(nrow(L)),usedL)) rows[[length(rows)+1L]] <-
    data.frame(left_peak=L$peak_id[i],right_peak=NA_integer_,status="DEATH_OR_UNMATCHED",cost=NA_real_)
  for(j in setdiff(seq_len(nrow(R)),usedR)) rows[[length(rows)+1L]] <-
    data.frame(left_peak=NA_integer_,right_peak=R$peak_id[j],status="BIRTH_OR_UNMATCHED",cost=NA_real_)
  if(length(rows)) do.call(rbind,rows) else data.frame(left_peak=integer(),right_peak=integer(),status=character(),cost=numeric())
}
p9_edges <- function(records,peaks,d,years=2018:2019) {
  ids <- names(records)[order(vapply(records,function(x) x$a,numeric(1)))];rows <- list();matches <- list();pairs <- list()
  for(year in years) for(k in seq_len(length(ids)-1L)) {
    l <- ids[k];r <- ids[k+1L];L <- peaks[peaks$scenario==l & peaks$year==year,,drop=FALSE]
    R <- peaks[peaks$scenario==r & peaks$year==year,,drop=FALSE]
    m <- p9_match(L,R,d);key <- paste(year,l,r,sep="/")
    if(nrow(m)) {m$year <- year;m$left <- l;m$right <- r;matches[[key]] <- m}
    leftglobal <- records[[l]]$annual[[as.character(year)]]$global_week
    rightglobal <- records[[r]]$annual[[as.character(year)]]$global_week
    good <- m[m$status=="MATCHED",,drop=FALSE]
    crossing <- FALSE
    if(nrow(good)>=2L) for(i in seq_len(nrow(good)-1L)) for(j in (i+1L):nrow(good)) {
      lh <- L$weekly_height[match(good$left_peak[c(i,j)],L$peak_id)]
      rh <- R$weekly_height[match(good$right_peak[c(i,j)],R$peak_id)]
      dl <- diff(lh);dr <- diff(rh);cross <- dl*dr<0
      pairs[[length(pairs)+1L]] <- data.frame(year=year,left=l,right=r,a_left=records[[l]]$a,a_right=records[[r]]$a,
        peak1_left=good$left_peak[i],peak2_left=good$left_peak[j],
        height1_left=lh[1],height2_left=lh[2],height1_right=rh[1],height2_right=rh[2],
        difference_left=dl,difference_right=dr,sign_change=cross,
        confirmed_sign_change=cross && all(m$status=="MATCHED"))
      crossing <- crossing || (cross && all(m$status=="MATCHED"))
    }
    jump <- abs(rightglobal-leftglobal)
    interp <- if(nrow(good)) max(abs(L$quadratic_week[match(good$left_peak,L$peak_id)]-
      R$quadratic_week[match(good$right_peak,R$peak_id)]),na.rm=TRUE) else NA_real_
    if(!is.finite(interp)) interp <- NA_real_
    uncertainty <- any(m$status!="MATCHED") || any(L$interpolation_uncertain) || any(R$interpolation_uncertain) ||
      any(L$global_tie) || any(R$global_tie)
    lwin <- L$peak_id[L$weekly_week==leftglobal & L$retained]
    rwin <- R$peak_id[R$weekly_week==rightglobal & R$retained]
    winner_exchange <- FALSE
    if(length(lwin)==1L && length(rwin)==1L && !uncertainty) {
      winner_exchange <- !any(good$left_peak==lwin & good$right_peak==rwin) &&
        any(good$left_peak==lwin) && any(good$right_peak==rwin)
    }
    confirmed_switch <- crossing && winner_exchange
    rapid <- jump>=d$rapid_week_change || (!is.na(interp) && interp>=d$rapid_interpolated_change)
    flat <- any(abs(c(L$curvature_relative,R$curvature_relative))<=d$flat_curvature_relative,na.rm=TRUE)
    trigger <- rapid || flat || crossing || uncertainty || nrow(L)!=nrow(R)
    rows[[length(rows)+1L]] <- data.frame(year=year,left=l,right=r,a_left=records[[l]]$a,a_right=records[[r]]$a,
      width=records[[r]]$a-records[[l]]$a,weekly_global_displacement=jump,
      eight_week_descriptive=jump>=d$descriptive_eight_week_threshold,
      matched_interpolated_displacement=interp,matched=sum(m$status=="MATCHED"),
      ambiguous=sum(m$status=="AMBIGUOUS"),births_deaths=sum(m$status %in% c("BIRTH_OR_UNMATCHED","DEATH_OR_UNMATCHED")),
      peak_count_change=nrow(L)!=nrow(R),flat=flat,confirmed_height_crossing=crossing,
      confirmed_weekly_switch=confirmed_switch,
      uncertainty=uncertainty,refine=trigger,
      classification=if(uncertainty) "UNRESOLVED" else if(confirmed_switch) "CONFIRMED_WEEKLY_SWITCH" else if(rapid) "PERSISTENT_RAPID_CANDIDATE" else "NO_RAPID_CANDIDATE")
  }
  list(edges=do.call(rbind,rows),matches=if(length(matches)) do.call(rbind,matches) else data.frame(),
    pairs=if(length(pairs)) do.call(rbind,pairs) else data.frame())
}
p9_refine <- function(records,baseline,check,profile,s) {
  d <- s$diagnostics;history <- list();unresolved <- list();added <- 0L
  for(level in 0:s$refinement$levels) {
    peaks <- do.call(rbind,lapply(records,function(r) do.call(rbind,lapply(check$years,function(y) p9_peaks(r,y,d)))))
    edges <- p9_edges(records,peaks,d,check$years);e <- edges$edges;e$level <- level
    history[[length(history)+1L]] <- e
    if(level==s$refinement$levels || !isTRUE(s$refinement$enabled)) break
    candidates <- e[e$refine,,drop=FALSE]
    keys <- paste(candidates$left,candidates$right)
    candidates <- candidates[!duplicated(keys),,drop=FALSE]
    new <- 0L
    for(i in seq_len(nrow(candidates))) {
      row <- candidates[i,];mid <- (row$a_left+row$a_right)/2
      if(mid %in% vapply(records,function(r) r$a,numeric(1))) next
      if(length(records)>=s$interventions$max_scenarios || added>=s$refinement$max_points) {
        unresolved[[length(unresolved)+1L]] <- data.frame(year=row$year,a_left=row$a_left,a_right=row$a_right,level=level,reason="SCENARIO_CAP")
        next
      }
      id <- sprintf("scenario-%04d",length(records)+1L)
      r <- p8_evaluate(baseline,check,profile,baseline$parameters$S0,mid,id,"a","refined",level+1L,paste(row$left,row$right,sep="/"))
      if(r$status!="VALID") stop("Phase 9 scenario failed: ",r$reason)
      records[[id]] <- r;added <- added+1L;new <- new+1L
    }
    if(new==0L) break
  }
  peaks <- do.call(rbind,lapply(records,function(r) do.call(rbind,lapply(check$years,function(y) p9_peaks(r,y,d)))))
  terminal <- p9_edges(records,peaks,d,check$years)
  if(any(terminal$edges$refine)) {
    t <- terminal$edges[terminal$edges$refine,,drop=FALSE]
    for(i in seq_len(nrow(t))) unresolved[[length(unresolved)+1L]] <- data.frame(year=t$year[i],a_left=t$a_left[i],a_right=t$a_right[i],
      level=max(vapply(history,function(x) max(x$level),numeric(1))),reason="REFINEMENT_LIMIT_OR_CAP")
  }
  list(records=records,peaks=peaks,terminal=terminal,history=do.call(rbind,history),
    unresolved=if(length(unresolved)) unique(do.call(rbind,unresolved)) else data.frame())
}
p9_branches <- function(records,peaks,matches,d) {
  ids <- names(records)[order(vapply(records,function(x) x$a,numeric(1)))];out <- list()
  for(year in 2018:2019) {
    prior <- character();nextid <- 0L
    for(k in seq_along(ids)) {
      p <- peaks[peaks$scenario==ids[k] & peaks$year==year & peaks$retained,,drop=FALSE]
      current <- setNames(rep(NA_character_,nrow(p)),as.character(p$peak_id))
      if(k>1L && nrow(matches)) {
        m <- matches[matches$year==year & matches$left==ids[k-1L] & matches$right==ids[k] & matches$status=="MATCHED",,drop=FALSE]
        for(j in seq_len(nrow(m))) current[as.character(m$right_peak[j])] <- prior[as.character(m$left_peak[j])]
      }
      for(j in which(is.na(current))) {nextid <- nextid+1L;current[j] <- paste0(year,"/branch-",nextid)}
      if(nrow(p)) {p$branch <- unname(current[as.character(p$peak_id)]);out[[length(out)+1L]] <- p}
      prior <- current
    }
  }
  b <- if(length(out)) do.call(rbind,out) else data.frame()
  if(!nrow(b)) return(b)
  b$dt_da <- NA_real_;b$difference_method <- "NA_UNMATCHED";b$grid_stable <- NA
  rank <- setNames(seq_along(ids),ids)
  for(id in unique(b$branch)) {
    ix <- which(b$branch==id);ix <- ix[order(b$a[ix])]
    if(length(ix)<2L) next
    for(k in seq_along(ix)) {
      j <- ix[k]
      left_ok <- k>1L && rank[b$scenario[ix[k-1L]]]==rank[b$scenario[j]]-1L
      right_ok <- k<length(ix) && rank[b$scenario[ix[k+1L]]]==rank[b$scenario[j]]+1L
      if(!left_ok && !right_ok) next
      neighbors <- if(left_ok && right_ok) c(ix[k-1L],ix[k+1L]) else if(right_ok) c(j,ix[k+1L]) else c(ix[k-1L],j)
      b$dt_da[j] <- diff(b$quadratic_week[neighbors])/diff(b$a[neighbors])
      b$difference_method[j] <- if(left_ok && right_ok) "CENTRAL_MATCHED" else "ONE_SIDED_MATCHED"
      b$grid_stable[j] <- if(left_ok && right_ok) {
        predicted <- b$quadratic_week[ix[k-1L]] +
          (b$a[j]-b$a[ix[k-1L]])/(b$a[ix[k+1L]]-b$a[ix[k-1L]]) *
          (b$quadratic_week[ix[k+1L]]-b$quadratic_week[ix[k-1L]])
        abs(b$quadratic_week[j]-predicted)<=d$stability_tolerance
      } else NA
    }
  }
  b
}
p9_run <- function(source,data,s=p9_defaults(),root=".") {
  baseline <- p8_source(source,data,s,root);check <- p9_settings(s,baseline)
  profile <- p8_profile(baseline,s);profile$kind <- "fitted"
  records <- list()
  for(a in check$a_values) {
    id <- sprintf("scenario-%04d",length(records)+1L)
    r <- p8_evaluate(baseline,check,profile,baseline$parameters$S0,a,id,"a","coarse")
    if(r$status!="VALID") stop("Phase 9 scenario failed: ",r$reason)
    records[[id]] <- r
  }
  base <- records[[which(check$a_values==1)]]$trajectory
  ref <- baseline$evaluation$states[seq_len(check$horizon),,drop=FALSE]
  if(!p6r_close(as.numeric(as.matrix(base[,names(ref)])),as.numeric(as.matrix(ref)),s$replay_tolerance)) stop("Exact baseline replay failed")
  x <- p9_refine(records,baseline,check,profile,s)
  x$branches <- p9_branches(x$records,x$peaks,x$terminal$matches,s$diagnostics)
  x$type <- "phase9_peak_sensitivity";x$status <- "FINITE_WEEKLY_DIAGNOSTICS_STORED"
  x$diagnostic_settings <- s;x$source_hash <- source$artifact_hash;x$source_identity <- source$identity
  x$source_ancestry <- baseline$source_ancestry;x$baseline_replay <- baseline$replay
  x$time_contract <- check
  x
}
p9_report <- function(result,out) {
  if(!identical(result$type,"phase9_peak_sensitivity")) stop("Saved Phase 9 result required")
  write <- function(x,name) utils::write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE)
  scenarios <- do.call(rbind,lapply(result$records,function(r) data.frame(scenario=r$id,a=r$a,s=r$s,stage=r$stage,level=r$level,parent=r$parent,status=r$status)))
  write(scenarios,"scenarios");write(result$peaks,"local_peaks_interpolation")
  write(result$branches,"branch_sensitivity");write(result$history,"refinement_history")
  write(result$terminal$edges,"terminal_intervals");write(result$terminal$matches,"peak_matching")
  write(result$terminal$pairs,"competing_heights");write(result$unresolved,"unresolved_intervals")
  grDevices::pdf(file.path(out,"phase9_diagnostics.pdf"),width=11,height=8)
  tryCatch({for(year in 2018:2019) {
    b <- if(nrow(result$branches)) result$branches[result$branches$year==year,,drop=FALSE] else data.frame()
    p <- result$peaks[result$peaks$year==year,,drop=FALSE]
    e <- result$terminal$edges[result$terminal$edges$year==year,,drop=FALSE]
    if(nrow(p)) {
      graphics::plot(p$a,p$weekly_week,pch=16,cex=.5,xlab="Transmission multiplier a",ylab="Weekly local peak week",main=paste(year,"weekly peaks and matched branches"))
      for(id in unique(b$branch)) {v <- b[b$branch==id,];v <- v[order(v$a),];graphics::lines(v$a,v$weekly_week,col=2)}
      if(any(is.finite(p$quadratic_week))) {
        graphics::plot(p$a,p$quadratic_week,pch=16,cex=.5,xlab="a",ylab="Estimated peak week",main=paste(year,"quadratic (black), cubic (blue), weekly (gray)"))
        graphics::points(p$a,p$cubic_week,col=4,pch=1);graphics::points(p$a,p$weekly_week,col="gray",pch=3)
      } else {graphics::plot.new();graphics::title(main=paste(year,"no valid interpolated peaks"))}
    } else {graphics::plot.new();graphics::title(main=paste(year,"no weekly local peaks"));graphics::plot.new()}
    finite <- is.finite(b$curvature_relative) & is.finite(b$dt_da)
    if(any(finite)) graphics::plot(b$curvature_relative[finite],b$dt_da[finite],pch=16,xlab="Relative local curvature",ylab="Estimated dt*/da",main=paste(year,"sensitivity vs curvature"))
    else {graphics::plot.new();graphics::title(main=paste(year,"no finite matched sensitivity estimates"))}
    pair <- result$terminal$pairs
    if(nrow(pair)) pair <- pair[pair$year==year,,drop=FALSE]
    if(nrow(pair)) graphics::plot(pair$a_left,pair$difference_left,pch=16,xlab="a",ylab="Competing height difference",main=paste(year,"matched competing peak heights")) else graphics::plot.new()
    graphics::plot(e$a_left,seq_len(nrow(e)),pch=ifelse(e$uncertainty,4,16),xlab="a",ylab="Interval",main=paste(year,"flagged intervals; x means uncertainty"))
    graphics::segments(e$a_left,seq_len(nrow(e)),e$a_right,seq_len(nrow(e)),col=ifelse(e$refine,2,1))
  }},finally=grDevices::dev.off())
  writeLines(c("# Phase 9 peak displacement and sensitivity diagnostics","",
    paste("Accepted Phase 7C source MD5:",result$source_hash),paste("Replay:",result$baseline_replay),
    paste("Scenarios:",length(result$records)),paste("Unresolved intervals:",nrow(result$unresolved)),"",
    "Observed weekly peaks are distinct from quadratic and shape-preserving cubic numerical estimates. Neither interpolation estimates an actual continuous-time SIRS peak.",
    "Finite differences use confidently matched branches only. One-sided estimates and ambiguous identities are labeled. Calendar years are analyzed without periodic wrapping; SIRS state continues across years.",
    "The eight-week Phase 8B threshold is descriptive and separate from refinement triggers. No causal, bifurcation, mathematical-discontinuity, or certified-convergence claim is made.",
    "Saved-only reporting uses the stored numerical diagnostics and performs no simulation.","",
    "## Predeclared settings","```",capture.output(dput(result$diagnostic_settings)),"```"),file.path(out,"phase9_summary.md"))
  invisible(result)
}
