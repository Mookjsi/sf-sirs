# Full-trajectory weekly peak identity diagnostics. No model equations are changed.
p10_defaults <- function() {
  s <- p9_defaults()
  s$refinement <- list(enabled=TRUE,levels=4L,max_points=80L,per_year=40L)
  s$interventions$max_scenarios <- 113L
  s$identity <- list(match_distance=12,max_log_height_change=log(2),height_weight=2,
    context_weight=1,margin=.5,near_height_fraction=.02,rapid_week_change=2,
    curvature_tolerance=1e-8)
  s
}
p10_settings <- function(s,b) {
  old <- s;s$refinement <- p9_defaults()$refinement
  check <- p9_settings(s,b);s <- old
  d <- s$identity
  if(!identical(sort(names(d)),sort(c("match_distance","max_log_height_change","height_weight",
    "context_weight","margin","near_height_fraction","rapid_week_change","curvature_tolerance")))) stop("Phase 10 identity settings must be predeclared")
  for(k in names(d)) scalar_number(d[[k]],k,0)
  integer_scalar(s$refinement$per_year,"per-year refinement cap",0L)
  if(s$refinement$max_points>s$interventions$max_scenarios-length(check$a_values)) stop("Refinement point cap exceeds scenario cap")
  check$settings <- old;check
}
p10_detect <- function(record,year,policy,d) {
  z <- record$trajectory;y <- z$incidence;valid <- which(!z$initialization_only)
  if(!length(valid) || any(!is.finite(y[valid]))) stop("Finite supported weekly incidence required")
  time <- seq_along(y);tol <- p4_tolerance(y[valid],policy)
  raw <- p4_weekly_peaks(y[valid],time[valid],policy)$local
  rows <- list()
  for(i in seq_len(nrow(raw))) {
    p <- raw[i,];at <- as.integer(p$first);last <- as.integer(p$last)
    if(z$year[at]!=year) next
    support <- at>min(valid) && last<max(valid)
    genuine <- support && !p$plateau
    left_floor <- if(at>min(valid)) min(y[min(valid):at]) else NA_real_
    right_floor <- if(last<max(valid)) min(y[last:max(valid)]) else NA_real_
    prominence <- if(support) max(0,p$height-max(left_floor,right_floor)) else NA_real_
    relative <- prominence/max(p$height,.Machine$double.eps)
    retained <- genuine && is.finite(prominence) && prominence>=policy$min_prominence && relative>=policy$min_relative_prominence
    q <- p10_interpolate(y,at,last,d)
    rows[[length(rows)+1L]] <- data.frame(scenario=record$id,year=year,a=record$a,peak_id=at,
      position=at,week=z$week[at],height=p$height,first=at,last=last,
      left_height=if(at>min(valid)) y[at-1L] else NA_real_,
      right_height=if(last<max(valid)) y[last+1L] else NA_real_,
      support=support,plateau=p$plateau,genuine=genuine,prominence=prominence,
      relative_prominence=relative,retained=retained,
      quadratic_position=q$position,quadratic_height=q$height,quadratic_status=q$status,
      curvature=q$curvature,cubic_position=q$cubic_position,cubic_height=q$cubic_height,
      cubic_status=q$cubic_status,stringsAsFactors=FALSE)
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p10_interpolate <- function(y,at,last,d) {
  out <- list(position=NA_real_,height=NA_real_,status="UNSUPPORTED",curvature=NA_real_,
    cubic_position=NA_real_,cubic_height=NA_real_,cubic_status="UNSUPPORTED")
  if(at!=last) {out$status <- out$cubic_status <- "PLATEAU";return(out)}
  if(at<=1L || at>=length(y)) return(out)
  q <- y[(at-1L):(at+1L)];curv <- q[1]-2*q[2]+q[3];out$curvature <- curv
  if(is.finite(curv) && curv < -d$curvature_tolerance*max(1,abs(q[2]))) {
    offset <- (q[1]-q[3])/(2*curv)
    if(is.finite(offset) && abs(offset)<=1+1e-8) {
      out$position <- at+offset;out$height <- q[2]-(q[3]-q[1])^2/(8*curv);out$status <- "OK"
    } else out$status <- "EXTREMUM_OUTSIDE_STENCIL"
  } else out$status <- "FLAT_OR_NONCONCAVE"
  if(at<=2L || at>=length(y)-1L) {out$cubic_status <- "BOUNDARY_STENCIL";return(out)}
  cubic <- tryCatch(stats::splinefun((at-2L):(at+2L),y[(at-2L):(at+2L)],method="monoH.FC"),error=function(e) NULL)
  if(is.null(cubic)) {out$cubic_status <- "SPLINE_FAILURE";return(out)}
  probe <- cubic(c(at-.5,at,at+.5))
  if(all(is.finite(probe)) && max(probe)-y[at]<=1e-7*max(1,y[at])) {
    out$cubic_position <- at;out$cubic_height <- y[at];out$cubic_status <- "KNOT_MAXIMUM_AUDIT"
  } else out$cubic_status <- "SHAPE_AUDIT_FAILURE"
  out
}
p10_full_count <- function(record,policy) {
  z <- record$trajectory;v <- which(!z$initialization_only)
  p <- p4_weekly_peaks(z$incidence[v],v,policy)$local
  sum(!p$plateau & !p$boundary)
}
p10_window <- function(record,year,peaks,policy) {
  z <- record$trajectory;ix <- which(z$year==year & !z$initialization_only)
  yy <- z$incidence[ix];tol <- p4_tolerance(yy,policy);top <- which(max(yy)-yy<=tol)
  at <- ix[top[1]];p <- peaks[peaks$position==at,,drop=FALSE]
  status <- if(length(top)>1L) "TIE" else if(!nrow(p)) {
    if(at==min(ix) || at==max(ix)) "WINDOW_SELECTION" else "NOT_LOCAL_MAXIMUM"
  } else if(p$plateau) "PLATEAU" else if(!p$support) "SUPPORT_BOUNDARY" else if(p$genuine) "GENUINE_LOCAL_PEAK" else "UNRESOLVED"
  data.frame(scenario=record$id,year=year,a=record$a,position=at,week=z$week[at],
    height=z$incidence[at],ties=paste(z$week[ix[top]],collapse=";"),status=status,
    peak_id=if(nrow(p)) p$peak_id[1] else NA_integer_,stringsAsFactors=FALSE)
}
p10_assignment <- function(L,R,d,reverse=FALSE) {
  if(reverse) {tmp <- L;L <- R;R <- tmp}
  L <- L[L$genuine,,drop=FALSE];R <- R[R$genuine,,drop=FALSE]
  if(!nrow(L) || !nrow(R)) return(data.frame(left_peak=integer(),right_peak=integer(),cost=numeric(),margin=numeric(),status=character(),direction=character()))
  distance <- abs(outer(L$position,R$position,`-`))
  ratio <- abs(outer(log(pmax(L$height,1e-300)),log(pmax(R$height,1e-300)),`-`))
  context <- abs(outer((L$left_height+L$right_height)/pmax(L$height,1e-300),
    (R$left_height+R$right_height)/pmax(R$height,1e-300),`-`))
  cost <- distance+d$height_weight*ratio+d$context_weight*context
  allowed <- distance<=d$match_distance & ratio<=d$max_log_height_change & is.finite(cost)
  safe <- cost;safe[!is.finite(safe)] <- 1e12
  solve <- function(mask) p5_circular_peak_match(L$position,R$position,d$match_distance,period=1e6,
    edge_allowed=mask,edge_cost=safe)
  fit <- solve(allowed);matched <- fit[fit$status=="MATCHED",,drop=FALSE]
  total <- function(tab) {v <- tab[tab$status=="MATCHED",,drop=FALSE];if(!nrow(v)) return(0)
    sum(cost[cbind(match(v$reference_x,L$position),match(v$target_x,R$position))])}
  base <- total(fit);rows <- list()
  for(k in seq_len(nrow(matched))) {
    i <- match(matched$reference_x[k],L$position);j <- match(matched$target_x[k],R$position)
    mask <- allowed;mask[i,j] <- FALSE;alt <- solve(mask)
    margin <- if(sum(alt$status=="MATCHED")==nrow(matched)) total(alt)-base else Inf
    rows[[length(rows)+1L]] <- data.frame(left_peak=L$peak_id[i],right_peak=R$peak_id[j],
      cost=cost[i,j],margin=margin,status=if(margin<=d$margin+1e-10) "AMBIGUOUS" else "MATCHED",
      direction=if(reverse) "reverse" else "forward")
  }
  for(i in setdiff(seq_len(nrow(L)),match(matched$reference_x,L$position))) rows[[length(rows)+1L]] <-
    data.frame(left_peak=L$peak_id[i],right_peak=NA_integer_,cost=NA_real_,margin=NA_real_,status="UNMATCHED_LEFT",direction=if(reverse) "reverse" else "forward")
  for(j in setdiff(seq_len(nrow(R)),match(matched$target_x,R$position))) rows[[length(rows)+1L]] <-
    data.frame(left_peak=NA_integer_,right_peak=R$peak_id[j],cost=NA_real_,margin=NA_real_,status="UNMATCHED_RIGHT",direction=if(reverse) "reverse" else "forward")
  out <- if(length(rows)) do.call(rbind,rows) else data.frame()
  if(reverse && nrow(out)) {tmp <- out$left_peak;out$left_peak <- out$right_peak;out$right_peak <- tmp}
  out
}
p10_candidates <- function(L,R,d) {
  L <- L[L$genuine,,drop=FALSE];R <- R[R$genuine,,drop=FALSE]
  rows <- list()
  for(i in seq_len(nrow(L))) for(j in seq_len(nrow(R))) {
    distance <- abs(L$position[i]-R$position[j])
    ratio <- abs(log(pmax(L$height[i],1e-300)/pmax(R$height[j],1e-300)))
    context <- abs((L$left_height[i]+L$right_height[i])/pmax(L$height[i],1e-300)-
      (R$left_height[j]+R$right_height[j])/pmax(R$height[j],1e-300))
    rows[[length(rows)+1L]] <- data.frame(left_peak=L$peak_id[i],right_peak=R$peak_id[j],
      distance=distance,log_height_change=ratio,context_change=context,
      cost=distance+d$height_weight*ratio+d$context_weight*context,
      admissible=distance<=d$match_distance && ratio<=d$max_log_height_change && is.finite(context))
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p10_identity_consistent <- function(forward,reverse) {
  key <- function(x) paste(x$left_peak[x$status=="MATCHED"],x$right_peak[x$status=="MATCHED"],sep="/")
  setequal(key(forward),key(reverse)) &&
    sum(forward$status=="AMBIGUOUS")==sum(reverse$status=="AMBIGUOUS")
}
p10_edges <- function(records,peaks,windows,s) {
  ids <- names(records)[order(vapply(records,function(r) r$a,numeric(1)))];rows <- list();matches <- list();pairs <- list();candidates <- list()
  for(year in sort(unique(windows$year))) for(k in seq_len(length(ids)-1L)) {
    l <- ids[k];r <- ids[k+1L];L <- peaks[peaks$scenario==l & peaks$year==year,,drop=FALSE]
    R <- peaks[peaks$scenario==r & peaks$year==year,,drop=FALSE]
    c <- p10_candidates(L,R,s$identity)
    if(nrow(c)) {c$year <- year;c$left <- l;c$right <- r;c$a_left <- records[[l]]$a;c$a_right <- records[[r]]$a;candidates[[length(candidates)+1L]] <- c}
    f <- p10_assignment(L,R,s$identity);v <- p10_assignment(L,R,s$identity,TRUE)
    consistent <- p10_identity_consistent(f,v)
    for(direction in list(f,v)) if(nrow(direction)) {direction$year <- year;direction$left <- l;direction$right <- r;direction$a_left <- records[[l]]$a;direction$a_right <- records[[r]]$a;direction$direction_consistent <- consistent;matches[[length(matches)+1L]] <- direction}
    wl <- windows[windows$scenario==l & windows$year==year,,drop=FALSE];wr <- windows[windows$scenario==r & windows$year==year,,drop=FALSE]
    good <- f[f$status=="MATCHED",,drop=FALSE]
    if(!consistent) good <- good[FALSE,,drop=FALSE]
    winner_l <- if(wl$status=="GENUINE_LOCAL_PEAK") wl$peak_id else NA_integer_
    winner_r <- if(wr$status=="GENUINE_LOCAL_PEAK") wr$peak_id else NA_integer_
    same <- any(good$left_peak==winner_l & good$right_peak==winner_r)
    pair_switch <- FALSE;near <- FALSE
    if(nrow(good)>=2L) for(i in seq_len(nrow(good)-1L)) for(j in (i+1L):nrow(good)) {
      x <- L[match(good$left_peak[c(i,j)],L$peak_id),,drop=FALSE]
      y <- R[match(good$right_peak[c(i,j)],R$peak_id),,drop=FALSE]
      dl <- x$height[1]-x$height[2];dr <- y$height[1]-y$height[2]
      tol <- s$peaks$tie_abs+s$peaks$tie_rel*max(c(x$height,y$height))
      crossing <- dl*dr<0 && abs(dl)>tol && abs(dr)>tol
      exchange <- crossing && !same && winner_l %in% x$peak_id && winner_r %in% y$peak_id &&
        x$peak_id[which.max(x$height)]==winner_l && y$peak_id[which.max(y$height)]==winner_r
      pair_switch <- pair_switch || exchange
      near <- near || min(abs(dl),abs(dr))/max(c(x$height,y$height))<=s$identity$near_height_fraction
      pairs[[length(pairs)+1L]] <- data.frame(year=year,left=l,right=r,a_left=records[[l]]$a,a_right=records[[r]]$a,
        branch1_left_peak=x$peak_id[1],branch2_left_peak=x$peak_id[2],
        branch1_right_peak=y$peak_id[1],branch2_right_peak=y$peak_id[2],
        branch1_left_week=x$week[1],branch2_left_week=x$week[2],
        branch1_right_week=y$week[1],branch2_right_week=y$week[2],
        branch1_left_height=x$height[1],branch2_left_height=x$height[2],
        branch1_right_height=y$height[1],branch2_right_height=y$height[2],
        difference_left=dl,difference_right=dr,weekly_order_reversal=crossing,
        global_exchange=exchange,identity_consistent=consistent,
        winner_left=winner_l,winner_right=winner_r,window_left=wl$status,window_right=wr$status,
        tie_left=wl$status=="TIE",tie_right=wr$status=="TIE",
        matching_uncertain=!consistent || any(f$status=="AMBIGUOUS") || any(v$status=="AMBIGUOUS"))
    }
    raw_change <- p10_full_count(records[[l]],s$peaks)!=p10_full_count(records[[r]],s$peaks)
    window_crossing <- sum(L$genuine)!=sum(R$genuine) && !raw_change
    filtered_change <- sum(L$retained)!=sum(R$retained)
    uncertain <- !consistent || any(f$status=="AMBIGUOUS") || any(v$status=="AMBIGUOUS") ||
      wl$status!="GENUINE_LOCAL_PEAK" || wr$status!="GENUINE_LOCAL_PEAK"
    movement <- abs(wr$position-wl$position)
    rapid <- same && movement>=s$identity$rapid_week_change
    flags <- c(winner_exchange=pair_switch,near_height=near,rapid_displacement=rapid,
      winner_identity_ambiguity=uncertain,raw_topology_change=raw_change,
      threshold_exclusion=filtered_change && !raw_change,window_crossing=window_crossing,window_selection=wl$status=="WINDOW_SELECTION" || wr$status=="WINDOW_SELECTION")
    priority <- if(pair_switch) 1L else if(near) 2L else if(rapid) 3L else if(uncertain && movement>=s$identity$rapid_week_change) 4L else if(raw_change) 5L else 99L
    class <- if(pair_switch && !uncertain) "CONFIRMED_WEEKLY_SWITCH" else if(same && !uncertain && movement>0) "PERSISTENT_BRANCH_DISPLACEMENT" else if(raw_change && !uncertain) "SAMPLED_TOPOLOGY_CANDIDATE" else if(any(flags) || movement>0) "UNRESOLVED" else "NO_DETECTED_MOVEMENT"
    rows[[length(rows)+1L]] <- data.frame(year=year,left=l,right=r,a_left=records[[l]]$a,a_right=records[[r]]$a,
      weekly_argmax_movement=movement,accumulated_movement=NA_real_,same_winner_branch=same,
      raw_left=sum(L$genuine),raw_right=sum(R$genuine),retained_left=sum(L$retained),retained_right=sum(R$retained),
      identity_consistent=consistent,matching_uncertain=uncertain,window_left=wl$status,window_right=wr$status,
      flag_winner_exchange=flags[1],flag_near_height=flags[2],flag_rapid_displacement=flags[3],
      flag_winner_identity_ambiguity=flags[4],flag_raw_topology_change=flags[5],
      flag_threshold_exclusion=flags[6],flag_window_crossing=flags[7],flag_window_selection=flags[8],priority=priority,
      refine=priority<99L,classification=class,stringsAsFactors=FALSE)
  }
  list(edges=do.call(rbind,rows),matches=if(length(matches)) do.call(rbind,matches) else data.frame(),
    pairs=if(length(pairs)) do.call(rbind,pairs) else data.frame(),
    candidates=if(length(candidates)) do.call(rbind,candidates) else data.frame())
}
p10_collect <- function(records,s) {
  peaks <- list();windows <- list()
  for(r in records) for(year in 2018:2019) {
    if(!any(r$trajectory$year==year & !r$trajectory$initialization_only)) next
    p <- p10_detect(r,year,s$peaks,s$identity)
    if(nrow(p)) peaks[[length(peaks)+1L]] <- p
    windows[[length(windows)+1L]] <- p10_window(r,year,p,s$peaks)
  }
  list(peaks=if(length(peaks)) do.call(rbind,peaks) else data.frame(),windows=do.call(rbind,windows))
}
p10_refine <- function(records,baseline,check,profile,s) {
  history <- list();decisions <- list();counts <- setNames(c(0L,0L),c("2018","2019"));added <- 0L
  for(level in 0:s$refinement$levels) {
    collected <- p10_collect(records,s);edge <- p10_edges(records,collected$peaks,collected$windows,s)
    e <- edge$edges;e$level <- level;history[[length(history)+1L]] <- e
    if(level==s$refinement$levels || !s$refinement$enabled) break
    candidates <- e[e$refine,,drop=FALSE]
    candidates <- candidates[order(candidates$priority,candidates$year,candidates$a_left),,drop=FALSE]
    new <- 0L
    for(i in seq_len(nrow(candidates))) {
      x <- candidates[i,];mid <- (x$a_left+x$a_right)/2
      existing <- any(abs(vapply(records,function(r) r$a,numeric(1))-mid)<1e-14)
      reason <- if(existing) "SHARED_MIDPOINT" else if(counts[as.character(x$year)]>=s$refinement$per_year) "YEAR_CAP" else
        if(added>=s$refinement$max_points || length(records)>=s$interventions$max_scenarios) "SCENARIO_CAP" else "SELECTED"
      decisions[[length(decisions)+1L]] <- data.frame(year=x$year,level=level,a_left=x$a_left,a_right=x$a_right,
        priority=x$priority,decision=reason,flags=paste(names(x)[vapply(x,function(v) is.logical(v) && length(v)==1L && isTRUE(v),logical(1)) & grepl("^flag_",names(x))],collapse=";"))
      if(reason!="SELECTED") next
      id <- sprintf("scenario-%04d",length(records)+1L)
      r <- p8_evaluate(baseline,check,profile,baseline$parameters$S0,mid,id,"a","refined",level+1L,paste(x$left,x$right,sep="/"))
      if(r$status!="VALID") stop("Phase 10 scenario failed: ",r$reason)
      records[[id]] <- r;counts[as.character(x$year)] <- counts[as.character(x$year)]+1L;added <- added+1L;new <- new+1L
    }
    if(!new) break
  }
  collected <- p10_collect(records,s);terminal <- p10_edges(records,collected$peaks,collected$windows,s)
  unresolved <- terminal$edges[terminal$edges$refine | terminal$edges$classification=="UNRESOLVED",,drop=FALSE]
  if(nrow(unresolved)) unresolved$reason <- ifelse(unresolved$refine,
    ifelse(counts[as.character(unresolved$year)]>=s$refinement$per_year,"YEAR_CAP",
      ifelse(added>=s$refinement$max_points || length(records)>=s$interventions$max_scenarios,"SCENARIO_CAP","REFINEMENT_LIMIT")),
    "IDENTITY_OR_WINDOW_UNRESOLVED")
  list(records=records,peaks=collected$peaks,windows=collected$windows,terminal=terminal,
    history=do.call(rbind,history),decisions=if(length(decisions)) do.call(rbind,decisions) else data.frame(),
    unresolved=unresolved,year_additions=counts)
}
p10_three_point <- function(x,y,at) {
  if(length(x)!=3L || any(!is.finite(y)) || any(diff(x)<=0)) return(NA_real_)
  (y[1]*(2*x[at]-x[2]-x[3])/((x[1]-x[2])*(x[1]-x[3]))+
   y[2]*(2*x[at]-x[1]-x[3])/((x[2]-x[1])*(x[2]-x[3]))+
   y[3]*(2*x[at]-x[1]-x[2])/((x[3]-x[1])*(x[3]-x[2])))
}
p10_branches <- function(records,peaks,matches,windows,s) {
  ids <- names(records)[order(vapply(records,function(r) r$a,numeric(1)))];rows <- list()
  for(year in 2018:2019) {
    prior <- character();next_id <- 0L
    for(k in seq_along(ids)) {
      p <- peaks[peaks$scenario==ids[k] & peaks$year==year & peaks$genuine,,drop=FALSE]
      current <- setNames(rep(NA_character_,nrow(p)),as.character(p$peak_id))
      if(k>1L && nrow(matches)) {
        m <- matches[matches$year==year & matches$left==ids[k-1L] & matches$right==ids[k] &
          matches$direction=="forward" & matches$status=="MATCHED" & matches$direction_consistent,,drop=FALSE]
        for(j in seq_len(nrow(m))) current[as.character(m$right_peak[j])] <- prior[as.character(m$left_peak[j])]
      }
      for(j in which(is.na(current))) {next_id <- next_id+1L;current[j] <- paste0(year,"/branch-",next_id)}
      if(nrow(p)) {p$branch <- unname(current[as.character(p$peak_id)]);rows[[length(rows)+1L]] <- p}
      prior <- current
    }
  }
  b <- if(length(rows)) do.call(rbind,rows) else data.frame()
  if(!nrow(b)) return(b)
  b$winner <- mapply(function(sc,yr,id) any(windows$scenario==sc & windows$year==yr & windows$peak_id==id & windows$status=="GENUINE_LOCAL_PEAK"),b$scenario,b$year,b$peak_id)
  b$accumulated_weekly_movement <- NA_real_
  for(branch in unique(b$branch)) {at <- which(b$branch==branch);at <- at[order(b$a[at])];b$accumulated_weekly_movement[at] <- b$position[at]-b$position[at[1]]}
  b$dt_da <- b$dt_dloga <- NA_real_;b$sensitivity_method <- "NA_INADEQUATE_IDENTITY_OR_INTERPOLATION"
  b$sensitivity_change_under_refinement <- NA_real_
  b$nested_sensitivity_change <- NA_real_
  b$response_numerator <- b$response_curvature_denominator <- b$surrogate_dt_da <- NA_real_
  rank <- setNames(seq_along(ids),ids)
  for(branch in unique(b$branch)) {
    ix <- which(b$branch==branch);ix <- ix[order(b$a[ix])]
    for(k in seq_along(ix)) {
      j <- ix[k];if(b$quadratic_status[j]!="OK") next
      left <- k>1L && rank[b$scenario[ix[k-1L]]]==rank[b$scenario[j]]-1L
      right <- k<length(ix) && rank[b$scenario[ix[k+1L]]]==rank[b$scenario[j]]+1L
      if(left && right && all(b$quadratic_status[ix[c(k-1L,k+1L)]]=="OK") &&
         all(b$winner[ix[c(k-1L,k+1L)]]==b$winner[j])) {
        q <- ix[c(k-1L,k,k+1L)];b$dt_da[j] <- p10_three_point(b$a[q],b$quadratic_position[q],2L)
        b$dt_dloga[j] <- p10_three_point(log(b$a[q]),b$quadratic_position[q],2L)
        b$sensitivity_method[j] <- "CENTRAL_THREE_POINT"
        coarse <- (b$quadratic_position[q[3]]-b$quadratic_position[q[1]])/(b$a[q[3]]-b$a[q[1]])
        b$sensitivity_change_under_refinement[j] <- b$dt_da[j]-coarse
        if(k>=3L && k<=length(ix)-2L) {
          far <- ix[c(k-2L,k,k+2L)]
          if(all(b$quadratic_status[far]=="OK") && all(b$winner[far]==b$winner[j]) &&
             rank[b$scenario[far[1]]]==rank[b$scenario[j]]-2L &&
             rank[b$scenario[far[3]]]==rank[b$scenario[j]]+2L)
            b$nested_sensitivity_change[j] <- b$dt_da[j]-p10_three_point(b$a[far],b$quadratic_position[far],2L)
        }
        anchor <- as.integer(b$position[j])
        slopes <- curves <- rep(NA_real_,3L)
        for(m in 1:3) {
          yy <- records[[b$scenario[q[m]]]]$trajectory$incidence
          if(anchor>1L && anchor<length(yy)) {
            slopes[m] <- (yy[anchor+1L]-yy[anchor-1L])/2
            curves[m] <- yy[anchor-1L]-2*yy[anchor]+yy[anchor+1L]
          }
        }
        if(all(is.finite(slopes)) && all(is.finite(curves)) && abs(curves[2])>s$identity$curvature_tolerance) {
          slope_prime <- p10_three_point(b$a[q],slopes,2L)
          curve_prime <- p10_three_point(b$a[q],curves,2L)
          b$response_numerator[j] <- slope_prime*curves[2]-slopes[2]*curve_prime
          b$response_curvature_denominator[j] <- curves[2]^2
          b$surrogate_dt_da[j] <- -b$response_numerator[j]/b$response_curvature_denominator[j]
        }
      } else if(left || right) {
        q <- if(left) ix[c(k-1L,k)] else ix[c(k,k+1L)]
        if(all(b$quadratic_status[q]=="OK") && all(b$winner[q]==b$winner[j])) {
          b$dt_da[j] <- diff(b$quadratic_position[q])/diff(b$a[q])
          b$dt_dloga[j] <- diff(b$quadratic_position[q])/diff(log(b$a[q]))
          b$sensitivity_method[j] <- if(left) "ONE_SIDED_LEFT" else "ONE_SIDED_RIGHT"
        }
      }
    }
  }
  b$dt_da_units <- "weeks per unit a";b$dt_dloga_units <- "weeks per unit log(a)"
  b
}
p10_run <- function(source,data,s=p10_defaults(),root=".") {
  baseline <- p8_source(source,data,s,root);check <- p10_settings(s,baseline)
  profile <- p8_profile(baseline,s);profile$kind <- "fitted";records <- list()
  for(a in check$a_values) {
    id <- sprintf("scenario-%04d",length(records)+1L)
    r <- p8_evaluate(baseline,check,profile,baseline$parameters$S0,a,id,"a","coarse")
    if(r$status!="VALID") stop("Phase 10 scenario failed: ",r$reason)
    records[[id]] <- r
  }
  base <- records[[which(check$a_values==1)]]$trajectory
  ref <- baseline$evaluation$states[seq_len(check$horizon),,drop=FALSE]
  if(!p6r_close(as.numeric(as.matrix(base[,names(ref)])),as.numeric(as.matrix(ref)),s$replay_tolerance)) stop("Exact baseline replay failed")
  x <- p10_refine(records,baseline,check,profile,s)
  x$branches <- p10_branches(x$records,x$peaks,x$terminal$matches,x$windows,s)
  if(nrow(x$terminal$pairs) && nrow(x$branches)) {
    lookup <- function(sc,yr,id,column) {v <- x$branches[x$branches$scenario==sc & x$branches$year==yr & x$branches$peak_id==id,column];if(length(v)) v[1] else NA}
    for(side in c("left","right")) for(num in 1:2) {
      key <- paste0("branch",num,"_",side,"_peak")
      x$terminal$pairs[[paste0("branch",num,"_",side)]] <- mapply(lookup,x$terminal$pairs[[side]],x$terminal$pairs$year,x$terminal$pairs[[key]],MoreArgs=list(column="branch"))
      x$terminal$pairs[[paste0("branch",num,"_",side,"_quadratic_position")]] <- mapply(lookup,x$terminal$pairs[[side]],x$terminal$pairs$year,x$terminal$pairs[[key]],MoreArgs=list(column="quadratic_position"))
      x$terminal$pairs[[paste0("branch",num,"_",side,"_quadratic_height")]] <- mapply(lookup,x$terminal$pairs[[side]],x$terminal$pairs$year,x$terminal$pairs[[key]],MoreArgs=list(column="quadratic_height"))
    }
  }
  x$type <- "phase10_peak_identity";x$status <- "FINITE_WEEKLY_DIAGNOSTICS_STORED"
  x$diagnostic_settings <- s;x$source_hash <- source$artifact_hash;x$source_identity <- source$identity
  x$source_ancestry <- baseline$source_ancestry;x$baseline_replay <- baseline$replay;x$time_contract <- check
  x
}
p10_report <- function(result,out) {
  if(!identical(result$type,"phase10_peak_identity")) stop("Saved Phase 10 result required")
  write <- function(x,name) utils::write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE)
  scenarios <- do.call(rbind,lapply(result$records,function(r) data.frame(scenario=r$id,a=r$a,s=r$s,stage=r$stage,level=r$level,parent=r$parent,status=r$status)))
  write(scenarios,"scenario_manifest");write(result$peaks,"full_support_local_peaks")
  write(result$windows,"annual_window_maxima");write(result$peaks[,intersect(c("scenario","year","a","peak_id","position","height","quadratic_position","quadratic_height","quadratic_status","curvature","cubic_position","cubic_height","cubic_status"),names(result$peaks)),drop=FALSE],"interpolation_audit")
  write(result$terminal$matches,"forward_reverse_matching");write(result$terminal$candidates,"admissible_match_candidates");write(result$branches,"branch_trajectories_sensitivity")
  write(result$terminal$pairs,"competing_branch_heights");write(result$terminal$edges,"interval_classifications")
  write(result$history,"refinement_history");write(result$decisions,"refinement_decisions");write(result$unresolved,"unresolved_intervals")
  grDevices::pdf(file.path(out,"phase10_diagnostics.pdf"),width=11,height=8)
  tryCatch({for(year in 2018:2019) {
    b <- result$branches[result$branches$year==year,,drop=FALSE]
    w <- result$windows[result$windows$year==year,,drop=FALSE]
    e <- result$terminal$edges[result$terminal$edges$year==year,,drop=FALSE]
    graphics::plot(w$a,w$week,pch=16,xlab="a",ylab="Annual argmax week",main=paste(year,"weekly annual maxima"))
    if(nrow(b)) {
      graphics::plot(b$a,b$week,pch=16,xlab="a",ylab="Weekly local peak week",main=paste(year,"tracked branches"))
      for(id in unique(b$branch)) {v <- b[b$branch==id,];v <- v[order(v$a),];graphics::lines(v$a,v$week,col=2)}
      graphics::plot(b$a,b$height,pch=16,xlab="a",ylab="Weekly peak height",main=paste(year,"branch heights"))
      graphics::plot(b$a,b$quadratic_position,pch=16,xlab="a",ylab="Chronological position",main=paste(year,"weekly and quadratic positions"))
      graphics::points(b$a,b$position,pch=3,col=4)
      graphics::plot(b$a,b$dt_da,pch=16,xlab="a",ylab="weeks per unit a",main=paste(year,"numerical sensitivity"))
      graphics::points(b$a,b$sensitivity_change_under_refinement,pch=3,col=2)
    }
    p <- result$terminal$pairs;p <- p[p$year==year,,drop=FALSE]
    if(nrow(p)) {graphics::plot(p$a_left,p$difference_left,pch=16,xlab="a",ylab="Signed weekly height difference",main=paste(year,"competing branches"));graphics::abline(h=0,lty=2)}
    graphics::plot(e$a_left,seq_len(nrow(e)),col=ifelse(e$classification=="UNRESOLVED",2,1),pch=16,
      xlab="a",ylab="Interval",main=paste(year,"terminal classifications and unresolved regions"))
    graphics::segments(e$a_left,seq_len(nrow(e)),e$a_right,seq_len(nrow(e)),col=ifelse(e$refine,2,1))
    selected <- e[e$refine,,drop=FALSE]
    if(!nrow(selected)) selected <- e[seq_len(min(1L,nrow(e))),,drop=FALSE]
    if(nrow(selected)) for(k in seq_len(min(4L,nrow(selected)))) {
      l <- result$records[[selected$left[k]]]$trajectory;r <- result$records[[selected$right[k]]]$trajectory
      ix <- which(l$year==year);ix <- seq(max(2L,min(ix)-2L),min(nrow(l),max(ix)+2L))
      graphics::plot(ix,l$incidence[ix],type="l",col=2,xlab="Chronological row",ylab="StoI",main=paste(year,"selected supported trajectories"))
      graphics::lines(ix,r$incidence[ix],col=4)
    }
  }},finally=grDevices::dev.off())
  writeLines(c("# Phase 10 peak identity diagnostics","",paste("Accepted source MD5:",result$source_hash),
    paste("Baseline replay:",result$baseline_replay),paste("Scenarios:",length(result$records)),
    paste("Unresolved terminal intervals:",nrow(result$unresolved)),
    "Weekly StoI, full chronological support, and separate annual-window selection are stored.",
    "Branch identity is an operational numerical assignment. Quadratic positions are numerical summaries; shape-preserving cubic positions audit sampled knots.",
    "Refinement and sensitivity are finite-grid diagnostics, not certified convergence or continuous-time SIRS peaks.",
    "No empirical mechanism is asserted by this saved-result summary.","",
    "## Predeclared settings","```",capture.output(dput(result$diagnostic_settings)),"```"),file.path(out,"phase10_summary.md"))
  invisible(result)
}
