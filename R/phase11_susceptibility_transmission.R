# Counterfactual timing diagnostics; the production SIRS transition is unchanged.
p11_defaults <- function() {
  h <- log(1.2)
  list(source_lambda=1,replay_tolerance=1e-7,grid=seq(-h,h,length.out=9L),
    derivative_steps=h/c(4,8,16),levels=2L,max_scenarios=180L,
    refinement=list(axis=16L,joint=16L,per_year=8L,rapid_weeks=2,
      near_height_fraction=.02),stability_tolerance=.5,
    peak=p10_defaults()$identity,peaks=p10_defaults()$peaks)
}
p11_validate <- function(s,s0,I0) {
  for(k in c("source_lambda","replay_tolerance","stability_tolerance")) scalar_number(s[[k]],k,0,strict_lower=TRUE)
  if(s$source_lambda!=1) stop("Phase 11 requires lambda=1")
  g <- s$grid;finite_numeric(g,"log grid");if(length(g)<3L || any(diff(g)<=0) || !any(g==0) ||
    !isTRUE(all.equal(g,-rev(g),tolerance=1e-12))) stop("Log grid must increase, be symmetric, and contain exact zero")
  st <- s$derivative_steps;finite_numeric(st,"derivative steps")
  if(any(st<=0) || any(diff(st)>=0) || max(st)>max(g)) stop("Positive descending nested derivative steps required")
  for(k in c("levels","max_scenarios")) integer_scalar(s[[k]],k,if(k=="levels") 0L else 1L)
  for(k in c("axis","joint","per_year")) integer_scalar(s$refinement[[k]],k,0L)
  for(k in c("rapid_weeks","near_height_fraction")) scalar_number(s$refinement[[k]],k,0)
  if(any(s0*exp(g)<=0 | s0*exp(g)>=1-I0)) stop("Configured susceptibility range violates valid simplex; provide an explicitly revised Phase 11 configuration")
  needed <- length(g)^2+4L*length(st)
  if(needed>s$max_scenarios) stop("Scenario cap below predeclared grid and derivative cells")
  if(needed+s$refinement$axis+s$refinement$joint>s$max_scenarios)
    stop("Scenario cap below predeclared grid and refinement allocation")
  invisible(s)
}
p11_key <- function(u,v) sprintf("%.14f/%.14f",u,v)
p11_plan <- function(s) {
  g <- s$grid;st <- s$derivative_steps;z <- expand.grid(u=g,v=g);z$stage <- "joint_grid"
  extra <- data.frame(u=c(st,-st,rep(0,2*length(st))),
    v=c(rep(0,2*length(st)),st,-st),stage="local_derivative")
  z <- rbind(z,extra);z$u[abs(z$u)<1e-14] <- 0;z$v[abs(z$v)<1e-14] <- 0
  z <- z[!duplicated(p11_key(z$u,z$v)),,drop=FALSE];z$id <- sprintf("scenario-%04d",seq_len(nrow(z)))
  z$axis <- ifelse(z$u==0 & z$v==0,"baseline",ifelse(z$v==0,"S",ifelse(z$u==0,"beta","joint")))
  z$s_multiplier <- exp(z$u);z$beta_multiplier <- exp(z$v);z$level <- 0L
  rownames(z) <- NULL;z
}
p11_scenario <- function(baseline,check,profile,row) {
  r <- p8_evaluate(baseline,check,profile,baseline$parameters$S0*exp(row$u),exp(row$v),
    row$id,row$axis,row$stage,row$level)
  r$u <- row$u;r$v <- row$v;r
}
p11_outcomes <- function(records,years,s) {
  peaks <- windows <- annual <- list()
  for(r in records) {
    if(r$status!="VALID") next
    z <- r$trajectory
    for(year in years) {
      ix <- which(z$year==year & !z$initialization_only)
      if(!length(ix)) next
      p <- p10_detect(r,year,s$peaks,s$peak)
      if(nrow(p)) {p$u <- r$u;p$v <- r$v;peaks[[length(peaks)+1L]] <- p}
      w <- p10_window(r,year,p,s$peaks);w$u <- r$u;w$v <- r$v
      if(max(z$incidence[ix])<=s$peaks$small_height) w$status <- "LOW_INCIDENCE"
      windows[[length(windows)+1L]] <- w
      annual[[length(annual)+1L]] <- data.frame(scenario=r$id,year=year,u=r$u,v=r$v,
        complete=identical(as.integer(z$week[ix]),1:52),weeks=length(ix),
        cumulative_incidence=sum(z$incidence[ix]),maximum=max(z$incidence[ix]),
        outcome_status=w$status,initialization_row=1L,initialization_year=z$year[1],
        initialization_week=z$week[1],elapsed_weeks_to_window=min(ix)-1L)
    }
  }
  bind <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
  list(peaks=bind(peaks),windows=bind(windows),annual=bind(annual))
}
p11_neighbor_pairs <- function(plan) {
  g <- sort(unique(plan$u[plan$v==0 & plan$stage=="joint_grid"]))
  rows <- list()
  for(axis in c("u","v")) for(fixed in g) {
    q <- plan[plan$stage=="joint_grid" & abs(plan[[if(axis=="u") "v" else "u"]]-fixed)<1e-12,,drop=FALSE]
    q <- q[order(q[[axis]]),,drop=FALSE]
    if(nrow(q)>1L) for(i in seq_len(nrow(q)-1L)) rows[[length(rows)+1L]] <-
      data.frame(axis=axis,left=q$id[i],right=q$id[i+1L],u_left=q$u[i],v_left=q$v[i],u_right=q$u[i+1L],v_right=q$v[i+1L])
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p11_axis_pairs <- function(plan) {
  rows <- list()
  for(axis in c("u","v")) {
    q <- plan[plan[[if(axis=="u") "v" else "u"]]==0,,drop=FALSE]
    q <- q[order(q[[axis]]),,drop=FALSE]
    if(nrow(q)>1L) for(i in seq_len(nrow(q)-1L)) rows[[length(rows)+1L]] <-
      data.frame(axis=axis,left=q$id[i],right=q$id[i+1L],u_left=q$u[i],v_left=q$v[i],u_right=q$u[i+1L],v_right=q$v[i+1L])
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p11_edges <- function(plan,peaks,years,d,pairs=p11_neighbor_pairs(plan)) {
  out <- list()
  for(year in years) for(i in seq_len(nrow(pairs))) {
    x <- pairs[i,];L <- peaks[peaks$scenario==x$left & peaks$year==year,,drop=FALSE]
    R <- peaks[peaks$scenario==x$right & peaks$year==year,,drop=FALSE]
    f <- p10_assignment(L,R,d);v <- p10_assignment(L,R,d,TRUE)
    consistent <- p10_identity_consistent(f,v)
    for(tab in list(f,v)) if(nrow(tab)) {
      tab$year <- year;tab$axis <- x$axis;tab$left <- x$left;tab$right <- x$right
      tab$direction_consistent <- consistent;out[[length(out)+1L]] <- tab
    }
  }
  if(length(out)) do.call(rbind,out) else data.frame()
}
p11_refine_axes <- function(plan,records,baseline,check,profile,years,s) {
  history <- decisions <- list();count <- 0L;per_year <- setNames(integer(length(years)),as.character(years))
  for(level in seq_len(s$levels)) {
    x <- p11_outcomes(records,years,s);pairs <- p11_axis_pairs(plan)
    edge <- p11_edges(plan,x$peaks,years,s$peak,pairs)
    candidate <- list()
    for(year in years) for(i in seq_len(nrow(pairs))) {
      p <- pairs[i,];a <- x$windows[x$windows$scenario==p$left & x$windows$year==year,,drop=FALSE]
      b <- x$windows[x$windows$scenario==p$right & x$windows$year==year,,drop=FALSE]
      if(nrow(a)!=1L || nrow(b)!=1L) next
      e <- edge[edge$year==year & edge$left==p$left & edge$right==p$right &
        edge$direction=="forward" & edge$status=="MATCHED" & edge$direction_consistent,,drop=FALSE]
      same <- a$status=="GENUINE_LOCAL_PEAK" && b$status=="GENUINE_LOCAL_PEAK" &&
        any(e$left_peak==a$peak_id & e$right_peak==b$peak_id)
      L <- x$peaks[x$peaks$scenario==p$left & x$peaks$year==year & x$peaks$genuine,,drop=FALSE]
      R <- x$peaks[x$peaks$scenario==p$right & x$peaks$year==year & x$peaks$genuine,,drop=FALSE]
      gap <- function(P) if(nrow(P)>=2L) {
        y <- sort(P$height,decreasing=TRUE);(y[1]-y[2])/max(y[1],.Machine$double.eps)
      } else Inf
      near <- min(gap(L),gap(R))<=s$refinement$near_height_fraction
      rapid <- same && abs(a$position-b$position)>=s$refinement$rapid_weeks
      clear <- a$status=="GENUINE_LOCAL_PEAK" && b$status=="GENUINE_LOCAL_PEAK"
      unclear <- !clear || !same
      priority <- if(clear && !same) 1L else if(near) 2L else if(rapid) 3L else if(unclear) 4L else 99L
      if(priority==99L) next
      row <- data.frame(level=level,year=year,axis=p$axis,left=p$left,right=p$right,
        log_left=if(p$axis=="u") p$u_left else p$v_left,
        log_right=if(p$axis=="u") p$u_right else p$v_right,
        priority=priority,near_height=near,rapid_displacement=rapid,identity_unresolved=unclear)
      candidate[[length(candidate)+1L]] <- row
    }
    if(!length(candidate)) break
    tab <- do.call(rbind,candidate);tab <- tab[order(tab$priority,tab$year,tab$axis,tab$log_left),,drop=FALSE]
    history[[length(history)+1L]] <- tab;added <- 0L
    for(i in seq_len(nrow(tab))) {
      z <- tab[i,];mid <- (z$log_left+z$log_right)/2
      u <- if(z$axis=="u") mid else 0;v <- if(z$axis=="v") mid else 0
      present <- any(abs(plan$u-u)<1e-12 & abs(plan$v-v)<1e-12)
      decision <- if(present) "SHARED_CELL" else if(per_year[as.character(z$year)]>=s$refinement$per_year) "YEAR_CAP" else
        if(count>=s$refinement$axis || nrow(plan)>=s$max_scenarios) "AXIS_OR_SCENARIO_CAP" else "SELECTED"
      decisions[[length(decisions)+1L]] <- data.frame(z,midpoint=mid,decision=decision)
      if(decision!="SELECTED") next
      row <- data.frame(u=u,v=v,stage="axis_refinement",id=sprintf("scenario-%04d",nrow(plan)+1L),
        axis=if(z$axis=="u") "S" else "beta",s_multiplier=exp(u),beta_multiplier=exp(v),level=level)
      plan <- rbind(plan,row);records[[row$id]] <- p11_scenario(baseline,check,profile,row)
      count <- count+1L;added <- added+1L;per_year[as.character(z$year)] <- per_year[as.character(z$year)]+1L
    }
    if(!added) break
  }
  list(plan=plan,records=records,history=if(length(history)) do.call(rbind,history) else data.frame(),
    decisions=if(length(decisions)) do.call(rbind,decisions) else data.frame(),axis_added=count)
}
p11_follow <- function(edges,year,left,right,peak) {
  e <- edges[edges$year==year & edges$left==left & edges$right==right &
    edges$direction=="forward" & edges$status=="MATCHED" & edges$direction_consistent & edges$left_peak==peak,,drop=FALSE]
  if(nrow(e)==1L) e$right_peak else NA_integer_
}
p11_loops <- function(plan,edges,peaks,years) {
  g <- sort(unique(plan$u[plan$stage=="joint_grid"]));rows <- list()
  id <- function(u,v) plan$id[which(abs(plan$u-u)<1e-12 & abs(plan$v-v)<1e-12)[1]]
  for(year in years) for(i in seq_len(length(g)-1L)) for(j in seq_len(length(g)-1L)) {
    a <- id(g[i],g[j]);b <- id(g[i+1L],g[j]);c <- id(g[i],g[j+1L]);d <- id(g[i+1L],g[j+1L])
    P <- peaks[peaks$scenario==a & peaks$year==year & peaks$genuine,,drop=FALSE]
    for(p in P$peak_id) {
      ab <- p11_follow(edges,year,a,b,p);ac <- p11_follow(edges,year,a,c,p)
      bd <- if(is.na(ab)) NA_integer_ else p11_follow(edges,year,b,d,ab)
      cd <- if(is.na(ac)) NA_integer_ else p11_follow(edges,year,c,d,ac)
      rows[[length(rows)+1L]] <- data.frame(year=year,u=g[i],v=g[j],start=a,start_peak=p,
        u_then_v=bd,v_then_u=cd,status=if(is.na(bd) || is.na(cd)) "UNRESOLVED_EDGE" else
          if(bd==cd) "CONSISTENT" else "INCONSISTENT_PATH")
    }
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p11_regions <- function(plan,peaks,windows,edges,axis_edges,loops) {
  if(!nrow(peaks)) return(list(branches=data.frame(),winners=data.frame()))
  nodes <- peaks[,c("year","scenario","peak_id","position","week","height","quadratic_position","quadratic_status"),drop=FALSE]
  nodes$key <- paste(nodes$year,nodes$scenario,nodes$peak_id,sep="/")
  if(anyDuplicated(nodes$key)) stop("Duplicate local peak key")
  adjacency <- vector("list",nrow(nodes))
  e <- rbind(edges,axis_edges)
  if(nrow(e)) {
    e <- e[e$direction=="forward" & e$status=="MATCHED" & e$direction_consistent,,drop=FALSE]
    for(i in seq_len(nrow(e))) {
      a <- match(paste(e$year[i],e$left[i],e$left_peak[i],sep="/"),nodes$key)
      b <- match(paste(e$year[i],e$right[i],e$right_peak[i],sep="/"),nodes$key)
      if(is.na(a) || is.na(b)) next
      adjacency[[a]] <- c(adjacency[[a]],b);adjacency[[b]] <- c(adjacency[[b]],a)
    }
  }
  component <- rep(NA_integer_,nrow(nodes));ncomp <- 0L
  for(start in seq_len(nrow(nodes))) if(is.na(component[start])) {
    ncomp <- ncomp+1L;queue <- start;component[start] <- ncomp;head <- 1L
    while(head<=length(queue)) {
      a <- queue[head];head <- head+1L
      for(b in adjacency[[a]]) if(is.na(component[b])) {component[b] <- ncomp;queue <- c(queue,b)}
    }
  }
  nodes$branch <- paste0(nodes$year,"/joint-branch-",component)
  bad <- rep(FALSE,nrow(nodes))
  if(nrow(loops)) for(i in which(loops$status!="CONSISTENT")) {
    u <- loops$u[i];v <- loops$v[i]
    g <- sort(unique(plan$u[plan$stage=="joint_grid"]));du <- min(diff(g));dv <- du
    corners <- plan$id[plan$stage=="joint_grid" &
      ((abs(plan$u-u)<1e-12 | abs(plan$u-(u+du))<1e-12) &
       (abs(plan$v-v)<1e-12 | abs(plan$v-(v+dv))<1e-12))]
    bad <- bad | (nodes$year==loops$year[i] & nodes$scenario %in% corners)
  }
  for(k in unique(component)) {
    ix <- which(component==k)
    if(anyDuplicated(paste(nodes$year[ix],nodes$scenario[ix],sep="/"))) bad[ix] <- TRUE
    if(any(bad[ix])) bad[ix] <- TRUE
  }
  nodes$identity_resolved <- !bad & peaks$genuine
  w <- windows[,c("year","scenario","peak_id","position","week","height","status"),drop=FALSE]
  w$key <- paste(w$year,w$scenario,w$peak_id,sep="/")
  at <- match(w$key,nodes$key);w$branch <- nodes$branch[at]
  w$identity_resolved <- !is.na(at) & nodes$identity_resolved[at]
  w$identity_resolved[is.na(w$identity_resolved)] <- FALSE
  list(branches=nodes,winners=w)
}
p11_derivative <- function(x,y,at=2L) p10_three_point(x,y,at)
p11_response <- function(plan,windows,peaks,edges,years,s) {
  base <- plan$id[plan$u==0 & plan$v==0];rows <- list();der <- list()
  find <- function(u,v) {i <- which(abs(plan$u-u)<1e-12 & abs(plan$v-v)<1e-12);if(length(i)) plan$id[i[1]] else NA_character_}
  win <- function(id,year) windows[windows$scenario==id & windows$year==year,,drop=FALSE]
  pos <- function(id,year,peak) {p <- peaks[peaks$scenario==id & peaks$year==year & peaks$peak_id==peak,,drop=FALSE];if(nrow(p) && p$quadratic_status[1]=="OK") p$quadratic_position[1] else NA_real_}
  for(year in years) for(axis in c("u","v")) {
    w0 <- win(base,year)
    q <- plan[plan$axis %in% c(if(axis=="u") "S" else "beta", "baseline"),,drop=FALSE]
    q <- q[order(q[[axis]]),,drop=FALSE]
    for(i in seq_len(nrow(q))) {
      w <- win(q$id[i],year);if(!nrow(w) || !nrow(w0)) next
      matched <- if(q$id[i]==base) TRUE else {
        chain <- q$id[seq.int(which(q$id==base),i)]
        p <- w0$peak_id;ok <- !is.na(p) && w0$status=="GENUINE_LOCAL_PEAK"
        if(length(chain)>1L) for(k in seq_len(length(chain)-1L)) {
          forward <- i>which(q$id==base)
          e <- edges[edges$year==year & edges$left==(if(forward) chain[k] else chain[k+1L]) &
            edges$right==(if(forward) chain[k+1L] else chain[k]) & edges$direction=="forward" &
            edges$status=="MATCHED" & edges$direction_consistent,,drop=FALSE]
          if(!nrow(e)) {ok <- FALSE;break}
          nextp <- if(forward) e$right_peak[e$left_peak==p] else e$left_peak[e$right_peak==p]
          if(length(nextp)!=1L) {ok <- FALSE;break};p <- nextp
        }
        ok && identical(as.integer(p),as.integer(w$peak_id))
      }
      status <- if(w$status!="GENUINE_LOCAL_PEAK" || w0$status!="GENUINE_LOCAL_PEAK") "UNDEFINED_TIMING" else
        if(matched) "PERSISTENT_WINNER" else "SWITCH_OR_UNRESOLVED"
      rows[[length(rows)+1L]] <- data.frame(year=year,axis=axis,log_change=q[[axis]][i],multiplier=exp(q[[axis]][i]),
        scenario=q$id[i],baseline=base,delta_week=w$position-w0$position,absolute_delta=abs(w$position-w0$position),
        baseline_peak=w0$peak_id,scenario_peak=w$peak_id,mechanism=status,
        timing_status=w$status,baseline_status=w0$status)
    }
    for(h in s$derivative_steps) {
      ids <- c(find(if(axis=="u") -h else 0,if(axis=="v") -h else 0),base,
        find(if(axis=="u") h else 0,if(axis=="v") h else 0))
      ww <- lapply(ids,win,year=year)
      valid <- all(vapply(ww,function(w) nrow(w)==1L && w$status=="GENUINE_LOCAL_PEAK",logical(1)))
      pp <- if(valid) vapply(seq_along(ids),function(i) pos(ids[i],year,ww[[i]]$peak_id),numeric(1)) else rep(NA_real_,3)
      if(valid) {
        L <- peaks[peaks$scenario==ids[1] & peaks$year==year,,drop=FALSE]
        C <- peaks[peaks$scenario==ids[2] & peaks$year==year,,drop=FALSE]
        R <- peaks[peaks$scenario==ids[3] & peaks$year==year,,drop=FALSE]
        left <- p10_assignment(L,C,s$peak);right <- p10_assignment(C,R,s$peak)
        valid <- p10_identity_consistent(left,p10_assignment(L,C,s$peak,TRUE)) &&
          p10_identity_consistent(right,p10_assignment(C,R,s$peak,TRUE)) &&
          any(left$left_peak==ww[[1]]$peak_id & left$right_peak==ww[[2]]$peak_id & left$status=="MATCHED") &&
          any(right$left_peak==ww[[2]]$peak_id & right$right_peak==ww[[3]]$peak_id & right$status=="MATCHED") && all(is.finite(pp))
      }
      der[[length(der)+1L]] <- data.frame(year=year,axis=axis,h=h,method="CENTRAL_THREE_POINT",
        sensitivity=if(valid) p11_derivative(c(-h,0,h),pp) else NA_real_,
        valid=valid,reason=if(valid) "PERSISTENT_INTERPOLATED_BRANCH" else "IDENTITY_TIE_OR_INTERPOLATION_INADEQUATE")
    }
  }
  response <- if(length(rows)) do.call(rbind,rows) else data.frame()
  derivative <- if(length(der)) do.call(rbind,der) else data.frame()
  if(nrow(derivative)) {derivative$nested_change <- NA_real_;derivative$stable <- FALSE
    for(year in years) for(axis in c("u","v")) {
      ix <- which(derivative$year==year & derivative$axis==axis);ix <- ix[order(-derivative$h[ix])]
      if(length(ix)>1L) for(k in seq_len(length(ix)-1L)) {
        derivative$nested_change[ix[k]] <- abs(derivative$sensitivity[ix[k]]-derivative$sensitivity[ix[k+1L]])
        derivative$stable[ix[k]] <- is.finite(derivative$nested_change[ix[k]]) && derivative$nested_change[ix[k]]<=s$stability_tolerance
      }
    }
  }
  list(response=response,derivative=derivative)
}
p11_interactions <- function(plan,windows,peaks,edges,loops,years,d) {
  g <- sort(unique(plan$u[plan$stage=="joint_grid"]));out <- list();global <- list()
  id <- function(u,v) plan$id[which(abs(plan$u-u)<1e-12 & abs(plan$v-v)<1e-12)[1]]
  for(year in years) for(u in g) for(v in g) {
    ids <- c(id(0,0),id(u,0),id(0,v),id(u,v))
    ww <- lapply(ids,function(a) windows[windows$scenario==a & windows$year==year,,drop=FALSE])
    if(any(vapply(ww,nrow,integer(1))!=1L)) next
    raw <- vapply(ww,function(w) w$position,numeric(1));delta <- raw[4]-raw[2]-raw[3]+raw[1]
    status <- vapply(ww,function(w) w$status,character(1))
    global[[length(global)+1L]] <- data.frame(year=year,u=u,v=v,T00=raw[1],Tu0=raw[2],T0v=raw[3],Tuv=raw[4],
      delta_global_week=delta,valid_timing=all(status=="GENUINE_LOCAL_PEAK"),
      mechanism=if(all(status=="GENUINE_LOCAL_PEAK")) "DESCRIPTIVE_WINNER_CONTRAST" else "UNDEFINED_TIMING")
    p0 <- ww[[1]]$peak_id;pu <- ww[[2]]$peak_id;pv <- ww[[3]]$peak_id;p11 <- ww[[4]]$peak_id
    matching <- all(status=="GENUINE_LOCAL_PEAK")
    if(matching) {
      nodes <- list(c(0,0),c(u,0),c(0,v),c(u,v));want <- c(p0,pu,pv,p11)
      for(k in 2:4) {
        L <- peaks[peaks$scenario==ids[1] & peaks$year==year,,drop=FALSE]
        R <- peaks[peaks$scenario==ids[k] & peaks$year==year,,drop=FALSE]
        f <- p10_assignment(L,R,d);rev <- p10_assignment(L,R,d,TRUE)
        if(!p10_identity_consistent(f,rev) ||
          !any(f$left_peak==p0 & f$right_peak==want[k] & f$status=="MATCHED")) matching <- FALSE
      }
      if(u!=0 && v!=0 && nrow(loops)) {
        in_box <- loops$year==year & loops$u>=min(0,u)-1e-12 & loops$u<max(0,u)-1e-12 &
          loops$v>=min(0,v)-1e-12 & loops$v<max(0,v)-1e-12
        if(any(loops$status[in_box]!="CONSISTENT")) matching <- FALSE
      }
    }
    T <- vapply(seq_along(ids),function(k) {
      p <- peaks[peaks$scenario==ids[k] & peaks$year==year & peaks$peak_id==want[k],,drop=FALSE]
      if(nrow(p) && p$quadratic_status[1]=="OK") p$quadratic_position[1] else NA_real_
    },numeric(1))
    matching <- matching && all(is.finite(T))
    out[[length(out)+1L]] <- data.frame(year=year,u=u,v=v,T00=T[1],Tu0=T[2],T0v=T[3],Tuv=T[4],
      branch00=p0,branch_u0=pu,branch_0v=pv,branch_uv=p11,
      delta_interaction=if(matching) T[4]-T[2]-T[3]+T[1] else NA_real_,
      valid=matching,reason=if(matching) "CONSISTENT_INTERPOLATED_BRANCH" else "IDENTITY_OR_TIMING_UNRESOLVED")
  }
  list(branch=if(length(out)) do.call(rbind,out) else data.frame(),global=if(length(global)) do.call(rbind,global) else data.frame())
}
p11_comparisons <- function(derivatives,years) {
  rows <- list()
  for(year in years) for(h in sort(unique(derivatives$h))) {
    a <- derivatives[derivatives$year==year & derivatives$axis=="u" & derivatives$h==h,,drop=FALSE]
    b <- derivatives[derivatives$year==year & derivatives$axis=="v" & derivatives$h==h,,drop=FALSE]
    if(nrow(a)!=1L || nrow(b)!=1L) next
    good <- is.finite(a$sensitivity) && is.finite(b$sensitivity)
    rows[[length(rows)+1L]] <- data.frame(year=year,u=0,v=0,h=h,
      E_S=a$sensitivity,E_beta=b$sensitivity,
      absolute_E_S=abs(a$sensitivity),absolute_E_beta=abs(b$sensitivity),
      absolute_difference=abs(a$sensitivity)-abs(b$sensitivity),
      opposite_signs=if(good) sign(a$sensitivity)!=sign(b$sensitivity) else NA,
      valid=good,stable=good && a$stable && b$stable,
      units="weeks per unit log change")
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p11_joint_derivatives <- function(plan,windows,peaks,edges,loops,years) {
  g <- sort(unique(plan$u[plan$stage=="joint_grid"]));rows <- list()
  id <- function(u,v) plan$id[which(plan$stage=="joint_grid" & abs(plan$u-u)<1e-12 & abs(plan$v-v)<1e-12)[1]]
  for(year in years) for(u in g) for(v in g) for(axis in c("u","v")) {
    j <- match(if(axis=="u") u else v,g)
    if(j==1L || j==length(g)) next
    coords <- g[(j-1L):(j+1L)]
    ids <- if(axis=="u") vapply(coords,function(a) id(a,v),character(1)) else
      vapply(coords,function(a) id(u,a),character(1))
    W <- lapply(ids,function(a) windows[windows$scenario==a & windows$year==year,,drop=FALSE])
    valid <- all(vapply(W,function(w) nrow(w)==1L && w$status=="GENUINE_LOCAL_PEAK",logical(1)))
    if(valid) for(k in 1:2) {
      e <- edges[edges$year==year & edges$left==ids[k] & edges$right==ids[k+1L] &
        edges$direction=="forward" & edges$status=="MATCHED" & edges$direction_consistent,,drop=FALSE]
      if(!any(e$left_peak==W[[k]]$peak_id & e$right_peak==W[[k+1L]]$peak_id)) valid <- FALSE
    }
    if(valid && nrow(loops)) {
      nearby <- loops$year==year & loops$u>=u-min(diff(g))-1e-12 & loops$u<=u+1e-12 &
        loops$v>=v-min(diff(g))-1e-12 & loops$v<=v+1e-12
      if(any(loops$status[nearby]!="CONSISTENT")) valid <- FALSE
    }
    T <- rep(NA_real_,3)
    if(valid) for(k in 1:3) {
      p <- peaks[peaks$scenario==ids[k] & peaks$year==year & peaks$peak_id==W[[k]]$peak_id,,drop=FALSE]
      if(nrow(p)==1L && p$quadratic_status=="OK") T[k] <- p$quadratic_position
    }
    valid <- valid && all(is.finite(T))
    rows[[length(rows)+1L]] <- data.frame(year=year,u=u,v=v,axis=axis,
      method="CENTRAL_THREE_POINT_UNEQUAL_GRID",sensitivity=if(valid) p11_derivative(coords,T) else NA_real_,
      valid=valid,reason=if(valid) "PERSISTENT_INTERPOLATED_WINNER" else "IDENTITY_OR_NUMERICAL_SUPPORT_INADEQUATE",
      nested_stability="NOT_ASSESSED_AT_JOINT_GRID")
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p11_joint_comparisons <- function(d) {
  if(!nrow(d)) return(data.frame())
  key <- unique(d[,c("year","u","v")]);rows <- list()
  for(i in seq_len(nrow(key))) {
    k <- key[i,];a <- d[d$year==k$year & d$u==k$u & d$v==k$v & d$axis=="u",,drop=FALSE]
    b <- d[d$year==k$year & d$u==k$u & d$v==k$v & d$axis=="v",,drop=FALSE]
    if(nrow(a)!=1L || nrow(b)!=1L) next
    good <- is.finite(a$sensitivity) && is.finite(b$sensitivity)
    rows[[length(rows)+1L]] <- data.frame(year=k$year,u=k$u,v=k$v,E_S=a$sensitivity,E_beta=b$sensitivity,
      absolute_E_S=abs(a$sensitivity),absolute_E_beta=abs(b$sensitivity),
      absolute_difference=abs(a$sensitivity)-abs(b$sensitivity),
      opposite_signs=if(good) sign(a$sensitivity)!=sign(b$sensitivity) else NA,
      valid=good,nested_stability="NOT_ASSESSED_AT_JOINT_GRID")
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p11_thresholds <- function(plan,windows,peaks,edges,years) {
  rows <- list()
  for(year in years) for(axis in c("u","v")) {
    q <- plan[plan[[if(axis=="u") "v" else "u"]]==0,,drop=FALSE]
    q <- q[order(q[[axis]]),,drop=FALSE];found <- FALSE;uncertain <- FALSE
    for(i in seq_len(nrow(q)-1L)) {
      a <- q$id[i];b <- q$id[i+1L]
      wa <- windows[windows$scenario==a & windows$year==year,,drop=FALSE]
      wb <- windows[windows$scenario==b & windows$year==year,,drop=FALSE]
      e <- edges[edges$year==year & edges$left==a & edges$right==b & edges$direction=="forward" &
        edges$status=="MATCHED" & edges$direction_consistent,,drop=FALSE]
      if(nrow(wa)!=1L || nrow(wb)!=1L || wa$status!="GENUINE_LOCAL_PEAK" || wb$status!="GENUINE_LOCAL_PEAK" || nrow(e)<2L) {
        uncertain <- TRUE;next
      }
      same <- any(e$left_peak==wa$peak_id & e$right_peak==wb$peak_id)
      if(same) next
      L <- peaks[peaks$scenario==a & peaks$year==year,,drop=FALSE]
      R <- peaks[peaks$scenario==b & peaks$year==year,,drop=FALSE]
      competitor <- e$left_peak[e$right_peak==wb$peak_id]
      target <- e$right_peak[e$left_peak==wa$peak_id]
      if(length(competitor)!=1L || length(target)!=1L) next
      d1 <- L$height[L$peak_id==wa$peak_id]-L$height[L$peak_id==competitor]
      d2 <- R$height[R$peak_id==target]-R$height[R$peak_id==wb$peak_id]
      if(length(d1)!=1L || length(d2)!=1L || d1<=0 || d2>=0) next
      found <- TRUE
      rows[[length(rows)+1L]] <- data.frame(year=year,axis=axis,log_left=q[[axis]][i],log_right=q[[axis]][i+1L],
        multiplier_left=exp(q[[axis]][i]),multiplier_right=exp(q[[axis]][i+1L]),
        left_branch=wa$peak_id,right_branch=wb$peak_id,height_difference_left=d1,
        height_difference_right=d2,timing_jump_weeks=wb$position-wa$position,
        status="SAMPLED_CONFIDENT_SWITCH_BRACKET")
    }
    if(!found) rows[[length(rows)+1L]] <- data.frame(year=year,axis=axis,log_left=NA_real_,log_right=NA_real_,
      multiplier_left=NA_real_,multiplier_right=NA_real_,left_branch=NA_integer_,right_branch=NA_integer_,
      height_difference_left=NA_real_,height_difference_right=NA_real_,timing_jump_weeks=NA_real_,
      status=if(uncertain) "UNRESOLVED_WITHIN_CONFIGURED_RANGE" else "NOT_DETECTED_WITHIN_CONFIGURED_RANGE")
  }
  do.call(rbind,rows)
}
p11_run <- function(source,data,s=p11_defaults(),root=".") {
  baseline <- p8_source(source,data,s,root)
  p11_validate(s,baseline$parameters$S0,baseline$parameters$I0)
  contract <- p10_defaults();contract$time$years <- 2018:2019
  contract$interventions$s_values <- baseline$parameters$S0*exp(s$grid)
  contract$interventions$a_values <- exp(s$grid)
  contract$interventions$max_scenarios <- s$max_scenarios
  check <- p8_settings(contract,baseline)
  if(!all(check$complete)) stop("Complete 2018 and 2019 reporting windows required")
  years <- intersect(2014:2019,check$years)
  markers <- p3_markers(baseline$data)[seq_len(check$horizon),,drop=FALSE]
  years <- sort(unique(c(2018:2019,vapply(2014:2019,function(y)
    if(identical(as.integer(markers$week[markers$year==y & seq_len(nrow(markers))>1]),1:52)) y else NA_integer_,integer(1)))))
  years <- years[!is.na(years)]
  profile <- p8_profile(baseline,contract);profile$kind <- "fitted"
  plan <- p11_plan(s);records <- list()
  for(i in seq_len(nrow(plan))) records[[plan$id[i]]] <- p11_scenario(baseline,check,profile,plan[i,])
  base_id <- plan$id[plan$u==0 & plan$v==0];base <- records[[base_id]]
  if(base$status!="VALID") stop("Baseline scenario failed: ",base$reason)
  ref <- baseline$evaluation$states[seq_len(check$horizon),,drop=FALSE]
  if(!p6r_close(as.numeric(as.matrix(base$trajectory[,names(ref)])),as.numeric(as.matrix(ref)),s$replay_tolerance)) stop("Exact baseline replay failed")
  refined <- p11_refine_axes(plan,records,baseline,check,profile,years,s)
  plan <- refined$plan;records <- refined$records
  x <- p11_outcomes(records,years,s);edges <- p11_edges(plan,x$peaks,years,s$peak)
  axis_edges <- p11_edges(plan,x$peaks,years,s$peak,p11_axis_pairs(plan))
  loops <- p11_loops(plan,edges,x$peaks,years)
  regions <- p11_regions(plan,x$peaks,x$windows,edges,axis_edges,loops)
  resp <- p11_response(plan,x$windows,x$peaks,axis_edges,years,s)
  inter <- p11_interactions(plan,x$windows,x$peaks,edges,loops,years,s$peak)
  comparisons <- p11_comparisons(resp$derivative,years)
  joint_derivatives <- p11_joint_derivatives(plan,x$windows,x$peaks,edges,loops,years)
  joint_comparisons <- p11_joint_comparisons(joint_derivatives)
  thresholds <- p11_thresholds(plan,x$windows,x$peaks,axis_edges,years)
  list(type="phase11_susceptibility_transmission",status="FINITE_DIAGNOSTICS_STORED",
    source_hash=source$artifact_hash,source_identity=source$identity,baseline_replay=baseline$replay,
    settings=s,plan=plan,records=records,years=years,peaks=x$peaks,windows=x$windows,annual=x$annual,
    edges=edges,axis_edges=axis_edges,loops=loops,branches=regions$branches,winner_regions=regions$winners,
    responses=resp$response,derivatives=resp$derivative,
    joint_derivatives=joint_derivatives,comparisons=comparisons,joint_comparisons=joint_comparisons,
    interactions=inter$branch,global_contrasts=inter$global,thresholds=thresholds,
    unresolved=loops[loops$status!="CONSISTENT",,drop=FALSE],refinement_history=refined$history,
    refinement_decisions=refined$decisions)
}
p11_report <- function(r,out) {
  if(!identical(r$type,"phase11_susceptibility_transmission")) stop("Saved Phase 11 result required")
  items <- list(scenario_manifest=r$plan,annual_outcomes=r$annual,local_peaks=r$peaks,
    annual_window_maxima=r$windows,neighbor_matching=r$edges,axis_neighbor_matching=r$axis_edges,
    grid_loop_consistency=r$loops,local_branch_regions=r$branches,global_winner_regions=r$winner_regions,
    finite_responses=r$responses,s_only_responses=r$responses[r$responses$axis=="u",,drop=FALSE],
    beta_only_responses=r$responses[r$responses$axis=="v",,drop=FALSE],
    local_sensitivities=r$derivatives,joint_grid_sensitivities=r$joint_derivatives,
    equal_change_comparisons=r$comparisons,joint_grid_comparisons=r$joint_comparisons,
    branch_interactions=r$interactions,
    global_timing_contrasts=r$global_contrasts,switching_thresholds=r$thresholds,
    refinement_history=r$refinement_history,refinement_decisions=r$refinement_decisions,unresolved_regions=r$unresolved)
  for(n in names(items)) utils::write.csv(items[[n]],file.path(out,paste0(n,".csv")),row.names=FALSE)
  grDevices::pdf(file.path(out,"phase11_diagnostics.pdf"),width=11,height=8)
  tryCatch({for(year in r$years) {
    x <- r$responses[r$responses$year==year,];d <- r$derivatives[r$derivatives$year==year,]
    graphics::plot(x$log_change,x$delta_week,col=ifelse(x$axis=="u",2,4),pch=16,
      xlab="log multiplier",ylab="Exact weekly timing change",main=paste(year,"axis responses"));graphics::abline(h=0,lty=2)
    graphics::plot(d$h,d$sensitivity,col=ifelse(d$axis=="u",2,4),pch=16,
      xlab="log step",ylab="Weeks per unit log change",main=paste(year,"local interpolated sensitivities"))
    b <- r$peaks[r$peaks$year==year & r$peaks$v==0,,drop=FALSE]
    if(nrow(b)) graphics::plot(b$u,b$height,pch=16,col=ifelse(b$retained,1,2),
      xlab="log susceptibility multiplier",ylab="Weekly local peak height",main=paste(year,"S-axis branch heights"))
    b <- r$peaks[r$peaks$year==year & r$peaks$u==0,,drop=FALSE]
    if(nrow(b)) graphics::plot(b$v,b$height,pch=16,col=ifelse(b$retained,1,2),
      xlab="log transmission multiplier",ylab="Weekly local peak height",main=paste(year,"beta-axis branch heights"))
    w <- r$windows[r$windows$year==year,];graphics::plot(w$u,w$v,col=as.integer(factor(w$position)),pch=16,
      xlab="log susceptibility multiplier",ylab="log transmission multiplier",main=paste(year,"computed winner timing"))
    z <- r$interactions[r$interactions$year==year & r$interactions$valid,]
    if(nrow(z)) graphics::plot(z$u,z$v,col=as.integer(cut(z$delta_interaction,breaks=8)),pch=16,
      xlab="u",ylab="v",main=paste(year,"valid branch interaction"))
    z <- r$loops[r$loops$year==year,]
    if(nrow(z)) graphics::plot(z$u,z$v,col=ifelse(z$status=="CONSISTENT",1,2),pch=16,
      xlab="u",ylab="v",main=paste(year,"loop identity audit"))
    z <- r$thresholds[r$thresholds$year==year & r$thresholds$status=="SAMPLED_CONFIDENT_SWITCH_BRACKET",]
    if(nrow(z)) {
      graphics::plot(z$log_left,seq_len(nrow(z)),pch=16,col=ifelse(z$axis=="u",2,4),
        xlab="log change",ylab="switch bracket",main=paste(year,"sampled switching brackets"))
      graphics::segments(z$log_left,seq_len(nrow(z)),z$log_right,seq_len(nrow(z)),col=ifelse(z$axis=="u",2,4))
    }
  }},finally=grDevices::dev.off())
  writeLines(c("# Phase 11 saved diagnostics","",paste("Source MD5:",r$source_hash),
    paste("Baseline replay:",r$baseline_replay),paste("Scenarios:",nrow(r$plan)),
    paste("Failed scenarios:",sum(vapply(r$records,function(x) x$status!="VALID",logical(1)))),
    paste("Unresolved grid paths:",nrow(r$unresolved)),
    "Sensitivities are weeks per unit log change. Positive is later; negative is earlier.",
    "A 1% relative increase has approximate timing response E * log(1.01), where the local stencil is valid.",
    "These are model-based conditional timing semi-elasticities, not historical causal shares.",
    "Sampled switching brackets and numerical stability do not certify continuous thresholds or convergence."),
    file.path(out,"phase11_summary.md"))
  invisible(r)
}
