# Peak definitions on a fixed noncircular elapsed-time interval.
p4_tolerance <- function(y,policy) policy$tie_abs+policy$tie_rel*max(abs(y))
p4_weekly_peaks <- function(y,time,policy) {
  finite_numeric(y,"weekly outcome",length(time));finite_numeric(time,"elapsed weeks")
  if(any(diff(time)!=1) || any(y<0)) stop("Peak input must be nonnegative and consecutive")
  tol <- p4_tolerance(y,policy);height <- max(y);ties <- time[height-y<=tol]
  rows <- list();i <- 1L
  while(i<=length(y)) {
    j <- i
    # Anchor each plateau at its first value to avoid tolerance-chain drift.
    while(j<length(y) && abs(y[j+1L]-y[i])<=tol) j <- j+1L
    left <- if(i==1L) -Inf else y[i-1L]
    right <- if(j==length(y)) -Inf else y[j+1L]
    if(max(y[i:j])>=max(left,right)+tol || (i==1L && j==length(y)))
      rows[[length(rows)+1L]] <- data.frame(first=time[i],last=time[j],peak_time=time[i],
        height=max(y[i:j]),plateau=j>i,boundary=i==1L || j==length(y),
        isolated=j==i && i>1L && j<length(y))
    i <- j+1L
  }
  local <- if(length(rows)) do.call(rbind,rows) else data.frame(first=numeric(),last=numeric(),
    peak_time=numeric(),height=numeric(),plateau=logical(),boundary=logical(),isolated=logical())
  ranked <- local[order(-local$height,local$peak_time),,drop=FALSE]
  second <- if(nrow(ranked)>=2L) ranked[2L,,drop=FALSE] else NULL
  list(height=height,time=ties[1],ties=ties,tolerance=tol,local=local,
    plateau=length(ties)>1L,boundary=any(ties %in% range(time)),small=height<=policy$small_height,
    runner_up_height=if(is.null(second)) NA_real_ else second$height,
    runner_up_time=if(is.null(second)) NA_real_ else second$peak_time,
    competitor_gap=if(is.null(second)) NA_real_ else abs(ranked$peak_time[1]-second$peak_time))
}
p4_windows <- function(y,time,early,late,policy) {
  tol <- p4_tolerance(y,policy)
  maximum <- function(ids) {
    at <- match(ids,time);h <- max(y[at]);list(height=h,times=ids[h-y[at]<=tol])
  }
  e <- maximum(early);l <- maximum(late);outside <- setdiff(time,c(early,late))
  other <- if(length(outside)) maximum(outside) else list(height=NA_real_,times=numeric())
  list(H_E=e$height,H_L=l$height,H_O=other$height,D=e$height-l$height,
    early_times=e$times,late_times=l$times,outside_times=other$times,
    outside_dominant=length(outside)>0L && other$height>max(e$height,l$height)+tol,
    outside_tied=length(outside)>0L && abs(other$height-max(e$height,l$height))<=tol,
    both_dominate_outside=!length(outside) || min(e$height,l$height)>other$height+tol,
    dominance_margin=if(length(outside)) min(e$height,l$height)-other$height else NA_real_)
}
p4_movement <- function(left,right,policy) {
  if(left$status!="VALID" || right$status!="VALID") return("COMPUTATION_FAILED")
  if(left$small || right$small) return("LOW_OUTCOME_TIMING_UNSTABLE")
  if(left$outside_dominant || right$outside_dominant) return("OUTSIDE_WINDOW_INTERVENTION")
  if(left$boundary || right$boundary) return("BOUNDARY_MAXIMUM")
  if(left$n_ties>1L || right$n_ties>1L) return("PLATEAU_OR_TIE")
  distance <- abs(left$peak_time-right$peak_time)
  if(distance==0) "STABLE_WEEK" else if(distance==1) "ADJACENT_WEEK_EXCHANGE" else
    if(distance>=policy$min_separation) "DISTANT_PEAK_EXCHANGE" else "LOCAL_MOVEMENT_UNRESOLVED"
}
p4_crossing_candidates <- function(table,tolerance) {
  rows <- list()
  for(i in seq_len(nrow(table))) {
    valid <- identical(table$status[i],"VALID") && is.finite(table$D[i])
    if(valid && abs(table$D[i])<=tolerance) {
      left <- if(i>1L && table$status[i-1L]=="VALID") i-1L else i
      right <- if(i<nrow(table) && table$status[i+1L]=="VALID") i+1L else i
      rows[[length(rows)+1L]] <- data.frame(kind="GRID_TIE",left=table$theta[left],right=table$theta[right],
        D_left=table$D[left],D_right=table$D[right],tie_theta=table$theta[i])
    }
    if(i<nrow(table) && valid && table$status[i+1L]=="VALID" && is.finite(table$D[i+1L]) &&
       abs(table$D[i])>tolerance && abs(table$D[i+1L])>tolerance && sign(table$D[i])!=sign(table$D[i+1L]))
      rows[[length(rows)+1L]] <- data.frame(kind="SIGN_CHANGE",left=table$theta[i],right=table$theta[i+1L],
        D_left=table$D[i],D_right=table$D[i+1L],tie_theta=NA_real_)
  }
  if(!length(rows)) return(data.frame(kind=character(),left=numeric(),right=numeric(),
    D_left=numeric(),D_right=numeric(),tie_theta=numeric()))
  result <- do.call(rbind,rows);result$id <- seq_len(nrow(result));result
}
p4_transmission_peaks <- function(source,cfg) {
  b <- source$basis;c <- constrained_coefficients(source$parameters$free_coordinates,b)
  settings <- cfg$transmission;mesh <- sort(unique(c(seq(0,52,by=settings$mesh),b$knots,52)))
  gprime <- function(x) drop(cyclic_basis_matrix(x,b,1L)%*%c)
  derivative <- gprime(mesh);tol <- settings$derivative_tol
  profile <- data.frame(x=mesh,seasonal_week_coordinate=mesh+1,
    beta=evaluate_spline_beta(mesh,c,source$parameters$beta0,b),gprime=derivative,
    beta_prime=evaluate_spline_beta(mesh,c,source$parameters$beta0,b,1L))
  if(max(abs(derivative))<=tol) return(list(status="CONSTANT_OR_NUMERICALLY_FLAT",profile=profile,
    all=data.frame(x=numeric(),beta=numeric(),prominence=numeric()),filtered=data.frame(),threshold_counts=data.frame(),
    limitation="No isolated peaks assigned to a flat curve"))
  roots <- numeric()
  # Sign brackets plus sampled near-zero values; tangential stationary points may be missed.
  for(i in seq_len(length(mesh)-1L)) {
    if(abs(derivative[i])<=tol) roots <- c(roots,mesh[i])
    if(derivative[i]*derivative[i+1L]<0) roots <- c(roots,stats::uniroot(gprime,mesh[i:(i+1L)],tol=settings$root_tol)$root)
  }
  roots <- sort(roots %% 52)
  if(length(roots)) roots <- roots[c(TRUE,diff(roots)>settings$root_tol*10)]
  if(length(roots)>1L && 52-tail(roots,1)+roots[1]<=settings$root_tol*10) roots <- head(roots,-1L)
  probe <- min(settings$mesh/10,b$spacing/100)
  rows <- list()
  for(x in roots) {
    left <- gprime(x-probe);right <- gprime(x+probe)
    if(left>tol && right< -tol) {
      # Prominence proxy uses the higher of the two adjacent circular basin floors.
      other <- roots[roots!=x]
      before <- if(length(other)) max((other-x) %% 52) else 26
      after <- if(length(other)) min((other-x) %% 52) else 26
      left_grid <- seq(x-(52-before),x,length.out=101)
      right_grid <- seq(x,x+after,length.out=101)
      base <- max(min(evaluate_spline_beta(left_grid,c,source$parameters$beta0,b)),
                  min(evaluate_spline_beta(right_grid,c,source$parameters$beta0,b)))
      height <- evaluate_spline_beta(x,c,source$parameters$beta0,b)
      rows[[length(rows)+1L]] <- data.frame(x=x,seasonal_week_coordinate=x+1,beta=height,
        prominence=max(0,height-base),gprime_left=left,gprime_right=right,
        beta_second=evaluate_spline_beta(x,c,source$parameters$beta0,b,2L))
    }
  }
  maxima <- if(length(rows)) do.call(rbind,rows) else data.frame(x=numeric(),beta=numeric(),prominence=numeric())
  filter_peaks <- function(prominence,separation) {
    ids <- which(maxima$prominence>=prominence);ids <- ids[order(-maxima$beta[ids],maxima$x[ids])];keep <- integer()
    for(i in ids) if(!length(keep) || all(pmin(abs(maxima$x[i]-maxima$x[keep]),52-abs(maxima$x[i]-maxima$x[keep]))>=separation)) keep <- c(keep,i)
    maxima[sort(keep),,drop=FALSE]
  }
  thresholds <- settings$thresholds
  if(is.null(thresholds)) thresholds <- data.frame(prominence=settings$prominence,min_separation=settings$min_separation)
  counts <- thresholds;counts$n_peaks <- vapply(seq_len(nrow(thresholds)),function(i)
    nrow(filter_peaks(thresholds$prominence[i],thresholds$min_separation[i])),integer(1))
  list(status="NUMERICAL_PEAKS",profile=profile,all=maxima,
    filtered=filter_peaks(settings$prominence,settings$min_separation),threshold_counts=counts,
    limitation="Sign/near-zero mesh detection; flat segments and tangential roots can be unresolved. Prominence is a sampled adjacent-basin proxy (101 points per side).")
}
