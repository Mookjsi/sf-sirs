# Phase 10C reads Phase 10 evidence without recomputing diagnostics.
p10c_defaults <- function() list(high_sensitivity=c(.98,1.01),max_new_scenarios=24L,
  max_per_year=12L,levels=2L)

p10c_validate <- function(s) {
  if(!is.numeric(s$high_sensitivity) || length(s$high_sensitivity)!=2L ||
     any(!is.finite(s$high_sensitivity)) || s$high_sensitivity[1]>=s$high_sensitivity[2])
    stop("Invalid high-sensitivity interval")
  for(n in c("max_new_scenarios","max_per_year","levels")) integer_scalar(s[[n]],n,0L)
  s
}

p10c_audit <- function(result,s=p10c_defaults()) {
  if(!identical(result$type,"phase10_peak_identity")) stop("Saved Phase 10 result required")
  s <- p10c_validate(s)
  e <- result$terminal$edges;u <- result$unresolved
  b <- result$branches;w <- result$windows;p <- result$terminal$pairs
  if(!is.data.frame(e) || !is.data.frame(u)) stop("Saved interval evidence required")
  rows <- vector("list",nrow(u))
  for(i in seq_len(nrow(u))) {
    x <- u[i,];edge <- e[e$year==x$year & abs(e$a_left-x$a_left)<1e-12 &
      abs(e$a_right-x$a_right)<1e-12,,drop=FALSE]
    if(nrow(edge)!=1L) stop("Unresolved interval has no unique terminal evidence")
    q <- p[p$year==x$year & abs(p$a_left-x$a_left)<1e-12 & abs(p$a_right-x$a_right)<1e-12,,drop=FALSE]
    wl <- w[w$scenario==x$left & w$year==x$year,,drop=FALSE]
    wr <- w[w$scenario==x$right & w$year==x$year,,drop=FALSE]
    if(nrow(wl)!=1L || nrow(wr)!=1L) stop("Missing annual-window evidence")
    m <- result$terminal$matches
    m <- m[m$year==x$year & m$left==x$left & m$right==x$right,,drop=FALSE]
    winner <- isTRUE(x$flag_winner_identity_ambiguity) || isTRUE(x$flag_winner_exchange) ||
      isTRUE(x$flag_rapid_displacement) || x$weekly_argmax_movement>0
    competitor <- nrow(q)>0L && any(q$weekly_order_reversal | q$matching_uncertain |
      abs(q$difference_left)/pmax(q$branch1_left_height,q$branch2_left_height,1e-300)<result$diagnostic_settings$identity$near_height_fraction |
      abs(q$difference_right)/pmax(q$branch1_right_height,q$branch2_right_height,1e-300)<result$diagnostic_settings$identity$near_height_fraction)
    high <- x$a_left<=s$high_sensitivity[2] && x$a_right>=s$high_sensitivity[1]
    consequential <- winner || competitor || (isTRUE(x$flag_raw_topology_change) && high)
    priority <- if(winner && competitor) 1L else if(winner) 2L else if(competitor && high) 3L else if(isTRUE(x$flag_raw_topology_change) && high) 4L else 9L
    rows[[i]] <- data.frame(year=x$year,a_left=x$a_left,a_right=x$a_right,width=x$a_right-x$a_left,
      original_classification=x$classification,unresolved_reason=if("reason" %in% names(x)) x$reason else "UNSPECIFIED",
      global_winner_involved=winner,competing_branches_involved=competitor,
      forward_reverse_consistent=isTRUE(x$identity_consistent),
      forward_matches=sum(m$direction=="forward" & m$status=="MATCHED"),
      reverse_matches=sum(m$direction=="reverse" & m$status=="MATCHED"),
      left_window_status=wl$status,right_window_status=wr$status,
      left_peak_support=if(wl$status=="GENUINE_LOCAL_PEAK") "SUPPORTED" else wl$status,
      right_peak_support=if(wr$status=="GENUINE_LOCAL_PEAK") "SUPPORTED" else wr$status,
      raw_peak_creation_loss=isTRUE(x$flag_raw_topology_change),
      high_sensitivity_intersection=high,scientifically_consequential=consequential,
      priority=priority,eligible=consequential && high,stringsAsFactors=FALSE)
  }
  audit <- if(length(rows)) do.call(rbind,rows) else data.frame()
  evidence <- lapply(sort(unique(e$year)),function(year) {
    z <- e[e$year==year,,drop=FALSE]; ww <- w[w$year==year,,drop=FALSE];bb <- b[b$year==year,,drop=FALSE]
    pp <- p[p$year==year,,drop=FALSE];aa <- audit[audit$year==year,,drop=FALSE]
    confirmed <- any(z$classification=="CONFIRMED_WEEKLY_SWITCH" & z$flag_winner_exchange & !z$matching_uncertain)
    displacement <- any(z$classification=="PERSISTENT_BRANCH_DISPLACEMENT" & z$same_winner_branch & !z$matching_uncertain)
    topology <- any(z$flag_raw_topology_change);window <- any(z$flag_window_selection)
    unresolved <- nrow(aa)>0L && any(aa$scientifically_consequential)
    mechanisms <- c(switch=confirmed,displacement=displacement,sampled_creation_loss=topology,window_selection=window)
    reversals <- if(nrow(pp)) sum(pp$weekly_order_reversal) else 0L
    exchanges <- if(nrow(pp)) sum(pp$global_exchange & !pp$matching_uncertain) else 0L
    label <- if(unresolved) "UNRESOLVED" else if(sum(mechanisms)>1L) "MIXED_MECHANISM" else if(confirmed) "CONFIRMED_SAMPLED_WEEKLY_SWITCH" else if(displacement) "CONFIDENT_TRACKED_DISPLACEMENT" else if(topology) "SAMPLED_PEAK_CREATION_LOSS" else if(window) "ANNUAL_WINDOW_SELECTION" else "UNRESOLVED"
    movement <- if(nrow(ww)>1L) abs(ww$position[which.max(ww$a)]-ww$position[which.min(ww$a)]) else NA_real_
    first <- bb[bb$winner & bb$a==min(ww$a),,drop=FALSE]
    last <- bb[bb$winner & bb$a==max(ww$a),,drop=FALSE]
    tracked <- if(nrow(first)==1L && nrow(last)==1L && identical(first$branch,last$branch))
      abs(last$position-first$position) else NA_real_
    data.frame(year=year,classification=label,confirmed_switch=confirmed,confident_displacement=displacement,
      sampled_peak_creation_loss=topology,annual_window_selection=window,
      consequential_unresolved=unresolved,total_annual_argmax_movement=movement,
      total_tracked_branch_movement=tracked,
      largest_adjacent_grid_movement=if(nrow(z)) max(z$weekly_argmax_movement) else NA_real_,
      identity_ambiguities=sum(z$matching_uncertain),
      ambiguity_intervals=paste(sprintf("[%.8f,%.8f]",z$a_left[z$matching_uncertain],z$a_right[z$matching_uncertain]),collapse=";"),
      signed_height_reversals=reversals,
      confirmed_global_exchanges=exchanges,
      crossing_evidence=if(exchanges>0L) "CONFIRMED_SAMPLED_CROSSING" else if(reversals>0L) "REVERSAL_WITHOUT_CONFIRMED_EXCHANGE" else "NO_SAMPLED_REVERSAL_DEMONSTRATED",
      supporting_intervals=paste(sprintf("[%.8f,%.8f]:%s",z$a_left[z$classification!="NO_DETECTED_MOVEMENT"],z$a_right[z$classification!="NO_DETECTED_MOVEMENT"],z$classification[z$classification!="NO_DETECTED_MOVEMENT"]),collapse=";"),stringsAsFactors=FALSE)
  })
  list(unresolved_audit=audit,mechanism_evidence=do.call(rbind,evidence),settings=s)
}

p10c_report <- function(result,out,s=p10c_defaults()) {
  x <- p10c_audit(result,s)
  utils::write.csv(x$unresolved_audit,file.path(out,"unresolved_audit.csv"),row.names=FALSE)
  utils::write.csv(x$mechanism_evidence,file.path(out,"mechanism_evidence.csv"),row.names=FALSE)
  b <- result$branches
  keep <- intersect(c("year","scenario","a","branch","peak_id","winner","position","quadratic_position",
    "quadratic_status","support","genuine","dt_da","dt_dloga","dt_da_units","dt_dloga_units",
    "sensitivity_method","sensitivity_change_under_refinement","nested_sensitivity_change",
    "response_numerator","response_curvature_denominator","surrogate_dt_da"),names(b))
  sensitivity <- b[,keep,drop=FALSE]
  if(nrow(sensitivity)) {
    sensitivity$interpolation_valid <- sensitivity$quadratic_status=="OK"
    sensitivity$remaining_numerical_uncertainty <- ifelse(is.na(sensitivity$dt_da),
      "IDENTITY_SUPPORT_OR_INTERPOLATION_INADEQUATE",ifelse(is.na(sensitivity$nested_sensitivity_change),
      "NO_NESTED_GRID_COMPARISON","FINITE_GRID_DIFFERENCE_ONLY"))
  }
  utils::write.csv(sensitivity,file.path(out,"sensitivity_robustness.csv"),row.names=FALSE)
  if(!is.null(result$phase10c)) {
    utils::write.csv(result$phase10c$decisions,file.path(out,"phase10c_refinement_decisions.csv"),row.names=FALSE)
    utils::write.csv(result$phase10c$history,file.path(out,"phase10c_refinement_history.csv"),row.names=FALSE)
  }
  lines <- c("# Phase 10C saved-result audit","",paste("Source Phase 10 hash:",result$source_hash),
    paste("Scenarios in source:",length(result$records)),paste("Unresolved intervals:",nrow(x$unresolved_audit)),
    paste("Consequential unresolved:",sum(x$unresolved_audit$scientifically_consequential)),"",
    "Priority 1: winner and competitor; 2: winner; 3: competitor in high-sensitivity region; 4: raw topology change in region; 9: other.",
    "Eligible refinement intervals are consequential and intersect the predeclared high-sensitivity region.",
    "Classifications describe sampled weekly evidence only. An unresolved consequential interval prevents a year-wide confirmed label.","",
    "## Year evidence")
  for(i in seq_len(nrow(x$mechanism_evidence))) {y <- x$mechanism_evidence[i,];lines <- c(lines,
    paste0("- ",y$year,": ",y$classification,"; annual movement ",y$total_annual_argmax_movement,
      " weeks; tracked branch movement ",y$total_tracked_branch_movement," weeks; ",y$crossing_evidence,"."))}
  writeLines(lines,file.path(out,"audit_summary.md"));invisible(x)
}

# A targeted run starts from a saved Phase 10 grid and keeps its diagnostic settings.
p10c_refine <- function(source,data,prior,s=p10c_defaults(),root=".") {
  s <- p10c_validate(s)
  if(!identical(prior$type,"phase10_peak_identity")) stop("Saved Phase 10 result required")
  if(!identical(unname(source$artifact_hash),unname(prior$source_hash))) stop("Phase 10 source hash changed")
  settings <- prior$diagnostic_settings
  baseline <- p8_source(source,data,settings,root)
  check <- p10_settings(settings,baseline)
  profile <- p8_profile(baseline,settings);profile$kind <- "fitted"
  records <- prior$records
  original <- length(records);history <- list();decisions <- list()
  counts <- setNames(c(0L,0L),c("2018","2019"))
  for(level in seq_len(s$levels)) {
    collected <- p10_collect(records,settings)
    terminal <- p10_edges(records,collected$peaks,collected$windows,settings)
    edges <- terminal$edges
    snapshot <- prior;snapshot$records <- records;snapshot$terminal <- terminal
    snapshot$unresolved <- edges[edges$refine | edges$classification=="UNRESOLVED",,drop=FALSE]
    if(nrow(snapshot$unresolved)) snapshot$unresolved$reason <- "TARGETED_REFINEMENT_PENDING"
    audit <- p10c_audit(snapshot,s)$unresolved_audit
    candidates <- audit[audit$eligible,,drop=FALSE]
    if(nrow(candidates)) candidates <- candidates[order(candidates$priority,candidates$year,candidates$a_left),,drop=FALSE]
    history[[length(history)+1L]] <- edges
    added <- 0L
    for(i in seq_len(nrow(candidates))) {
      x <- candidates[i,];mid <- (x$a_left+x$a_right)/2
      existing <- any(abs(vapply(records,function(r) r$a,numeric(1))-mid)<1e-14)
      reason <- if(existing) "ALREADY_PRESENT" else if(counts[as.character(x$year)]>=s$max_per_year) "YEAR_CAP" else
        if(length(records)-original>=s$max_new_scenarios || length(records)>=settings$interventions$max_scenarios) "SCENARIO_CAP" else "SELECTED"
      decisions[[length(decisions)+1L]] <- data.frame(level=level,year=x$year,a_left=x$a_left,a_right=x$a_right,
        midpoint=mid,priority=x$priority,decision=reason,stringsAsFactors=FALSE)
      if(reason!="SELECTED") next
      id <- sprintf("scenario-%04d",max(as.integer(sub("scenario-","",names(records))))+1L)
      r <- p8_evaluate(baseline,check,profile,baseline$parameters$S0,mid,id,"a","phase10c_refined",level,
        paste(x$a_left,x$a_right,sep="/"))
      if(r$status!="VALID") stop("Phase 10C scenario failed: ",r$reason)
      records[[id]] <- r;counts[as.character(x$year)] <- counts[as.character(x$year)]+1L;added <- added+1L
    }
    if(!added) break
  }
  collected <- p10_collect(records,settings);terminal <- p10_edges(records,collected$peaks,collected$windows,settings)
  result <- prior;result$records <- records;result$peaks <- collected$peaks;result$windows <- collected$windows
  result$terminal <- terminal;result$branches <- p10_branches(records,collected$peaks,terminal$matches,collected$windows,settings)
  result$unresolved <- terminal$edges[terminal$edges$refine | terminal$edges$classification=="UNRESOLVED",,drop=FALSE]
  if(nrow(result$unresolved)) result$unresolved$reason <- if(length(records)-original>=s$max_new_scenarios) "PHASE10C_SCENARIO_CAP" else "PHASE10C_REMAINING_UNRESOLVED"
  result$type <- "phase10_peak_identity";result$status <- "PHASE10C_TARGETED_GRID_STORED"
  result$phase10c <- list(settings=s,source_phase10_hash=prior$artifact_hash,new_scenarios=length(records)-original,
    decisions=if(length(decisions)) do.call(rbind,decisions) else data.frame(),
    history=if(length(history)) do.call(rbind,history) else data.frame(),year_additions=counts)
  result
}
