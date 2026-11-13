# Fine-grid diagnostic on the accepted Phase 7C endpoint; no fitted model changes.
p8b_defaults <- function() {
  s <- p8_defaults()
  grid <- (40:60)/50
  s$time$years <- 2014:2019
  s$time$include_incomplete <- FALSE
  s$interventions$s_values <- NULL
  s$interventions$s_multipliers <- grid
  s$interventions$a_values <- grid
  s$interventions$joint <- FALSE
  s$interventions$max_scenarios <- 160L
  s$annual$enabled <- FALSE
  s$refinement <- list(enabled=TRUE,levels=4L,max_points=118L)
  s
}
p8b_settings <- function(s,baseline) {
  if(!identical(s$time$years,2014:2019) || isTRUE(s$time$include_incomplete) ||
     !is.null(s$time$horizon_rows) || isTRUE(s$annual$enabled) || isTRUE(s$interventions$joint))
    stop("Phase 8B requires full supported horizon, complete 2014-2019, fitted forcing, and independent axes")
  grid <- (40:60)/50
  if(!identical(s$interventions$a_values,grid) || !identical(s$interventions$s_multipliers,grid) ||
      !is.null(s$interventions$s_values)) stop("Phase 8B prespecified 0.80:0.02:1.20 grids required")
  check <- p8_settings(s,baseline)
  if(!all(check$complete)) stop("Phase 8B requires all six complete annual windows")
  if(!identical(check$s_values,baseline$parameters$S0*grid)) stop("Saved S0 ratio grid changed")
  check
}
p8b_pair_differences <- function(edge) {
  t <- edge$matching$table
  if(!nrow(t) || !all(c("status","left_peak","right_peak") %in% names(t))) return(data.frame())
  t <- t[t$status=="MATCHED",,drop=FALSE]
  if(nrow(t)<2L) return(data.frame())
  L <- edge$matching$L;R <- edge$matching$R
  rows <- list()
  for(i in seq_len(nrow(t)-1L)) for(j in (i+1L):nrow(t)) {
    h <- function(tab,id) tab$height[match(id,tab$peak_id)]
    dl <- h(L,t$left_peak[i])-h(L,t$left_peak[j])
    dr <- h(R,t$right_peak[i])-h(R,t$right_peak[j])
    rows[[length(rows)+1L]] <- data.frame(left=edge$left,right=edge$right,year=edge$year,axis=edge$axis,
      theta_left=edge$row$theta_left,theta_right=edge$row$theta_right,
      branch1_left_peak=t$left_peak[i],branch1_right_peak=t$right_peak[i],
      branch2_left_peak=t$left_peak[j],branch2_right_peak=t$right_peak[j],
      branch1_height_left=h(L,t$left_peak[i]),branch1_height_right=h(R,t$right_peak[i]),
      branch2_height_left=h(L,t$left_peak[j]),branch2_height_right=h(R,t$right_peak[j]),
      height_difference_left=dl,height_difference_right=dr,
      sign_change=(dl*dr<0),matching_ambiguous=isTRUE(edge$matching$ambiguous))
  }
  do.call(rbind,rows)
}
p8b_edge <- function(left,right,axis,year,policy,stage="fine",level=0L) {
  e <- p8_edge(left,right,axis,year,policy,stage,level)
  pairs <- p8b_pair_differences(e)
  t <- e$matching$table
  displacement <- if(nrow(t) && "circular_distance" %in% names(t)) {
    v <- t$circular_distance[t$status=="MATCHED"]
    if(length(v)) max(v) else NA_real_
  } else NA_real_
  e$row$max_matched_branch_displacement <- displacement
  e$row$height_crossing <- nrow(pairs)>0 && any(pairs$sign_change & !pairs$matching_ambiguous)
  e$row$raw_count_change <- isTRUE(e$row$raw_peaks_left!=e$row$raw_peaks_right)
  e$row$refinement_trigger <- isTRUE(e$row$global_circular_displacement>=8) ||
    e$row$height_crossing || e$row$raw_count_change || isTRUE(e$matching$ambiguous)
  e$row$pair_height_differences <- nrow(pairs)
  e$pair_differences <- pairs
  e
}
p8b_group <- function(records,axis) {
  ids <- names(records)[vapply(records,function(r) r$axis==axis,logical(1))]
  values <- vapply(records[ids],function(r) if(axis=="s") r$s else r$a,numeric(1))
  ids[order(values)]
}
p8b_edges <- function(records,years,policy,level=0L) {
  edges <- list()
  for(axis in c("a","s")) {
    ids <- p8b_group(records,axis)
    if(length(ids)<2L) next
    for(i in seq_len(length(ids)-1L)) for(year in years)
      edges[[length(edges)+1L]] <- p8b_edge(records[[ids[i]]],records[[ids[i+1L]]],axis,year,policy,
        if(level==0L) "fine" else "refined",level)
  }
  edges
}
p8b_refine <- function(records,years,policy,settings,evaluator) {
  history <- list();unresolved <- list();current <- p8b_edges(records,years,policy)
  history <- c(history,current)
  if(!settings$refinement$enabled || settings$refinement$levels==0L)
    return(list(records=records,history=history,terminal=current,unresolved=data.frame(),status="NOT_REQUESTED"))
  added <- 0L
  for(level in seq_len(settings$refinement$levels)) {
    candidate <- Filter(function(e) isTRUE(e$row$refinement_trigger),current)
    if(!length(candidate)) break
    keys <- vapply(candidate,function(e) paste(e$axis,e$left,e$right,sep="/"),character(1))
    candidate <- candidate[!duplicated(keys)]
    new_ids <- character()
    for(e in candidate) {
      if(added>=settings$refinement$max_points || length(records)>=settings$interventions$max_scenarios) {
        unresolved[[length(unresolved)+1L]] <- data.frame(axis=e$axis,theta_left=e$row$theta_left,
          theta_right=e$row$theta_right,left=e$left,right=e$right,level=level,reason="SCENARIO_CAP")
        next
      }
      theta <- (e$row$theta_left+e$row$theta_right)/2
      if(theta %in% c(e$row$theta_left,e$row$theta_right)) next
      id <- sprintf("scenario-%04d",length(records)+1L)
      left <- records[[e$left]]
      records[[id]] <- evaluator(if(e$axis=="s") theta else left$s,
        if(e$axis=="a") theta else left$a,id,e$axis,level,paste(e$left,e$right,sep="/"))
      added <- added+1L;new_ids <- c(new_ids,id)
    }
    if(!length(new_ids)) break
    all_edges <- p8b_edges(records,years,policy,level)
    current <- Filter(function(e) e$left %in% new_ids || e$right %in% new_ids,all_edges)
    history <- c(history,current)
  }
  terminal <- p8b_edges(records,years,policy,settings$refinement$levels)
  list(records=records,history=history,terminal=terminal,unresolved=if(length(unresolved)) do.call(rbind,unresolved) else data.frame(),
    status=if(length(unresolved)) "SCENARIO_CAP_UNRESOLVED" else "REFINEMENT_COMPLETED")
}
p8b_critical <- function(fine,terminal,unresolved) {
  rows <- list()
  for(e in fine) {
    if(!isTRUE(e$row$refinement_trigger)) next
    t <- Filter(function(x) x$axis==e$axis && x$year==e$year &&
      x$row$theta_left>=e$row$theta_left && x$row$theta_right<=e$row$theta_right,terminal)
    persist <- Filter(function(x) isTRUE(x$row$refinement_trigger),t)
    chosen <- if(length(persist)) persist[[which.min(vapply(persist,function(x) x$row$parameter_interval,numeric(1)))]] else e
    cap <- nrow(unresolved)>0 && any(unresolved$axis==e$axis & unresolved$theta_left>=e$row$theta_left &
      unresolved$theta_right<=e$row$theta_right)
    rows[[length(rows)+1L]] <- data.frame(year=e$year,axis=e$axis,initial_left=e$row$theta_left,
      initial_right=e$row$theta_right,narrowest_left=chosen$row$theta_left,narrowest_right=chosen$row$theta_right,
      initial_mechanism=e$row$mechanism,terminal_mechanism=chosen$row$mechanism,
      persists_under_refinement=length(persist)>0,unresolved_by_cap=cap)
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
p8b_run <- function(source,data,s=p8b_defaults(),root=".") {
  baseline <- p8_source(source,data,s,root);check <- p8b_settings(s,baseline)
  profile <- p8_profile(baseline,s);profile$kind <- "fitted"
  if(42L>s$interventions$max_scenarios) stop("Fine grid exceeds scenario cap")
  records <- list();add <- function(sval,a,id,axis,level=0L,parent="")
    p8_evaluate(baseline,check,profile,sval,a,id,axis,if(level==0L) "fine" else "refined",level,parent)
  for(axis in c("a","s")) for(value in if(axis=="a") check$a_values else check$s_values) {
    id <- sprintf("scenario-%04d",length(records)+1L)
    records[[id]] <- add(if(axis=="s") value else baseline$parameters$S0,if(axis=="a") value else 1,id,axis)
  }
  if(any(vapply(records,function(r) r$status!="VALID",logical(1)))) stop("Fine-grid scenario failed; inspect valid simplex and simulator")
  ref <- baseline$evaluation$states[seq_len(check$horizon),,drop=FALSE]
  for(axis in c("a","s")) {
    r <- records[[p8b_group(records,axis)[11L]]]
    if(!p6r_close(as.numeric(as.matrix(r$trajectory[,names(ref)])),as.numeric(as.matrix(ref)),s$replay_tolerance))
      stop("Exact baseline replay failed")
  }
  fine <- p8b_edges(records,check$years,s$peaks)
  refined <- p8b_refine(records,check$years,s$peaks,s,add)
  track <- p8_track(refined$records,p8_series(refined$records),check$years,s$peaks,"all_computed")
  list(type="phase8b_fine_grid",status="FINITE_WEEKLY_DIAGNOSTICS_STORED",config=source$config,
    diagnostic_settings=s,source_hash=source$artifact_hash,source_identity=source$identity,
    source_config=source$config,source_ancestry=baseline$source_ancestry,
    baseline_replay=baseline$replay,baseline_parameters=baseline$parameters,
    time_contract=check,records=refined$records,fine_edges=fine,edge_history=refined$history,
    terminal_edges=refined$terminal,branches=track$branches,
    critical=p8b_critical(fine,refined$terminal,refined$unresolved),
    cap_unresolved=refined$unresolved,refinement_status=refined$status)
}
p8b_report <- function(result,out) {
  if(!identical(result$type,"phase8b_fine_grid")) stop("Saved Phase 8B result required")
  write <- function(x,name) utils::write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE)
  rows <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
  records <- result$records
  scenarios <- rows(lapply(records,function(r) data.frame(scenario=r$id,axis=r$axis,
    parameter=if(r$axis=="a") r$a else r$s,ratio=if(r$axis=="a") r$a else r$s/result$baseline_parameters$S0,
    s=r$s,a=r$a,stage=r$stage,level=r$level,parent=r$parent,status=r$status,reason=r$reason)))
  annual <- rows(lapply(records,function(r) {
    t <- r$annual_summary;t$scenario <- r$id;t$axis <- r$axis
    t$parameter <- if(r$axis=="a") r$a else r$s;t$stage <- r$stage;t$level <- r$level;t
  }))
  local <- list();for(r in records) for(year in names(r$annual)) {
    t <- r$annual[[year]]$local
    if(!nrow(t)) next
    t$scenario <- r$id;t$year <- as.integer(year);t$axis <- r$axis
    t$parameter <- if(r$axis=="a") r$a else r$s;t$stage <- r$stage
    local[[length(local)+1L]] <- t
  }
  edges <- result$edge_history
  transitions <- rows(lapply(edges,function(e) e$row))
  terminal_keys <- vapply(result$terminal_edges,function(e) paste(e$axis,e$year,e$left,e$right),character(1))
  if(nrow(transitions)) transitions$terminal <- vapply(edges,function(e)
    paste(e$axis,e$year,e$left,e$right) %in% terminal_keys,logical(1))
  pairs <- rows(lapply(edges,function(e) e$pair_differences))
  matching <- rows(lapply(edges,function(e) {
    t <- e$matching$table
    if(!nrow(t)) return(NULL)
    t$left <- e$left;t$right <- e$right;t$axis <- e$axis;t$year <- e$year
    t$theta_left <- e$row$theta_left;t$theta_right <- e$row$theta_right
    t$matching_ambiguous <- isTRUE(e$matching$ambiguous)
    t
  }))
  branch <- result$branches
  if(nrow(pairs) && nrow(branch)) {
    lookup <- function(scenario,year,peak) {
      x <- branch$branch[branch$scenario==scenario & branch$year==year & branch$peak_id==peak]
      if(length(x)) x[1] else NA_character_
    }
    pairs$branch1 <- mapply(lookup,pairs$left,pairs$year,pairs$branch1_left_peak)
    pairs$branch2 <- mapply(lookup,pairs$left,pairs$year,pairs$branch2_left_peak)
  }
  counts <- list()
  for(year in result$time_contract$years) for(axis in c("a","s")) {
    t <- rows(lapply(Filter(function(e) e$year==year && e$axis==axis,result$terminal_edges),function(e) e$row))
    if(!nrow(t)) next
    capped <- vapply(seq_len(nrow(t)),function(i) nrow(result$cap_unresolved)>0 &&
      any(result$cap_unresolved$axis==axis & result$cap_unresolved$theta_left<=t$theta_left[i] &
        result$cap_unresolved$theta_right>=t$theta_right[i]),logical(1))
    counts[[length(counts)+1L]] <- data.frame(year=year,axis=axis,
      max_adjacent_global_week_difference=max(t$global_circular_displacement,na.rm=TRUE),
      max_adjacent_matched_branch_displacement=if(all(is.na(t$max_matched_branch_displacement))) NA_real_ else max(t$max_matched_branch_displacement,na.rm=TRUE),
      switching_candidates=sum(t$mechanism=="PEAK_SWITCHING"),
      displacement_candidates=sum(t$mechanism=="PEAK_DISPLACEMENT"),
      creation_annihilation_candidates=sum(t$mechanism=="PEAK_CREATION_ANNIHILATION"),
      unresolved_transitions=sum(t$mechanism=="AMBIGUOUS_UNRESOLVED" | capped),
      intervals_unresolved_by_cap=sum(capped))
  }
  write(scenarios,"fine_grid_scenarios");write(annual,"fine_grid_annual_peaks")
  write(rows(local),"fine_grid_local_peaks");write(branch,"fine_grid_branches")
  write(transitions,"fine_grid_transitions");write(result$critical,"critical_intervals")
  write(rows(counts),"mechanism_counts")
  write(pairs,"competing_peak_heights");write(matching,"fine_grid_peak_matching")
  write(result$cap_unresolved,"cap_unresolved_intervals")
  grDevices::pdf(file.path(out,"phase8b_diagnostics.pdf"),width=11,height=8)
  tryCatch({
    for(axis in c("a","s")) for(year in result$time_contract$years) {
      a <- annual[annual$axis==axis & annual$year==year & annual$status=="ANALYZED",,drop=FALSE]
      a <- a[order(a$parameter),,drop=FALSE]
      if(!nrow(a)) next
      graphics::plot(a$parameter,a$global_week,xlab=if(axis=="a") "Transmission multiplier" else "Initial susceptible fraction",
        ylab="Global incidence peak week",ylim=c(1,52),pch=ifelse(a$stage=="fine",16,1),
        main=paste(year,axis,"global weekly peaks; open points are refined"))
      b <- branch[branch$year==year & grepl(paste0(":",axis,":"),branch$group,fixed=FALSE),,drop=FALSE]
      if(nrow(b)) {
        graphics::plot(b$theta,b$week,xlab=axis,ylab="Local peak week",ylim=c(1,52),
          main=paste(year,axis,"retained local branches"),type="n")
        ids <- unique(b$branch)
        for(i in seq_along(ids)) {v <- b[b$branch==ids[i],,drop=FALSE];v <- v[order(v$theta),]
          graphics::lines(v$theta,v$week,type="b",col=(i-1L)%%8L+1L,pch=16,cex=.5)}
      }
      p <- pairs[pairs$year==year & pairs$axis==axis,,drop=FALSE]
      if(nrow(p)) {
        yrange <- range(c(p$branch1_height_left,p$branch1_height_right,p$branch2_height_left,p$branch2_height_right))
        graphics::plot(range(c(p$theta_left,p$theta_right)),yrange,type="n",xlab=axis,ylab="Peak height",
          main=paste(year,axis,"competing matched heights; segments connect computed endpoints"))
        for(i in seq_len(nrow(p))) {
          graphics::segments(p$theta_left[i],p$branch1_height_left[i],p$theta_right[i],p$branch1_height_right[i],col=2)
          graphics::segments(p$theta_left[i],p$branch2_height_left[i],p$theta_right[i],p$branch2_height_right[i],col=4)
        }
      }
      t <- transitions[transitions$axis==axis & transitions$year==year & transitions$terminal,,drop=FALSE]
      if(nrow(t)) {
        mid <- (t$theta_left+t$theta_right)/2
        labels <- unique(t$mechanism)
        graphics::plot(mid,seq_along(mid),col=match(t$mechanism,labels),pch=19,xlab=axis,
          ylab="Adjacent interval",main=paste(year,axis,"terminal classifications"))
        graphics::segments(t$theta_left,seq_along(mid),t$theta_right,seq_along(mid),col=match(t$mechanism,labels))
        graphics::legend("topright",labels,col=seq_along(labels),pch=19,bty="n",cex=.7)
      }
      candidates <- Filter(function(e) e$year==year && e$axis==axis && isTRUE(e$row$refinement_trigger),result$terminal_edges)
      for(e in candidates) {
        left <- records[[e$left]]$trajectory;right <- records[[e$right]]$trajectory
        left <- left[left$year==year & !left$initialization_only,];right <- right[right$year==year & !right$initialization_only,]
        graphics::plot(left$week,left$incidence,type="l",col=2,xlab="Week",ylab="Incident infections",
          ylim=range(c(left$incidence,right$incidence)),main=paste(year,axis,sprintf("[%.12g, %.12g]",e$row$theta_left,e$row$theta_right),e$row$mechanism))
        graphics::lines(right$week,right$incidence,col=4)
        graphics::legend("topright",c("left computed curve","right computed curve"),col=c(2,4),lty=1,bty="n")
      }
    }
  },finally=grDevices::dev.off())
  writeLines(c("# Phase 8B fine-grid diagnostics","",
    paste("Accepted Phase 7C source MD5:",result$source_hash),
    paste("Replay:",result$baseline_replay),
    paste("Refinement:",result$refinement_status),
    paste("Scenarios:",length(records)),
    "Complete 2014–2019 annual windows are simulated within the full saved horizon; state continues across years.",
    "The a and s/s0 grids use 0.80:0.02:1.20, with exact baseline 1. Each scenario starts from the saved initial state with only its declared intervention changed.",
    "Classification is finite weekly evidence. Ambiguous identities, near ties, boundary peaks, and unresolved capped intervals retain explicit uncertainty.",
    "No mathematical discontinuity, continuous-time bifurcation, causal explanation, or fitted-model change is asserted.",
    "Height difference sign changes and all intermediate numerical values are in competing_peak_heights.csv; transition intervals and refinement status are in critical_intervals.csv.",
    "The PDF connects computed peak summaries to show branch identity; it never interpolates epidemic trajectories.",
    "Saved-only reporting performs no simulation, fitting, peak detection, or refinement.","",
    "## Exact diagnostic settings","```",capture.output(dput(result$diagnostic_settings)),"```"),
    file.path(out,"phase8b_summary.md"))
  invisible(result)
}
