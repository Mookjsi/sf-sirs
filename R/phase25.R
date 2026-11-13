# Targeted audit of the saved Phase 24 weekly intervention trajectories.
p25_verify <- function(cfg,source_cfg,p24cfg) {
  if(!identical(cfg$padding_weeks,p24cfg$padding_weeks)) stop("Phase 24 support window changed")
  report <- read_phase2_result(cfg$phase24_report)
  if(!identical(report$metadata$phase,24L) || !identical(report$metadata$stage,"report") ||
     !identical(report$result$status,"REPORT_BUILT")) stop("Wrong Phase 24 report")
  paths <- report$metadata$input_paths
  if(anyNA(tools::md5sum(paths)) ||
     !identical(unname(tools::md5sum(paths)),unname(report$metadata$input_hashes)))
    stop("Phase 24 report dependency hash mismatch")
  if(!identical(unname(tools::md5sum(cfg$phase24_adjudicate)),report$result$adjudicate_hash))
    stop("Phase 24 adjudication identity mismatch")
  x <- p24_read_stage(cfg$phase24_adjudicate,"adjudicate",p24cfg)
  if(!identical(x$status,"ADJUDICATED") || !all(x$replay$passed)) stop("Phase 24 replay gate failed")
  files <- c("interval_classifications.csv","all_peak_branches.csv","intervention_grid.csv","state_swaps.csv")
  for(f in files) if(!file.exists(file.path(cfg$phase24_dir,f))) stop("Missing Phase 24 report file: ",f)
  published_grid <- utils::read.csv(file.path(cfg$phase24_dir,"intervention_grid.csv"))
  published_swaps <- utils::read.csv(file.path(cfg$phase24_dir,"state_swaps.csv"))
  if(!identical(published_grid$scenario,x$grid$scenario) ||
     max(abs(published_grid$v-x$grid$v))>1e-12 ||
     !identical(published_swaps$swap,x$swaps$swap) ||
     max(abs(published_swaps$S-x$swaps$S),na.rm=TRUE)>cfg$replay_tolerance)
    stop("Published Phase 24 grid or state swaps differ from saved adjudication")
  edges <- utils::read.csv(file.path(cfg$phase24_dir,files[1]))
  targets <- edges[edges$reference=="JOINT" & edges$season==2017L &
    edges$classification=="UNRESOLVED" &
    ((edges$path=="BETA_ONLY" & edges$x_left < -.10 & edges$x_right > -.10) |
     (edges$path=="JOINT" & edges$x_left > 0 & edges$x_right > .08)),]
  if(nrow(targets)!=2L || !setequal(targets$path,c("BETA_ONLY","JOINT")))
    stop("Target interval identity changed")
  targets <- targets[match(c("BETA_ONLY","JOINT"),targets$path),]
  targets$target <- c("A","B")
  for(i in seq_len(nrow(targets))) {
    t <- targets[i,]
    a <- x$grid[match(c(t$left,t$right),x$grid$scenario),]
    if(anyNA(a$scenario) || any(a$reference!="JOINT" | a$season!=2017L | a$path!=t$path) ||
       max(abs(a$v-c(t$x_left,t$x_right)))>1e-12 || abs(a$u[1]-a$u[2])>1e-12)
      stop("Target endpoint coordinates changed")
    targets$u[i] <- a$u[1]
  }
  list(saved=x,targets=targets,hashes=tools::md5sum(c(cfg$phase24_report,
    cfg$phase24_adjudicate,file.path(cfg$phase24_dir,files),paths)))
}

p25_independent_peaks <- function(y,tie_abs=1e-8,tie_rel=1e-10) {
  tol <- tie_abs+tie_rel*max(abs(y))
  which(vapply(seq_along(y),function(i)
    i>1L && i<length(y) && y[i]>y[i-1L]+tol && y[i]>y[i+1L]+tol,logical(1)))
}

p25_run <- function(score,inputs,clock,reference,target,path,u,v,id,cfg) {
  z <- p24_scenario(score,inputs,clock,2017L,NULL,u,v,padding=cfg$padding_weeks)
  if(z$status!="VALID") stop("Phase 25 scenario invalid: ",id)
  p <- p24_peaks(z,id,reference,2017L,path,u,v,0L)
  p$target <- target
  p$boundary_status <- ifelse(!p$support,"SUPPORT_BOUNDARY",
    ifelse(!p$in_window,"OUTSIDE_PRIMARY_WINDOW","INTERIOR"))
  a <- z$states$expected_cases
  curve <- data.frame(target=target,reference=reference,scenario=id,path=path,u=u,v=v,
    season_week=seq_len(length(a))-1L,expected=a,
    in_window=seq_along(a)>=2L & seq_along(a)<=53L,
    support=seq_along(a)>1L & seq_along(a)<length(a))
  list(grid=data.frame(target=target,reference=reference,scenario=id,path=path,u=u,v=v,
    selected_week=z$selected_position,selected_height=z$selected_height,
    peak_count=sum(p$genuine & p$in_window),selected_branch=NA_character_,
    match_status="UNTRACKED"),peaks=p,curve=curve,raw=z)
}

p25_replay <- function(checked,scores,inputs,clock,cfg) {
  rows <- list();counter <- 0L
  for(ti in seq_len(nrow(checked$targets))) for(side in c("left","right")) {
    t <- checked$targets[ti,];id <- t[[side]]
    old <- checked$saved$grid[checked$saved$grid$scenario==id,]
    if(nrow(old)!=1L) stop("Missing saved endpoint")
    counter <- counter+1L
    z <- p25_run(scores$JOINT,inputs,clock,"JOINT",t$target,t$path,t$u,
                 old$v,id,cfg)
    op <- checked$saved$peaks[checked$saved$peaks$scenario==id,]
    np <- z$peaks
    op <- op[order(op$position),];np <- np[order(np$position),]
    same <- nrow(np)==nrow(op) && identical(np$position,op$position) &&
      identical(np$genuine,op$genuine) && identical(np$in_window,op$in_window)
    err <- if(same) max(abs(np$height-op$height)) else Inf
    old_curve <- checked$saved$curves[checked$saved$curves$scenario==id,]
    curve_err <- if(nrow(old_curve)==52L) max(abs(z$curve$expected[2:53]-old_curve$expected)) else Inf
    serr <- max(abs(c(z$grid$selected_week-old$selected_position,
                      z$grid$selected_height-old$selected_height,curve_err)))
    published <- utils::read.csv(file.path(cfg$phase24_dir,"all_peak_branches.csv"))
    published <- published[published$scenario==id,]
    published <- published[order(published$position),]
    branch_ok <- same && identical(op$branch,published$branch)
    rows[[counter]] <- data.frame(target=t$target,scenario=id,side=side,
      peak_count=nrow(np),max_peak_height_error=err,selected_error=serr,
      branch_ids=paste(op$branch,collapse=";"),branch_verified=branch_ok,
      passed=same && branch_ok && err<=cfg$replay_tolerance && serr<=cfg$replay_tolerance)
  }
  ans <- do.call(rbind,rows)
  if(!all(ans$passed)) stop("Phase 24 endpoint replay failed")
  ans
}

p25_track <- function(runs,cfg) {
  g <- do.call(rbind,lapply(runs,`[[`,"grid"));p <- do.call(rbind,lapply(runs,`[[`,"peaks"))
  g <- g[order(g$v),];p$match_status <- "NEW";p$branch <- NA_character_
  next_id <- 0L;previous <- NULL
  for(k in seq_len(nrow(g))) {
    ix <- which(p$scenario==g$scenario[k]);cur <- p[ix,,drop=FALSE]
    if(!is.null(previous) && nrow(cur) && nrow(previous)) {
      a <- p10_assignment(previous,cur,cfg$identity)
      b <- p10_assignment(previous,cur,cfg$identity,TRUE)
      for(j in seq_len(nrow(a))) if(a$status[j]=="MATCHED" && !is.na(a$right_peak[j])) {
        reverse <- b[b$left_peak==a$left_peak[j] & b$right_peak==a$right_peak[j],]
        right <- match(a$right_peak[j],cur$peak_id)
        left <- match(a$left_peak[j],previous$peak_id)
        if(nrow(reverse)==1L && reverse$status=="MATCHED" && !is.na(left)) {
          cur$branch[right] <- previous$branch[left];cur$match_status[right] <- "MATCHED"
        } else cur$match_status[right] <- "AMBIGUOUS"
      }
      if(any(a$status=="AMBIGUOUS")) cur$match_status[match(na.omit(a$right_peak[a$status=="AMBIGUOUS"]),cur$peak_id)] <- "AMBIGUOUS"
    }
    for(j in which(is.na(cur$branch))) {
      next_id <- next_id+1L;cur$branch[j] <- paste0("tracked-",next_id)
    }
    p[ix,] <- cur
    eligible <- cur[cur$genuine & cur$in_window,]
    selected <- eligible[eligible$position==g$selected_week[k],]
    if(nrow(selected)==1L) {
      g$selected_branch[k] <- selected$branch
      g$match_status[k] <- selected$match_status
    } else g$match_status[k] <- "WINDOW_OR_AMBIGUOUS"
    previous <- cur
  }
  list(grid=g,peaks=p)
}

p25_gap <- function(tracked,cfg) {
  g <- tracked$grid;p <- tracked$peaks;rows <- list()
  for(k in seq_len(nrow(g))) {
    q <- p[p$scenario==g$scenario[k] & p$genuine & p$in_window,]
    early <- q[q$position<=cfg$early_max_week,]
    late <- q[q$position>=cfg$late_min_week,]
    if(nrow(early)!=1L || nrow(late)!=1L) next
    rows[[length(rows)+1L]] <- data.frame(target=g$target[k],reference=g$reference[k],
      scenario=g$scenario[k],v=g$v[k],early_branch=early$branch,
      late_branch=late$branch,early_week=early$position,late_week=late$position,
      early_height=early$height,late_height=late$height,
      signed_gap=early$height-late$height,
      trackable=all(c(early$match_status,late$match_status)!="AMBIGUOUS"))
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}

p25_transitions <- function(tracked,gaps,cfg) {
  g <- tracked$grid;p <- tracked$peaks;out <- list()
  for(k in seq_len(nrow(g)-1L)) {
    a <- g[k,];b <- g[k+1L,]
    ga <- gaps[gaps$scenario==a$scenario,];gb <- gaps[gaps$scenario==b$scenario,]
    exchange <- nrow(ga)==1L && nrow(gb)==1L && ga$trackable && gb$trackable &&
      ga$early_branch==gb$early_branch && ga$late_branch==gb$late_branch &&
      ga$signed_gap*gb$signed_gap<0
    changed <- a$peak_count!=b$peak_count || !identical(a$selected_branch,b$selected_branch) ||
      a$selected_week!=b$selected_week || exchange ||
      any(p$match_status[p$scenario %in% c(a$scenario,b$scenario)]=="AMBIGUOUS")
    if(!changed) next
    late_a <- p[p$scenario==a$scenario & p$genuine & p$position>=cfg$late_min_week,]
    late_b <- p[p$scenario==b$scenario & p$genuine & p$position>=cfg$late_min_week,]
    kind <- if(exchange) "DOMINANCE_EXCHANGE" else if(nrow(late_a)>0L && nrow(late_b)==0L) {
      if(any(p$match_status[p$scenario==b$scenario]=="AMBIGUOUS")) "UNRESOLVED" else "BRANCH_LOSS"
    } else if(any(p$match_status[p$scenario %in% c(a$scenario,b$scenario)]=="AMBIGUOUS")) "UNRESOLVED"
    else if(identical(a$selected_branch,b$selected_branch) && a$peak_count==b$peak_count)
      "CONTINUOUS_DISPLACEMENT" else "UNRESOLVED"
    # A candidate just outside the primary window is a window selection only if it is supported.
    if(kind=="BRANCH_LOSS" && any(p$scenario==b$scenario & p$position>52 & p$support)) kind <- "WINDOW_SELECTION"
    loss_detail <- if(kind=="BRANCH_LOSS") {
      nearby <- p[p$scenario==b$scenario & p$genuine,]
      if(nrow(nearby) && min(abs(outer(late_a$position,nearby$position,`-`)))<=1)
        "POSSIBLE_MERGE_IN_WEEKLY_SAMPLES" else "SAMPLED_LOCAL_MAXIMUM_LOSS"
    } else if(kind=="WINDOW_SELECTION") "CROSSED_PRIMARY_WINDOW" else NA_character_
    out[[length(out)+1L]] <- data.frame(target=a$target,reference=a$reference,
      left=a$scenario,right=b$scenario,v_left=a$v,v_right=b$v,width=b$v-a$v,
      left_count=a$peak_count,right_count=b$peak_count,
      left_selected=a$selected_week,right_selected=b$selected_week,
      D_left=if(nrow(ga)) ga$signed_gap else NA_real_,
      D_right=if(nrow(gb)) gb$signed_gap else NA_real_,
      classification=kind,
      loss_detail=loss_detail,
      evidence_scope=if(kind=="BRANCH_LOSS") "SAMPLED_WEEKLY_PEAK_LOSS" else
        if(kind=="WINDOW_SELECTION") "FIXED_PRIMARY_WINDOW" else "WEEKLY_TRAJECTORY")
  }
  if(length(out)) do.call(rbind,out) else data.frame()
}

p25_refine <- function(score,inputs,clock,reference,t,cfg) {
  runs <- list();left <- t$x_left;right <- t$x_right
  for(v in c(left,right)) {
    id <- sprintf("%s-%s-%08d",t$target,reference,length(runs)+1L)
    runs[[length(runs)+1L]] <- p25_run(score,inputs,clock,reference,t$target,t$path,t$u,v,id,cfg)
  }
  stable <- 0L;previous <- NULL
  for(level in seq_len(cfg$maximum_bisections)) {
    track <- p25_track(runs,cfg);gap <- p25_gap(track,cfg);tr <- p25_transitions(track,gap,cfg)
    signature <- if(nrow(tr)) paste(tr$classification,collapse=";") else "NONE"
    stable <- if(identical(signature,previous)) stable+1L else 0L
    previous <- signature
    critical <- if(nrow(tr)) tr else data.frame(v_left=left,v_right=right,width=right-left)
    if(level>cfg$minimum_bisections && all(critical$width<=cfg$transition_width) &&
       stable>=cfg$stable_rounds) break
    candidates <- if(level<=cfg$minimum_bisections) {
      vals <- sort(vapply(runs,function(z) z$grid$v,numeric(1)))
      (head(vals,-1L)+tail(vals,-1L))/2
    } else (critical$v_left+critical$v_right)/2
    candidates <- unique(candidates)
    old <- vapply(runs,function(z) z$grid$v,numeric(1))
    candidates <- candidates[vapply(candidates,function(v) all(abs(old-v)>1e-14),logical(1))]
    if(!length(candidates)) break
    for(v in candidates) {
      id <- sprintf("%s-%s-%08d",t$target,reference,length(runs)+1L)
      runs[[length(runs)+1L]] <- p25_run(score,inputs,clock,reference,t$target,t$path,t$u,v,id,cfg)
    }
  }
  track <- p25_track(runs,cfg);gap <- p25_gap(track,cfg);tr <- p25_transitions(track,gap,cfg)
  list(grid=track$grid,peaks=track$peaks,gaps=gap,transitions=tr,
    curves=do.call(rbind,lapply(runs,`[[`,"curve")),
    converged=nrow(tr)>0L && all(tr$width<=cfg$transition_width) && stable>=cfg$stable_rounds,
    levels=level)
}

p25_sensitivity <- function(result) {
  rows <- list()
  for(i in seq_len(nrow(result$grid))) {
    g <- result$grid[i,]
    curve <- result$curves[result$curves$scenario==g$scenario,]
    ix <- p25_independent_peaks(curve$expected)
    independent <- curve$season_week[ix]
    saved <- result$peaks$position[result$peaks$scenario==g$scenario &
      result$peaks$genuine]
    rows[[i]] <- data.frame(target=g$target,reference=g$reference,scenario=g$scenario,
      v=g$v,independent_peaks=paste(independent,collapse=";"),
      primary_peaks=paste(saved,collapse=";"),
      detector_agrees=identical(as.integer(independent),as.integer(saved)),
      simulator_tolerance="NOT_APPLICABLE_DISCRETE_WEEKLY_RECURRENCE",
      within_week_available=FALSE)
  }
  do.call(rbind,rows)
}

p25_verdict <- function(result) {
  tr <- result$transitions
  kinds <- unique(tr$classification)
  verdict <- if("UNRESOLVED" %in% kinds) "UNRESOLVED" else
    if("DOMINANCE_EXCHANGE" %in% kinds) "DOMINANCE_EXCHANGE" else
    if("BRANCH_LOSS" %in% kinds) "BRANCH_LOSS" else
    if("WINDOW_SELECTION" %in% kinds) "WINDOW_SELECTION" else
    if(length(kinds)==1L && kinds=="CONTINUOUS_DISPLACEMENT") "CONTINUOUS_DISPLACEMENT" else "UNRESOLVED"
  if(!result$converged || any(result$sensitivity$detector_agrees==FALSE)) verdict <- "UNRESOLVED"
  data.frame(target=result$grid$target[1],reference=result$grid$reference[1],
    classification=verdict,sequence=paste(tr$classification,collapse=";"),
    transition_count=nrow(tr),max_bracket_width=if(nrow(tr)) max(tr$width) else NA_real_,
    converged=result$converged,refinement_levels=result$levels,
    weekly_only=TRUE)
}

p25_plot <- function(result,out,target,reference) {
  g <- result$grid;p <- result$peaks;tr <- result$transitions
  draw <- function(filename,expr) {
    grDevices::pdf(file.path(out,filename),width=9,height=7)
    tryCatch(force(expr),finally=grDevices::dev.off())
  }
  draw(sprintf("%s_%s_mechanism.pdf",target,reference),{
    op <- graphics::par(mfrow=c(2,1),mar=c(4,4,2,1))
    for(field in c("position","height")) {
      yr <- range(p[[field]],finite=TRUE)
      graphics::plot(range(g$v),yr,type="n",xlab="Log beta multiplier",
        ylab=if(field=="position") "Season week" else "Expected cases",
        main=paste(target,reference,field))
      if(nrow(tr)) for(i in seq_len(nrow(tr)))
        graphics::rect(tr$v_left[i],yr[1],tr$v_right[i],yr[2],
          col=grDevices::adjustcolor(if(tr$classification[i]=="DOMINANCE_EXCHANGE")
            "firebrick" else "grey60",alpha.f=.2),border=NA)
      for(b in unique(p$branch)) {
        q <- p[p$branch==b & p$genuine,];q <- q[order(q$v),]
        graphics::lines(q$v,q[[field]],col="steelblue")
        graphics::points(q$v,q[[field]],pch=ifelse(q$match_status=="AMBIGUOUS",4,1),cex=.55)
      }
      selected <- p[p$branch %in% g$selected_branch & p$scenario %in% g$scenario,]
      selected <- selected[selected$position==g$selected_week[match(selected$scenario,g$scenario)],]
      graphics::points(selected$v,selected[[field]],pch=16,col="firebrick",cex=.65)
    }
    graphics::par(op)
  })
  draw(sprintf("%s_%s_before_after.pdf",target,reference),{
    if(!nrow(tr)) {graphics::plot.new();return(invisible(NULL))}
    op <- graphics::par(mfrow=c(min(nrow(tr),4L),1),mar=c(4,4,2,1))
    for(i in seq_len(min(nrow(tr),4L))) {
      a <- result$curves[result$curves$scenario==tr$left[i] & result$curves$in_window,]
      b <- result$curves[result$curves$scenario==tr$right[i] & result$curves$in_window,]
      graphics::plot(a$season_week,a$expected,type="l",col="steelblue",
        ylim=range(c(a$expected,b$expected)),xlab="Season week",ylab="Expected cases",
        main=sprintf("%s: %.9f to %.9f",tr$classification[i],tr$v_left[i],tr$v_right[i]))
      graphics::lines(b$season_week,b$expected,col="firebrick")
      graphics::legend("topright",c("before","after"),lty=1,
        col=c("steelblue","firebrick"),bty="n")
    }
    graphics::par(op)
  })
}

p25_audit <- function(cfg,source_cfg,p24cfg,out) {
  checked <- p25_verify(cfg,source_cfg,p24cfg)
  inputs <- p22_inputs(source_cfg);clock <- p24_clock(inputs)
  scores <- p24_scores(inputs,source_cfg,cfg$replay_tolerance)$scores
  if(!all(c("HIST","JOINT","JOINT_COMPETING") %in% names(scores)))
    stop("A required saved reference could not be replayed")
  endpoint <- p25_replay(checked,scores,inputs,clock,cfg)
  write <- function(x,name) utils::write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE)
  write(endpoint,"endpoint_replay")
  runs <- list()
  for(i in seq_len(nrow(checked$targets))) for(ref in c("JOINT","HIST","JOINT_COMPETING")) {
    t <- checked$targets[i,]
    key <- paste(t$target,ref,sep="_")
    runs[[key]] <- p25_refine(scores[[ref]],inputs,clock,ref,t,cfg)
    runs[[key]]$sensitivity <- p25_sensitivity(runs[[key]])
  }
  combine <- function(field) {
    z <- lapply(runs,`[[`,field);z <- Filter(function(x) nrow(x)>0L,z)
    if(length(z)) do.call(rbind,z) else data.frame()
  }
  grid <- combine("grid");peaks <- combine("peaks");gaps <- combine("gaps")
  trans <- combine("transitions");curves <- combine("curves")
  sensitivity <- combine("sensitivity")
  verdicts <- do.call(rbind,lapply(runs,p25_verdict))
  robustness <- do.call(rbind,lapply(split(verdicts,verdicts$target),function(z)
    data.frame(target=z$target[1],all_references_agree=length(unique(z$classification))==1L,
      all_converged=all(z$converged),classifications=paste(paste(z$reference,z$classification,
        sep=":"),collapse=";"))))
  write(grid,"refined_grid");write(peaks,"all_tracked_peaks")
  write(gaps,"paired_peak_height_gaps");write(trans,"transition_brackets")
  write(verdicts,"mechanism_verdicts");write(sensitivity,"numerical_sensitivity")
  write(robustness,"reference_robustness")
  write(curves,"weekly_trajectories")
  write(data.frame(path=names(checked$hashes),md5=unname(checked$hashes)),"artifact_hashes")
  for(key in names(runs)) {
    z <- runs[[key]];parts <- strsplit(key,"_",fixed=TRUE)[[1]]
    p25_plot(z,out,parts[1],paste(parts[-1],collapse="_"))
  }
  lines <- c("# Phase 25 weekly peak audit","",
    "The fitted models, states, forcing, calendar, observation window and simulator are reused from Phase 24.",
    "Results below describe simulated intervention paths, not the historical cause.","",
    "## Reproduction","",
    "Run `Rscript tests/test_phase25.R` and `Rscript scripts/phase25.R --out artifacts/phase25-audit-001` from the repository root.",
    "The output directory must be new. The full experiment is run only by the audit command.","",
    sprintf("Endpoint absolute tolerance: %.1e; minimum bisections: %d; width target: %.1e; stable rounds: %d; maximum bisections: %d.",
      cfg$replay_tolerance,cfg$minimum_bisections,cfg$transition_width,cfg$stable_rounds,cfg$maximum_bisections),
    "A width alone does not establish convergence; the classification signature must also persist.","",
    "## Verdicts",capture.output(print(verdicts,row.names=FALSE)),"",
    "## Reference robustness",capture.output(print(robustness,row.names=FALSE)),"",
    "## Transition brackets",capture.output(print(trans,row.names=FALSE)),"",
    "The weekly SIRS simulator has no integration tolerance or within-week incidence. Branch loss can establish sampled weekly peak loss only; a continuous-time fold cannot be claimed.",
    "Neighboring weekly values for every detected peak are in all_tracked_peaks.csv; complete support trajectories are in weekly_trajectories.csv.",
    "The Phase 24 fitted historical peak shift remains -5 weeks versus the observed -13 weeks. Counterfactual paths do not change that adequacy result.","",
    "## Artifact hashes",capture.output(print(checked$hashes)))
  writeLines(lines,file.path(out,"phase25_report.md"))
  list(status="AUDIT_SAVED",config_version=cfg$version,endpoint=endpoint,
       verdicts=verdicts,hashes=checked$hashes)
}
