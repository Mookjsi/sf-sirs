# Independent analysis of the saved Phase 27 M2 experiments.
m2_read <- function(dir,stage,files) {
  path <- file.path(dir,"result.rds")
  x <- read_phase2_result(path)
  status <- readLines(file.path(dir,"STATUS"),warn=FALSE)
  expected <- c(replay="REPLAY_PASSED",mechanism="MECHANISM_SAVED",
    sensitivity="SENSITIVITY_SAVED")[[stage]]
  if(!identical(x$metadata$phase,27L) || !identical(x$metadata$stage,stage) ||
     !identical(x$result$status,expected) || !identical(status[1L],"STORAGE_COMPLETED"))
    stop("Incomplete or wrong Phase 27 stage: ",dir)
  manifest <- dget(file.path(dir,"manifest.R"))
  if(!identical(manifest$reference_status,"USER_SELECTED_M2") ||
     !identical(manifest$endpoint_fit_md5,"32167fcf5e58f76f8c21d8a013166d1b") ||
     !identical(manifest$endpoint_parameter_hash,"da822ce1b48d242f1585fbf4d1811eec") ||
     !identical(unname(tools::md5sum(manifest$endpoint_fit)),manifest$endpoint_fit_md5) ||
     !identical(unname(tools::md5sum(manifest$config$evaluation)),
       "0faf32bab031b68433258758850c949d"))
    stop("Unexpected M2 endpoint: ",dir)
  tables <- lapply(files,function(f) utils::read.csv(file.path(dir,paste0(f,".csv")),
    stringsAsFactors=FALSE))
  names(tables) <- files
  list(metadata=x$metadata,result=x$result,manifest=manifest,tables=tables,dir=dir)
}
m2_validate_columns <- function(x,columns,name) {
  missing <- setdiff(columns,names(x))
  if(length(missing)) stop(name," missing columns: ",paste(missing,collapse=", "))
  invisible(TRUE)
}
m2_cell_signature <- function(stages,stage) {
  files <- c("R/phase27.R","R/m2_hypotheses.R","R/sirs_core.R",
    "config/phase27.R","config/phase26.R")
  list(schema="m2-full-season-cell-1",stage=stage,
    fit_hash=stages$replay$manifest$endpoint_fit_md5,
    parameter_hash=stages$replay$manifest$endpoint_parameter_hash,
    evaluation_hash=unname(tools::md5sum(stages$replay$manifest$config$evaluation)),
    input_hashes=stages$mechanism$manifest$input_hashes,
    source_hashes=tools::md5sum(files),config=phase27_config())
}
m2_cell_read <- function(path,signature) {
  if(!file.exists(path)) return(NULL)
  x <- readRDS(path)
  compatible <- identical(x$signature,signature)
  previous_modules <- if(identical(signature$stage,"H1_MATCHED_BETA_ONLY"))
    "53014353ef72a6107dc015064502038b" else if(identical(signature$stage,
      "H3_STATE_MEMORY")) c("050f030d2c46bcb8f93673741be7a3c9",
        "247d16b637799e4ebcfce6995af38f48") else character()
  previous_module <- unname(x$signature$source_hashes["R/m2_hypotheses.R"])
  if(!compatible && length(previous_module)==1L &&
     previous_module %in% previous_modules) {
    prior <- signature
    prior$source_hashes["R/m2_hypotheses.R"] <- previous_module
    compatible <- identical(x$signature,prior)
    if(compatible) message("Reusing numerical checkpoint after reporting fix: ",path)
  }
  if(!compatible) stop("Incompatible completed cell: ",path)
  x$value
}
m2_cell_write <- function(path,signature,value) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- paste0(path,".partial")
  saveRDS(list(signature=signature,value=value),temporary)
  if(!file.rename(temporary,path)) stop("Could not save cell: ",path)
  invisible(path)
}
m2_inventory <- function(paths,out) {
  specs <- list(replay=c("replay_checks"),mechanism=c("peak_path_grid",
    "peak_branches","paired_height_gaps","mechanism_intervals","peak_curves"),
    sensitivity=c("paired_sensitivity","sensitivity_variants","joint_comparison"))
  stages <- lapply(names(specs),function(s) m2_read(paths[[s]],s,specs[[s]]))
  names(stages) <- names(specs)
  hashes <- vapply(stages,function(s) s$manifest$endpoint_parameter_hash,"")
  if(length(unique(hashes))!=1L) stop("Phase 27 endpoint hashes differ")
  checks <- stages$replay$tables$replay_checks
  m2_validate_columns(checks,c("passed","error"),"replay_checks")
  if(!all(checks$passed) || any(!is.finite(checks$error))) stop("Saved replay checks failed")
  m <- stages$mechanism$tables;s <- stages$sensitivity$tables
  m2_validate_columns(m$peak_path_grid,c("scenario","target","policy","path","v",
    "selected_week","selected_branch","right_support"),"peak_path_grid")
  m2_validate_columns(m$mechanism_intervals,c("policy","season","path","v_left","v_right",
    "left","right","classification"),"mechanism_intervals")
  m2_validate_columns(m$peak_branches,c("scenario","position","height","genuine","support"),"peak_branches")
  m2_validate_columns(m$paired_height_gaps,c("scenario","signed_gap","trackable"),"paired_height_gaps")
  m2_validate_columns(m$peak_curves,c("scenario","position","expected"),"peak_curves")
  m2_validate_columns(s$paired_sensitivity,c("reference","policy","design","year","h",
    "paired_eligible","exclusion","E_S","E_beta"),"paired_sensitivity")
  m2_validate_columns(s$sensitivity_variants,c("reference","policy","design","year","h",
    "axis","sign","status","selected_week","quadratic_time"),"sensitivity_variants")
  code <- "R/phase27.R";current <- unname(tools::md5sum(code))
  archived <- unname(stages$mechanism$manifest$code_hashes[
    grepl("/R/phase27[.]R$",names(stages$mechanism$manifest$code_hashes))])
  if(length(archived)!=1L) stop("Archived Phase 27 source hash is ambiguous")
  archived_all <- vapply(stages,function(a) {
    h <- unname(a$manifest$code_hashes[grepl("/R/phase27[.]R$",
      names(a$manifest$code_hashes))])
    if(length(h)!=1L) stop("Archived source hash missing")
    h
  },"")
  if(length(unique(archived_all))!=1L) stop("Saved stages used different Phase 27 source hashes")
  inv <- do.call(rbind,lapply(names(stages),function(k) {
    a <- stages[[k]]
    data.frame(stage=k,directory=normalizePath(a$dir),status=a$result$status,
      reference=a$manifest$reference_status,fit_md5=a$manifest$endpoint_fit_md5,
      parameter_hash=a$manifest$endpoint_parameter_hash,
      evaluation_md5=unname(tools::md5sum(a$manifest$config$evaluation)),
      code_commit=a$manifest$code_commit,archived_phase27_md5=archived,
      current_phase27_md5=current,source_match=identical(current,archived),
      evidence_scope=if(identical(current,archived)) "CURRENT_CODE" else "HISTORICAL_CODE",
      stringsAsFactors=FALSE)
  }))
  inv$selected <- TRUE
  candidates <- unlist(lapply(names(stages),function(k)
    Sys.glob(file.path("artifacts",paste0("phase27-",k,"-*")))))
  candidates <- setdiff(normalizePath(candidates,mustWork=FALSE),inv$directory)
  extra <- lapply(candidates,function(dir) {
    manifest_path <- file.path(dir,"manifest.R")
    if(!file.exists(manifest_path)) return(NULL)
    a <- tryCatch(dget(manifest_path),error=function(e) NULL)
    if(is.null(a) || is.null(a$stage)) return(NULL)
    code_hash <- unname(a$code_hashes[grepl("/R/phase27[.]R$",names(a$code_hashes))])
    if(length(code_hash)!=1L) code_hash <- NA_character_
    data.frame(stage=a$stage,directory=dir,
      status=if(file.exists(file.path(dir,"STATUS")))
        readLines(file.path(dir,"STATUS"),n=2L,warn=FALSE)[1L] else "MISSING",
      reference=a$reference_status,fit_md5=a$endpoint_fit_md5,
      parameter_hash=a$endpoint_parameter_hash,
      evaluation_md5=if(file.exists(a$config$evaluation))
        unname(tools::md5sum(a$config$evaluation)) else NA_character_,
      code_commit=a$code_commit,archived_phase27_md5=code_hash,
      current_phase27_md5=current,source_match=identical(current,code_hash),
      evidence_scope="UNSELECTED_RUN",selected=FALSE,stringsAsFactors=FALSE)
  })
  extra <- Filter(Negate(is.null),extra)
  if(length(extra)) inv <- rbind(inv,do.call(rbind,extra))
  p27_write(inv,out,"existing_experiment_inventory")
  current_sources <- c("R/phase27.R","R/phase25.R","R/sirs_core.R",
    "R/m2_hypotheses.R","scripts/m2_hypotheses.R")
  p27_write(data.frame(path=current_sources,
    current_md5=unname(tools::md5sum(current_sources))),out,"current_source_hashes")
  p27_write(data.frame(component=c("recurrence","intervention_schedule",
    "peak_matching","interpolation","classification","other_stages"),
    archived_source_available=FALSE,current_code_md5=current,
    archived_code_md5=archived,
    audit_status="HISTORICAL_SOURCE_UNAVAILABLE_FOR_EXACT_DIFF",
    stringsAsFactors=FALSE),out,"source_version_audit")
  commit <- suppressWarnings(system2("git",c("rev-parse","HEAD"),stdout=TRUE))
  dirty <- suppressWarnings(system2("git",c("status","--porcelain"),stdout=TRUE))
  source_files <- c(current_sources,"config/phase27.R","config/phase26.R")
  provenance <- data.frame(schema_version="m2-full-season-1",
    code_commit=commit[1L],dirty=length(dirty)>0L,
    file=source_files,md5=unname(tools::md5sum(source_files)),
    fit_md5=stages$replay$manifest$endpoint_fit_md5,
    parameter_hash=stages$replay$manifest$endpoint_parameter_hash,
    evaluation_md5=unname(tools::md5sum(stages$replay$manifest$config$evaluation)))
  p27_write(provenance,out,"run_provenance")
  stages
}
m2_reclassify <- function(e,g,p,curves) {
  a <- g[match(e$left,g$scenario),];b <- g[match(e$right,g$scenario),]
  if(anyNA(a$scenario) || anyNA(b$scenario)) stop("Mechanism edge has missing scenario")
  revised <- e$classification;reason <- rep("SAVED_CLASSIFICATION_RETAINED",nrow(e))
  for(i in seq_len(nrow(e))) {
    if(!a$right_support[i] || !b$right_support[i]) {
      revised[i] <- "UNRESOLVED_SUPPORT";reason[i] <- "MISSING_CHRONOLOGICAL_RIGHT_SUPPORT"
    } else if(a$selected_week[i] %in% c(1L,a$window_length[i]) ||
               b$selected_week[i] %in% c(1L,b$window_length[i])) {
      revised[i] <- "BOUNDARY_CANDIDATE";reason[i] <- "SELECTED_WEEK_AT_WINDOW_EDGE"
      extended_winner <- function(scenario,window) {
        q <- p[p$scenario==scenario & p$genuine & p$support &
          p$position>=1L,]
        if(!nrow(q)) return(NA)
        winner <- q[q$height==max(q$height),]
        if(nrow(winner)!=1L) return(NA)
        winner$position>window
      }
      wa <- extended_winner(a$scenario[i],a$window_length[i])
      wb <- extended_winner(b$scenario[i],b$window_length[i])
      if(!is.na(wa) && !is.na(wb) && wa!=wb &&
         a$selected_week[i]!=b$selected_week[i]) {
        revised[i] <- "WINDOW_SELECTION_DEMONSTRATED"
        reason[i] <- "SUPPORTED_EXTENDED_WINNER_CROSSES_WINDOW"
      }
    } else if(e$classification[i]=="DOMINANCE_EXCHANGE" &&
              (!is.finite(e$early_minus_late_left[i]) ||
               !is.finite(e$early_minus_late_right[i]) ||
               e$early_minus_late_left[i]*e$early_minus_late_right[i]>=0 ||
               is.na(e$left_branch[i]) || is.na(e$right_branch[i]) ||
               e$left_branch[i]==e$right_branch[i])) {
      revised[i] <- "UNRESOLVED_MATCHING"
      reason[i] <- "EXCHANGE_EVIDENCE_INCOMPLETE"
    } else if(e$left_count[i]!=e$right_count[i]) {
      revised[i] <- if(e$left_count[i]>e$right_count[i]) "SAMPLED_PEAK_LOSS" else "SAMPLED_PEAK_BIRTH"
      reason[i] <- "WEEKLY_LOCAL_MAXIMUM_COUNT_CHANGE"
    } else if(e$classification[i]=="SAMPLED_PEAK_LOSS") {
      revised[i] <- "UNRESOLVED_MATCHING";reason[i] <- "NO_PEAK_COUNT_CHANGE"
    }
  }
  e$old_classification <- e$classification;e$classification <- revised;e$reason <- reason
  e$timing_change <- a$selected_week!=b$selected_week
  e$absolute_movement <- abs(a$selected_week-b$selected_week)
  e$actual_window <- ifelse(e$season<=2019L,"EPIDEMIC_SEASON",
    ifelse(e$path=="BETA_ONLY","CALENDAR_YEAR","EPIDEMIC_SEASON"))
  e$origin <- ifelse(e$path=="BETA_ONLY","2013_ORIGIN","MATCHED_SEASON_ORIGIN")
  gap_scenarios <- unique(curves$scenario)
  e$curve_support <- e$left %in% gap_scenarios & e$right %in% gap_scenarios
  e
}
m2_parent_index <- function(left,grid) {
  pmax(1L,pmin(length(grid)-1L,findInterval(left,grid,rightmost.closed=TRUE)))
}
m2_h1 <- function(stages,out) {
  t <- stages$mechanism$tables
  e <- m2_reclassify(t$mechanism_intervals,t$peak_path_grid,t$peak_branches,t$peak_curves)
  p27_write(e,out,"classification_audit")
  # Each parent is an adjacent pair in the declared unrefined mechanism grid.
  base <- sort(unique(p22_mechanism_grid(list(mechanism_grid=phase27_config()$mechanism_grid))))
  target <- p25_verify(phase25_config(),phase22_config(),phase24_config())$targets
  grids <- list(BETA_ONLY=base)
  for(path in c("MATCHED_BETA_ONLY","MATCHED_JOINT")) {
    target_row <- target[target$path==if(path=="MATCHED_BETA_ONLY") "BETA_ONLY" else "JOINT",]
    grids[[path]] <- sort(unique(c(base,-base,target_row$x_left,target_row$x_right)))
  }
  e$parent <- vapply(seq_len(nrow(e)),function(i) {
    b <- grids[[e$path[i]]]
    m2_parent_index(e$v_left[i],b)
  },integer(1))
  keys <- unique(e[c("policy","season","path","parent")])
  rows <- lapply(seq_len(nrow(keys)),function(i) {
    k <- keys[i,];q <- e[e$policy==k$policy & e$season==k$season &
      e$path==k$path & e$parent==k$parent,]
    changing <- q[q$timing_change,,drop=FALSE]
    known <- changing$classification %in% c("CONTINUOUS_DISPLACEMENT","DOMINANCE_EXCHANGE")
    b <- grids[[k$path]]
    data.frame(k,parent_left=b[k$parent],parent_right=b[k$parent+1L],
      timing_event=nrow(changing)>0L,displacement=any(changing$classification=="CONTINUOUS_DISPLACEMENT"),
      exchange=any(changing$classification=="DOMINANCE_EXCHANGE"),
      unresolved=any(!known),mixed=length(unique(changing$classification))>1L,
      absolute_movement=sum(changing$absolute_movement),
      displacement_movement=sum(changing$absolute_movement[
        changing$classification=="CONTINUOUS_DISPLACEMENT"]),
      exchange_movement=sum(changing$absolute_movement[
        changing$classification=="DOMINANCE_EXCHANGE"]),
      stringsAsFactors=FALSE)
  })
  parents <- do.call(rbind,rows)
  p27_write(parents,out,"h1_event_summary")
  coverage <- aggregate(t$paired_height_gaps$scenario,
    t$paired_height_gaps[c("policy","season","path")],length)
  names(coverage)[4L] <- "exact_early_late_pairs"
  total <- aggregate(t$peak_path_grid$scenario,
    data.frame(policy=t$peak_path_grid$policy,
      season=as.integer(t$peak_path_grid$target),path=t$peak_path_grid$path),length)
  names(total)[4L] <- "scenarios"
  coverage <- merge(total,coverage,by=c("policy","season","path"),all.x=TRUE)
  coverage$exact_early_late_pairs[is.na(coverage$exact_early_late_pairs)] <- 0L
  coverage$fraction <- coverage$exact_early_late_pairs/coverage$scenarios
  p27_write(coverage,out,"h1_gap_coverage")
  list(edges=e,parents=parents,coverage=coverage)
}
m2_h2 <- function(stages,out) {
  x <- stages$sensitivity$tables$paired_sensitivity
  v <- stages$sensitivity$tables$sensitivity_variants
  x$magnitude_beta_gt_S <- x$paired_eligible & abs(x$E_beta)>abs(x$E_S)
  x$sign_beta_label <- ifelse(!x$paired_eligible,"INELIGIBLE",ifelse(x$E_beta<0,"NEGATIVE",ifelse(x$E_beta>0,"POSITIVE","ZERO")))
  x$sign_S_label <- ifelse(!x$paired_eligible,"INELIGIBLE",ifelse(x$E_S<0,"NEGATIVE",ifelse(x$E_S>0,"POSITIVE","ZERO")))
  p27_write(x,out,"h2_magnitude_summary")
  p27_write(x[c("reference","policy","design","year","h","paired_eligible",
    "exclusion","sign_beta_label","sign_S_label","E_beta","E_S")],out,"h2_direction_summary")
  keys <- c("reference","policy","design","year","h")
  rows <- lapply(seq_len(nrow(x)),function(i) {
    a <- x[i,];q <- v[v$reference==a$reference & v$policy==a$policy &
      v$design==a$design & v$year==a$year & abs(v$h-a$h)<1e-12,]
    qS <- q[q$axis=="S",];qB <- q[q$axis=="beta",]
    data.frame(a[keys],status=if(a$paired_eligible) "ELIGIBLE" else a$exclusion,
      S_stencil_changed=nrow(qS)==2L && length(unique(qS$selected_week))>1L,
      beta_stencil_changed=nrow(qB)==2L && length(unique(qB$selected_week))>1L,
      discrete_S=a$discrete_S,discrete_beta=a$discrete_beta,
      centroid_S=a$cog_S,centroid_beta=a$cog_beta)
  })
  audit <- do.call(rbind,rows)
  for(i in seq_len(nrow(audit))) {
    q <- x[x$reference==audit$reference[i] & x$policy==audit$policy[i] &
      x$design==audit$design[i] & x$year==audit$year[i] & x$paired_eligible,]
    audit$step_range_S[i] <- if(nrow(q)>1L) diff(range(q$E_S)) else NA_real_
    audit$step_range_beta[i] <- if(nrow(q)>1L) diff(range(q$E_beta)) else NA_real_
  }
  p27_write(audit,out,"h2_step_audit")
  matched <- x[x$reference=="M2" & x$design=="MATCHED_SEASON_ORIGIN",]
  selected <- do.call(rbind,lapply(2013:2022,function(y)
    do.call(rbind,lapply(c("raw","blend"),function(policy) {
      q <- matched[matched$year==y & matched$policy==policy,]
      eligible <- q[q$paired_eligible & is.finite(q$E_S) & is.finite(q$E_beta),]
      if(nrow(eligible)) {
        eligible <- eligible[order(eligible$h),]
        chosen <- eligible[(nrow(eligible)+1L)%/%2L,]
      } else chosen <- q[1L,]
      data.frame(season_id=sprintf("%d/%02d",y,(y+1L)%%100L),policy=policy,
        tested_steps=nrow(q),eligible_steps=nrow(eligible),selected_h=if(nrow(eligible)) chosen$h else NA_real_,
        selected_E_S=if(nrow(eligible)) chosen$E_S else NA_real_,
        selected_E_beta=if(nrow(eligible)) chosen$E_beta else NA_real_,
        magnitude_order=if(nrow(eligible)) if(abs(chosen$E_beta)>abs(chosen$E_S))
          "BETA_LARGER" else "S_LARGER" else "INELIGIBLE",
        direction=if(nrow(eligible)) if(sign(chosen$E_beta)==sign(chosen$E_S))
          if(chosen$E_beta<0) "BOTH_ADVANCE" else "BOTH_DELAY" else "OPPOSITE"
          else "INELIGIBLE",
        selection_rule="MEDIAN_ELIGIBLE_LOG_STEP",
        step_range_S=if(nrow(eligible)>1L) diff(range(eligible$E_S)) else NA_real_,
        step_range_beta=if(nrow(eligible)>1L) diff(range(eligible$E_beta)) else NA_real_)
    }))))
  p27_write(selected,out,"h2_season_summary")
  list(estimates=x,audit=audit)
}
m2_season_table <- function(clock) {
  rows <- lapply(2013:2022,function(y) {
    ix <- p27_window(clock,y,TRUE)
    at <- ix[1L]-1L
    dates <- as.Date(mapply(p24_date,clock$year[ix],clock$week[ix]))
    data.frame(season_id=sprintf("%d/%02d",y,(y+1L)%%100L),
      year=y,start_date=p24_date(clock$year[ix[1L]],clock$week[ix[1L]]),
      end_date=p24_date(clock$year[tail(ix,1L)],clock$week[tail(ix,1L)]),
      date_source="ISO_WEEK_DERIVED_FROM_PANEL_YEAR_WEEK",
      observations=length(ix),week_53=any(clock$week[ix]==53L),
      omitted_week_53=any(as.integer(diff(dates))==14L),
      restart_state_date=p24_date(clock$year[at],clock$week[at]),
      first_affected_transition_date=p24_date(clock$year[ix[1L]],clock$week[ix[1L]]),
      left_support=at>=1L,right_support=tail(ix,1L)<nrow(clock),
      window_status="COMPLETE",stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}
m2_plan <- function(stages,out) {
  inputs <- p22_inputs(phase22_config());clock <- p26_clock(inputs)
  season <- m2_season_table(clock);p27_write(season,out,"season_coverage")
  mechanism <- stages$mechanism$tables$peak_path_grid
  sensitivity <- stages$sensitivity$tables$paired_sensitivity
  state_dirs <- Sys.glob("artifacts/m2-hypotheses-state-*")
  state_dirs <- state_dirs[vapply(state_dirs,function(d)
    file.exists(file.path(d,"STATUS")) &&
      identical(readLines(file.path(d,"STATUS"),n=1L),"STORAGE_COMPLETED") &&
      file.exists(file.path(d,"state_memory.csv")),logical(1))]
  state <- if(length(state_dirs)) do.call(rbind,lapply(state_dirs,function(d)
    utils::read.csv(file.path(d,"state_memory.csv")))) else data.frame()
  rows <- lapply(2013:2022,function(y) lapply(c("raw","blend"),function(policy) {
    m <- mechanism[mechanism$path=="MATCHED_BETA_ONLY" &
      as.integer(mechanism$target)==y & mechanism$policy==policy,]
    h <- sensitivity[sensitivity$reference=="M2" &
      sensitivity$design=="MATCHED_SEASON_ORIGIN" & sensitivity$year==y &
      sensitivity$policy==policy,]
    scan <- if(nrow(state)) state[state$policy==policy & state$year==y,] else state
    local <- if(nrow(scan)) any(grepl("/Iminus/",scan$scenario,fixed=TRUE)) else FALSE
    path <- if(nrow(scan)) any(grepl("/path/",scan$scenario,fixed=TRUE)) else FALSE
    data.frame(season_id=season$season_id[season$year==y],policy=policy,
      h1_matched_beta_only=if(nrow(m)) "COMPLETED_HISTORICAL_SIMULATION" else "IMPLEMENTED_NOT_RUN",
      h2_matched_S_beta=if(nrow(h)==length(phase27_config()$sensitivity_steps))
        if(any(h$paired_eligible)) "COMPLETED" else "COMPLETED_NUMERICALLY_INELIGIBLE"
        else "IMPLEMENTED_NOT_RUN",
      h3_local_I=if(local) "COMPLETED_HISTORICAL_SIMULATION" else "IMPLEMENTED_NOT_RUN",
      h3_finite_I_path=if(path) "COMPLETED_HISTORICAL_SIMULATION" else "IMPLEMENTED_NOT_RUN")
  }))
  coverage <- do.call(rbind,unlist(rows,recursive=FALSE))
  p27_write(coverage,out,"coverage_matrix")
  designs <- list(
    list(name="H1_MATCHED_BETA_ONLY",column="h1_matched_beta_only",
      perturbation="effective_beta_in_season * exp(v)"),
    list(name="H2_MATCHED_S_BETA",column="h2_matched_S_beta",
      perturbation="initial_S * exp(u), R compensation; effective_beta * exp(v)"),
    list(name="H3_LOCAL_I",column="h3_local_I",
      perturbation="initial_I * exp(z), saved-mass R compensation"),
    list(name="H3_FINITE_I_PATH",column="h3_finite_I_path",
      perturbation="initial_I * exp(z), z in [-log(1.2),log(1.2)]"))
  inventory <- do.call(rbind,lapply(designs,function(d) data.frame(
    experiment=d$name,season_id=coverage$season_id,policy=coverage$policy,
    origin="SAVED_PRESEASON_STATE",window="WEEK_27_THROUGH_FOLLOWING_WEEK_26",
    perturbation=d$perturbation,status=coverage[[d$column]],
    simulation_provenance=ifelse(grepl("HISTORICAL",coverage[[d$column]]),
      "HISTORICAL_EXPERIMENT",ifelse(coverage[[d$column]]=="COMPLETED",
        "HISTORICAL_PHASE27_H2","NEW_SIMULATION_REQUIRED")))))
  p27_write(inventory,out,"design_inventory")
  writeLines(c("# M2 full-season execution plan",
    "H2 matched-season cells are reused from the completed Phase 27 sensitivity artifact.",
    "The historical 2017/18 H1 numerical cell may be reused after grid and window checks.",
    "The new H1 command fills remaining seasons and retains its historical simulation source.",
    "H3 prior partial outputs are inventory evidence; run the full-season command for a common grid.",
    paste("Missing H1 season-policy cells:",sum(coverage$h1_matched_beta_only=="IMPLEMENTED_NOT_RUN")),
    paste("Missing H2 season-policy cells:",sum(coverage$h2_matched_S_beta=="IMPLEMENTED_NOT_RUN")),
    paste("Missing H3 local season-policy cells:",sum(coverage$h3_local_I=="IMPLEMENTED_NOT_RUN"))),
    file.path(out,"plan.md"))
  invisible(coverage)
}
m2_h1_full <- function(stages,out,progress=function(...) invisible(NULL)) {
  c22 <- phase22_config();c26 <- phase26_config();c27 <- phase27_config()
  inputs <- p22_inputs(c22);ref <- p27_reference(inputs,c22,c26,c27)
  clock <- p26_clock(inputs);season <- m2_season_table(clock)
  p27_write(season,out,"season_coverage")
  base <- p22_mechanism_grid(list(mechanism_grid=c27$mechanism_grid))
  target <- p25_verify(phase25_config(),c22,phase24_config())$targets
  target <- target[target$path=="BETA_ONLY",]
  vals0 <- sort(unique(c(base,-base,target$x_left,target$x_right)))
  old <- stages$mechanism$tables
  archived_config <- stages$mechanism$manifest$config
  knobs <- c("refinement_width","refinement_rounds","refinement_cap","padding_weeks")
  compatible_config <- all(vapply(knobs,function(k) identical(archived_config[[k]],c27[[k]]),logical(1)))
  grid_hash <- stages$mechanism$manifest$input_hashes[c27$mechanism_grid]
  compatible_grid <- length(grid_hash)==1L && !is.na(grid_hash) &&
    identical(unname(grid_hash),unname(tools::md5sum(c27$mechanism_grid)))
  all <- list()
  signature <- m2_cell_signature(stages,"H1_MATCHED_BETA_ONLY")
  for(policy in c("raw","blend")) {
    score <- if(policy=="raw") ref$score else p26_score(ref$fit$endpoint,
      inputs,p22_problem(inputs),c26,"M2",ref$fit$penalty,blend=TRUE)
    for(y in 2013:2022) {
      progress("H1",paste(policy,y))
      cell_path <- file.path(out,"cells",sprintf("h1-%s-%d.rds",policy,y))
      cached <- m2_cell_read(cell_path,signature)
      if(!is.null(cached)) {
        all[[paste(policy,y)]] <- cached
        progress("H1_REUSED_CHECKPOINT",paste(policy,y))
        next
      }
      saved <- old$peak_path_grid[old$peak_path_grid$policy==policy &
        old$peak_path_grid$path=="MATCHED_BETA_ONLY" &
        as.integer(old$peak_path_grid$target)==y,,drop=FALSE]
      if(nrow(saved) && compatible_config && compatible_grid &&
         all(vapply(vals0,function(v) any(abs(saved$v-v)<1e-12),logical(1))) &&
         all(saved$window_length==season$observations[season$year==y])) {
        key <- function(t) t[t$policy==policy & t$path=="MATCHED_BETA_ONLY" &
          as.integer(if("target" %in% names(t)) t$target else t$season)==y,,drop=FALSE]
        all[[paste(policy,y)]] <- list(grid=saved,peaks=key(old$peak_branches),
          curves=key(old$peak_curves),edges=key(old$mechanism_intervals),
          gaps=key(old$paired_height_gaps),policy=policy,source="HISTORICAL_SIMULATION")
        m2_cell_write(cell_path,signature,all[[paste(policy,y)]])
        next
      }
      vals <- vals0;cache <- list();round <- 0L
      repeat {
        missing <- vals[!vapply(vals,function(v) any(vapply(cache,function(r)
          abs(r$grid$v-v)<1e-14,logical(1))),logical(1))]
        for(v in missing) {
          id <- sprintf("%s-MATCHED_BETA_ONLY-%d-%016.12f",policy,y,v)
          run <- p27_state_profile(score,inputs,y,v=v)
          r <- p27_peak_record(run$expected,clock,y,id,c27$padding_weeks,TRUE)
          r$grid$v <- v;r$grid$path <- "MATCHED_BETA_ONLY"
          if(nrow(r$peaks)) {r$peaks$v <- v;r$peaks$path <- "MATCHED_BETA_ONLY"}
          r$curve$v <- v;r$curve$path <- "MATCHED_BETA_ONLY"
          cache[[id]] <- r
        }
        tracked <- p25_track(cache,phase25_config())
        gaps <- p27_height_gaps(tracked$grid,tracked$peaks,phase25_config())
        edges <- p27_classify(tracked$grid,tracked$peaks,gaps)
        trigger <- edges$classification %in% c("DOMINANCE_EXCHANGE","UNRESOLVED",
          "SAMPLED_PEAK_LOSS","WINDOW_SELECTION") |
          edges$left_week!=edges$right_week | edges$left_count!=edges$right_count |
          (is.finite(edges$early_minus_late_left) &
           is.finite(edges$early_minus_late_right) &
           edges$early_minus_late_left*edges$early_minus_late_right<=0)
        candidate <- (edges$v_left+edges$v_right)/2
        candidate <- setdiff(candidate[trigger & edges$width>c27$refinement_width],vals)
        round <- round+1L
        if(!length(candidate) || round>c27$refinement_rounds ||
           length(vals)>=c27$refinement_cap) break
        vals <- sort(c(vals,head(candidate,c27$refinement_cap-length(vals))))
      }
      all[[paste(policy,y)]] <- list(grid=tracked$grid,peaks=tracked$peaks,
        curves=p27_bind(lapply(cache,`[[`,"curve")),edges=edges,gaps=gaps,
        policy=policy,source="NEW_SIMULATION")
      m2_cell_write(cell_path,signature,all[[paste(policy,y)]])
    }
  }
  table <- function(key) {
    rows <- lapply(all,function(a) {
      q <- a[[key]]
      if(!nrow(q)) return(NULL)
      if(!"u" %in% names(q) && key=="curves") q$u <- 0
      q$policy <- a$policy
      q$simulation_source <- a$source
      q
    })
    rows <- Filter(Negate(is.null),rows)
    columns <- unique(unlist(lapply(rows,names),use.names=FALSE))
    rows <- lapply(rows,function(q) {
      for(column in setdiff(columns,names(q))) q[[column]] <- NA
      q[columns]
    })
    p27_bind(rows)
  }
  tables <- list(peak_path_grid=table("grid"),peak_branches=table("peaks"),
    peak_curves=table("curves"),mechanism_intervals=table("edges"),
    paired_height_gaps=table("gaps"))
  for(name in names(tables)) p27_write(tables[[name]],out,name)
  stages$mechanism$tables <- tables
  result <- m2_h1(stages,out)
  parents <- result$parents
  edge_movement <- aggregate(result$edges$absolute_movement,
    result$edges[c("policy","season","path","classification")],sum)
  names(edge_movement)[5L] <- "absolute_movement"
  p27_write(edge_movement,out,"h1_movement_by_mechanism")
  grid <- tables$peak_path_grid
  summaries <- lapply(seq_len(nrow(season)*2L),function(i) {
    y <- season$year[(i+1L)%/%2L];policy <- c("raw","blend")[(i-1L)%%2L+1L]
    q <- parents[parents$season==y & parents$policy==policy &
      parents$path=="MATCHED_BETA_ONLY",,drop=FALSE]
    g <- grid[as.integer(grid$target)==y & grid$policy==policy &
      grid$path=="MATCHED_BETA_ONLY",,drop=FALSE]
    g <- g[order(g$v),]
    identified <- q$timing_event & !q$unresolved & !q$mixed
    data.frame(season_id=season$season_id[season$year==y],policy=policy,
      shift_events=sum(q$timing_event),displacement_events=sum(identified & q$displacement),
      exchange_events=sum(identified & q$exchange),mixed_events=sum(q$mixed),
      unresolved_events=sum(q$unresolved),no_shift=sum(q$timing_event)==0L,
      displacement_movement=sum(q$displacement_movement[identified]),
      exchange_movement=sum(q$exchange_movement[identified]),
      mixed_parent_movement=sum(q$absolute_movement[q$mixed]),
      unresolved_parent_movement=sum(q$absolute_movement[q$unresolved]),
      total_absolute_movement=sum(q$absolute_movement),
      signed_endpoint_movement=tail(g$selected_week,1L)-g$selected_week[1L],
      identified_domain=sum(q$parent_right[!q$unresolved]-q$parent_left[!q$unresolved]),
      unresolved_domain=sum(q$parent_right[q$unresolved]-q$parent_left[q$unresolved]),
      displacement_event_fraction=if(sum(identified))
        sum(identified & q$displacement)/sum(identified) else NA_real_)
  })
  summary <- do.call(rbind,summaries)
  p27_write(summary,out,"h1_season_summary")
  equal <- aggregate(summary$displacement_event_fraction,
    list(policy=summary$policy),function(x) if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE))
  names(equal)[2L] <- "equal_season_displacement_fraction"
  p27_write(equal,out,"h1_equal_season_summary")
  invisible(result)
}
m2_state_case <- function(score,inputs,clock,year,policy,id,u=0,z=0,
                          donor_state=NULL,padding=13L) {
  run <- p27_state_profile(score,inputs,year,u,z,donor_state=donor_state)
  peak <- p27_peak_for_sensitivity(run$expected,clock,year,id,season=TRUE)
  record <- p27_peak_record(run$expected,clock,year,id,padding,season=TRUE)
  ix <- run$window;w <- run$expected[ix]
  trajectory <- run$states
  trajectory$global_row <- run$restart+trajectory$row-1L
  trajectory <- trajectory[trajectory$global_row %in% ix,]
  data <- data.frame(policy=policy,year=year,scenario=id,u=u,z=z,
    S_restart=run$initial_state["S"],I_restart=run$initial_state["I"],
    R_restart=run$initial_state["R"],saved_mass=run$saved_mass,
    peak_week=record$grid$selected_week,quadratic_time=peak$time,
    peak_status=peak$status,peak_height=max(w),
    centroid=sum(seq_along(w)*w)/sum(w),burden=sum(w),
    stringsAsFactors=FALSE)
  list(data=data,record=record,trajectory=trajectory,peak=peak)
}
m2_interaction_row <- function(state,policy,year) {
  prefix <- paste(policy,year,sep="/")
  ids <- paste0(prefix,"/",c("baseline","swap/S","swap/I","swap/SI"))
  q <- state[match(ids,state$scenario),,drop=FALSE]
  valid <- !anyNA(q$scenario) && all(q$peak_status=="ELIGIBLE") &&
    all(is.finite(q$quadratic_time))
  data.frame(policy=policy,year=year,
    donor=if(year==2015L) 2017L else 2015L,
    delta_S=q$S_restart[2L]-q$S_restart[1L],
    delta_I=q$I_restart[3L]-q$I_restart[1L],
    Delta_S=if(valid) q$quadratic_time[2L]-q$quadratic_time[1L] else NA_real_,
    Delta_I=if(valid) q$quadratic_time[3L]-q$quadratic_time[1L] else NA_real_,
    interaction=if(valid) q$quadratic_time[4L]-q$quadratic_time[2L]-
      q$quadratic_time[3L]+q$quadratic_time[1L] else NA_real_,
    reason=if(valid) "" else paste(unique(c(if(anyNA(q$scenario)) "MISSING_CELL",
      q$peak_status[!is.na(q$peak_status) & q$peak_status!="ELIGIBLE"])),collapse=";"))
}
m2_state_memory <- function(stages,out,years=2013:2022,
                            progress=function(...) invisible(NULL)) {
  c22 <- phase22_config();c26 <- phase26_config();c27 <- phase27_config()
  inputs <- p22_inputs(c22);ref <- p27_reference(inputs,c22,c26,c27)
  clock <- p26_clock(inputs);cases <- list();runs <- list();scans <- list();refinement_failures <- list()
  signature <- m2_cell_signature(stages,"H3_STATE_MEMORY")
  for(policy in c("raw","blend")) {
    score <- if(policy=="raw") ref$score else p26_score(ref$fit$endpoint,
      inputs,p22_problem(inputs),c26,"M2",ref$fit$penalty,blend=TRUE)
    add <- function(year,id,u=0,z=0,donor_state=NULL) {
      r <- tryCatch(m2_state_case(score,inputs,clock,year,policy,id,u,z,donor_state),
        error=function(e) {
          if(!grepl("Infeasible matched state",conditionMessage(e),fixed=TRUE)) stop(e)
          at <- p27_window(clock,year,TRUE)[1L]-1L
          saved <- unlist(score$states[at,c("S","I","R")],use.names=TRUE)
          candidate <- if(is.null(donor_state))
            c(S=unname(saved["S"])*exp(u),I=unname(saved["I"])*exp(z)) else donor_state
          list(data=data.frame(policy=policy,year=year,scenario=id,u=u,z=z,
            S_restart=unname(candidate["S"]),I_restart=unname(candidate["I"]),
            R_restart=sum(saved)-sum(candidate),saved_mass=sum(saved),
            peak_week=NA_real_,quadratic_time=NA_real_,
            peak_status=paste0("INFEASIBLE: ",conditionMessage(e)),
            peak_height=NA_real_,centroid=NA_real_,burden=NA_real_),
            peak=list(status="INFEASIBLE",time=NA_real_),record=NULL,
            trajectory=NULL)
        })
      if(!is.null(r$record)) runs[[id]] <<- r
      cases[[length(cases)+1L]] <<- r$data
      r
    }
    for(year in years) {
      progress("SCAN",paste(policy,year))
      cell_path <- file.path(out,"cells",sprintf("h3-scan-%s-%d.rds",policy,year))
      cached <- m2_cell_read(cell_path,signature)
      if(!is.null(cached)) {
        runs <- c(runs,cached$runs);cases <- c(cases,cached$cases)
        scans <- c(scans,cached$scans)
        progress("SCAN_REUSED_CHECKPOINT",paste(policy,year))
        next
      }
      previous_runs <- names(runs);previous_cases <- length(cases);previous_scans <- length(scans)
      tag <- paste(policy,year,sep="/")
      base <- add(year,paste0(tag,"/baseline"))
      for(h in c27$sensitivity_steps) {
        minus <- add(year,paste0(tag,"/Iminus/",h),z=-h)
        plus <- add(year,paste0(tag,"/Iplus/",h),z=h)
        eligible <- all(vapply(list(base,minus,plus),function(r)
          r$peak$status=="ELIGIBLE",logical(1)))
        if(eligible) {
          f <- p10_assignment(base$peak$peaks,plus$peak$peaks,p10_defaults()$identity)
          g <- p10_assignment(base$peak$peaks,minus$peak$peaks,p10_defaults()$identity)
          eligible <- any(f$status=="MATCHED" & f$left_peak==base$peak$peak &
            f$right_peak==plus$peak$peak) &&
            any(g$status=="MATCHED" & g$left_peak==base$peak$peak &
            g$right_peak==minus$peak$peak)
        }
        scans[[length(scans)+1L]] <- data.frame(policy=policy,year=year,h=h,
          eligible=eligible,reason=if(eligible) "" else "PEAK_OR_BRANCH_INELIGIBLE",
          E_I=if(eligible) (plus$peak$time-minus$peak$time)/(2*h) else NA_real_,
          baseline_week=base$data$peak_week,minus_week=minus$data$peak_week,
          plus_week=plus$data$peak_week,
          baseline_peak_id=if(is.null(base$peak[["peak"]])) NA_integer_ else base$peak[["peak"]],
          minus_peak_id=if(is.null(minus$peak[["peak"]])) NA_integer_ else minus$peak[["peak"]],
          plus_peak_id=if(is.null(plus$peak[["peak"]])) NA_integer_ else plus$peak[["peak"]],
          baseline_peak_status=base$peak$status,minus_peak_status=minus$peak$status,
          plus_peak_status=plus$peak$status,
          discrete_I=(plus$data$peak_week-minus$data$peak_week)/(2*h),
          centroid_I=(plus$data$centroid-minus$data$centroid)/(2*h))
      }
      zgrid <- seq(-log(1.2),log(1.2),length.out=17L)
      for(z in zgrid) {
        id <- paste0(tag,"/path/",sprintf("%.12f",z))
        if(is.null(runs[[id]])) add(year,id,z=z)
      }
      m2_cell_write(cell_path,signature,list(runs=runs[setdiff(names(runs),previous_runs)],
        cases=cases[seq.int(previous_cases+1L,length(cases))],
        scans=scans[seq.int(previous_scans+1L,length(scans))]))
    }
    for(recipient in c(2015L,2017L)) {
      donor <- if(recipient==2015L) 2017L else 2015L
      a <- p27_window(clock,recipient,TRUE)[1L]-1L
      b <- p27_window(clock,donor,TRUE)[1L]-1L
      r <- unlist(score$states[a,c("S","I")],use.names=TRUE)
      d <- unlist(score$states[b,c("S","I")],use.names=TRUE)
      for(axis in c("S","I","SI")) {
        mix <- c(S=if(grepl("S",axis,fixed=TRUE)) d["S"] else r["S"],
          I=if(grepl("I",axis,fixed=TRUE)) d["I"] else r["I"])
        names(mix) <- c("S","I")
        id <- paste(policy,recipient,"swap",axis,sep="/")
        tryCatch(add(recipient,id,donor_state=mix),error=function(e) {
          cases[[length(cases)+1L]] <<- data.frame(policy=policy,year=recipient,
            scenario=id,u=NA_real_,z=NA_real_,S_restart=mix["S"],
            I_restart=mix["I"],R_restart=NA_real_,saved_mass=sum(score$states[a,c("S","I","R")]),
            peak_week=NA_real_,quadratic_time=NA_real_,peak_status=paste0("INFEASIBLE: ",conditionMessage(e)),
            peak_height=NA_real_,centroid=NA_real_,burden=NA_real_)
        })
      }
    }
  }
  scan_columns <- c("policy","year","h","eligible","reason","E_I",
    "baseline_week","minus_week","plus_week","baseline_peak_id",
    "minus_peak_id","plus_peak_id","baseline_peak_status",
    "minus_peak_status","plus_peak_status","discrete_I","centroid_I")
  scans <- lapply(scans,function(q) {
    for(column in setdiff(scan_columns,names(q))) q[[column]] <- NA_integer_
    q[1L,scan_columns,drop=FALSE]
  })
  state <- do.call(rbind,cases);scan <- do.call(rbind,scans)
  p27_write(scan,out,"state_memory_scan")
  local <- do.call(rbind,lapply(years,function(year)
    do.call(rbind,lapply(c("raw","blend"),function(policy) {
      q <- scan[scan$year==year & scan$policy==policy,]
      eligible <- q[q$eligible & is.finite(q$E_I),]
      if(nrow(eligible)) eligible <- eligible[order(eligible$h),]
      chosen <- if(nrow(eligible)) eligible[(nrow(eligible)+1L)%/%2L,] else q[1L,]
      data.frame(season_id=sprintf("%d/%02d",year,(year+1L)%%100L),policy=policy,
        tested_steps=nrow(q),eligible_steps=nrow(eligible),
        selected_h=if(nrow(eligible)) chosen$h else NA_real_,
        selected_E_I=if(nrow(eligible)) chosen$E_I else NA_real_,
        selection_rule="MEDIAN_ELIGIBLE_LOG_STEP",
        step_range_I=if(nrow(eligible)>1L) diff(range(eligible$E_I)) else NA_real_,
        reason=if(nrow(eligible)) "" else paste(unique(q$reason),collapse=";"))
    }))))
  p27_write(local,out,"state_local_summary")
  # Track branches separately by policy and recipient season.
  branches <- list();intervals <- list()
  for(policy in c("raw","blend")) for(year in years) {
    progress("REFINE",paste(policy,year))
    cell_path <- file.path(out,"cells",sprintf("h3-refine-%s-%d.rds",policy,year))
    cached <- m2_cell_read(cell_path,signature)
    if(!is.null(cached)) {
      runs <- c(runs,cached$runs);cases <- c(cases,cached$cases)
      branches[[length(branches)+1L]] <- cached$branches
      intervals[[length(intervals)+1L]] <- cached$intervals
      progress("REFINE_REUSED_CHECKPOINT",paste(policy,year))
      next
    }
    previous_runs <- names(runs);previous_cases <- length(cases)
    score <- if(policy=="raw") ref$score else p26_score(ref$fit$endpoint,
      inputs,p22_problem(inputs),c26,"M2",ref$fit$penalty,blend=TRUE)
    ids <- names(runs)[grepl(paste0("^",policy,"/",year,"/path/"),names(runs))]
    round <- 0L
    repeat {
      cache <- lapply(ids,function(id) {
        r <- runs[[id]]$record;r$grid$v <- runs[[id]]$data$z
        r$grid$path <- "I_ONLY"
        if(nrow(r$peaks)) {r$peaks$v <- r$grid$v;r$peaks$path <- "I_ONLY"}
        r
      })
      names(cache) <- ids
      tracked <- p25_track(cache,phase25_config())
      gaps <- p27_height_gaps(tracked$grid,tracked$peaks,phase25_config())
      edges <- p27_classify(tracked$grid,tracked$peaks,gaps)
      trigger <- edges$left_week!=edges$right_week | edges$left_count!=edges$right_count |
        edges$classification %in% c("DOMINANCE_EXCHANGE","UNRESOLVED","SAMPLED_PEAK_LOSS","WINDOW_SELECTION")
      candidate <- (edges$v_left+edges$v_right)/2
      candidate <- candidate[trigger & edges$width>c27$refinement_width]
      round <- round+1L
      if(!length(candidate) || round>c27$refinement_rounds ||
         length(ids)>=c27$refinement_cap) break
      candidate <- head(candidate,c27$refinement_cap-length(ids))
      for(value in candidate) {
        id <- paste0(policy,"/",year,"/path/",sprintf("%.15f",value))
        if(!is.null(runs[[id]])) next
        r <- tryCatch(m2_state_case(score,inputs,clock,year,policy,id,z=value),
          error=function(e) e)
        if(inherits(r,"error")) {
          if(!grepl("Infeasible matched state",conditionMessage(r),fixed=TRUE)) stop(r)
          refinement_failures[[length(refinement_failures)+1L]] <- data.frame(
            policy=policy,year=year,z=value,reason=conditionMessage(r))
        } else {
          runs[[id]] <- r;cases[[length(cases)+1L]] <- r$data
          ids <- c(ids,id)
        }
      }
    }
    branches[[length(branches)+1L]] <- transform(tracked$grid,policy=policy)
    revised <- m2_reclassify(edges,tracked$grid,tracked$peaks,
      p27_bind(lapply(cache,`[[`,"curve")))
    intervals[[length(intervals)+1L]] <- transform(revised,policy=policy)
    m2_cell_write(cell_path,signature,list(runs=runs[setdiff(names(runs),previous_runs)],
      cases=if(length(cases)>previous_cases) cases[seq.int(previous_cases+1L,length(cases))] else list(),
      branches=branches[[length(branches)]],intervals=intervals[[length(intervals)]]))
  }
  state <- do.call(rbind,cases)
  p27_write(state,out,"state_memory")
  p27_write(if(length(refinement_failures)) do.call(rbind,refinement_failures) else
    data.frame(policy=character(),year=integer(),z=numeric(),reason=character()),
    out,"state_refinement_failures")
  p27_write(p27_bind(lapply(names(runs),function(id)
    transform(runs[[id]]$trajectory,scenario=id))),out,"state_trajectories")
  p27_write(p27_bind(branches),out,"state_branches")
  p27_write(p27_bind(intervals),out,"state_intervals")
  effects <- list()
  for(policy in c("raw","blend")) for(year in years) {
    q <- state[state$policy==policy & state$year==year &
      grepl("/path/",state$scenario,fixed=TRUE),]
    lo <- q[which.min(q$z),];hi <- q[which.max(q$z),]
    baseline <- state[state$policy==policy & state$year==year &
      state$scenario==paste(policy,year,"baseline",sep="/"),]
    valid <- nrow(baseline)==1L && all(c(lo$peak_status,baseline$peak_status,
      hi$peak_status)=="ELIGIBLE")
    difference <- function(a,b,column) if(valid) a[[column]]-b[[column]] else NA_real_
    effects[[length(effects)+1L]] <- data.frame(policy=policy,year=year,
      low_z=lo$z,high_z=hi$z,
      I_low_ratio=exp(lo$z),I_high_ratio=exp(hi$z),
      baseline_time=baseline$quadratic_time,
      low_to_baseline=difference(baseline,lo,"quadratic_time"),
      baseline_to_high=difference(hi,baseline,"quadratic_time"),
      timing_change=difference(hi,lo,"quadratic_time"),
      discrete_change=difference(hi,lo,"peak_week"),
      centroid_change=hi$centroid-lo$centroid,
      height_change=hi$peak_height-lo$peak_height,
      burden_change=hi$burden-lo$burden,
      peak_status=if(valid) "ELIGIBLE" else "INELIGIBLE_PEAK_OR_SUPPORT",
      at_least_one_week=if(valid) abs(hi$quadratic_time-lo$quadratic_time)>=1 else NA,
      at_least_two_weeks=if(valid) abs(hi$quadratic_time-lo$quadratic_time)>=2 else NA,
      at_least_four_weeks=if(valid) abs(hi$quadratic_time-lo$quadratic_time)>=4 else NA)
  }
  p27_write(do.call(rbind,effects),out,"state_effects")
  interactions <- list()
  for(policy in c("raw","blend")) for(year in c(2015L,2017L)) {
    interactions[[length(interactions)+1L]] <- m2_interaction_row(state,policy,year)
  }
  p27_write(do.call(rbind,interactions),out,"state_interactions")
  list(state=state,scan=scan)
}
m2_report <- function(evidence,out,state_dir=NULL) {
  h1 <- utils::read.csv(file.path(evidence,"h1_event_summary.csv"))
  edges <- utils::read.csv(file.path(evidence,"classification_audit.csv"))
  h2 <- utils::read.csv(file.path(evidence,"h2_magnitude_summary.csv"))
  inv <- utils::read.csv(file.path(evidence,"existing_experiment_inventory.csv"))
  state <- if(!is.null(state_dir)) utils::read.csv(file.path(state_dir,"state_memory.csv")) else NULL
  effects <- if(!is.null(state_dir)) utils::read.csv(file.path(state_dir,"state_effects.csv")) else NULL
  interactions <- if(!is.null(state_dir))
    utils::read.csv(file.path(state_dir,"state_interactions.csv")) else NULL
  local_h3 <- if(!is.null(state_dir) && file.exists(file.path(state_dir,
    "state_local_summary.csv"))) utils::read.csv(file.path(state_dir,
    "state_local_summary.csv")) else NULL
  full_h1 <- file.exists(file.path(evidence,"h1_season_summary.csv"))
  season_h1 <- if(full_h1) utils::read.csv(file.path(evidence,"h1_season_summary.csv")) else NULL
  movement_h1 <- if(file.exists(file.path(evidence,"h1_movement_by_mechanism.csv")))
    utils::read.csv(file.path(evidence,"h1_movement_by_mechanism.csv")) else NULL
  season_h2 <- if(file.exists(file.path(evidence,"h2_season_summary.csv")))
    utils::read.csv(file.path(evidence,"h2_season_summary.csv")) else NULL
  full_h3 <- !is.null(effects) && length(unique(effects$year))==10L &&
    all(2013:2022 %in% effects$year)
  coverage <- do.call(rbind,lapply(2013:2022,function(y)
    do.call(rbind,lapply(c("raw","blend"),function(policy) {
      qh1 <- if(full_h1) season_h1[season_h1$season_id==
        sprintf("%d/%02d",y,(y+1L)%%100L) & season_h1$policy==policy,] else data.frame()
      qh2 <- h2[h2$reference=="M2" & h2$design=="MATCHED_SEASON_ORIGIN" &
        h2$year==y & h2$policy==policy,]
      qh3 <- if(!is.null(effects)) effects[effects$year==y &
        effects$policy==policy,] else data.frame()
      data.frame(season_id=sprintf("%d/%02d",y,(y+1L)%%100L),policy=policy,
        h1_matched_beta_only=if(nrow(qh1)) "COMPLETED" else "IMPLEMENTED_NOT_RUN",
        h2_matched_S_beta=if(nrow(qh2)==5L)
          if(any(qh2$paired_eligible)) "COMPLETED" else "COMPLETED_NUMERICALLY_INELIGIBLE"
          else "IMPLEMENTED_NOT_RUN",
        h3_local_I=if(nrow(qh3)) "COMPLETED" else "IMPLEMENTED_NOT_RUN",
        h3_finite_I_path=if(nrow(qh3)) if(!"peak_status" %in% names(qh3) ||
          qh3$peak_status=="ELIGIBLE") "COMPLETED" else
          "COMPLETED_NUMERICALLY_INELIGIBLE" else "IMPLEMENTED_NOT_RUN")
    }))))
  p27_write(coverage,out,"coverage_matrix")
  group <- unique(h2[c("reference","policy","design")])
  counts <- lapply(seq_len(nrow(group)),function(i) {
    k <- group[i,];q <- h2[h2$reference==k$reference & h2$policy==k$policy &
      h2$design==k$design,]
    data.frame(k,attempted=nrow(q),eligible=sum(q$paired_eligible),
      beta_larger=sum(q$magnitude_beta_gt_S),negative_S=sum(q$E_S<0,na.rm=TRUE),
      negative_beta=sum(q$E_beta<0,na.rm=TRUE))
  })
  counts <- do.call(rbind,counts)
  p27_write(counts,out,"h2_counts")
  exceptions <- h2[h2$reference=="M2" & h2$design=="MATCHED_SEASON_ORIGIN" &
    h2$paired_eligible & !h2$magnitude_beta_gt_S,
    c("policy","design","year","h","E_S","E_beta")]
  p27_write(exceptions,out,"h2_magnitude_exceptions")
  groups <- unique(h1[c("policy","path")])
  bounds <- do.call(rbind,lapply(seq_len(nrow(groups)),function(i) {
    k <- groups[i,];q <- h1[h1$policy==k$policy & h1$path==k$path & h1$timing_event,]
    data.frame(k,events=nrow(q),
      resolved_displacement=sum(q$displacement & !q$unresolved & !q$mixed),
      possible_displacement=sum(q$displacement | q$unresolved | q$mixed),
      unresolved=sum(q$unresolved),mixed=sum(q$mixed))
  }))
  p27_write(bounds,out,"h1_bounds")
  exchanges <- edges[edges$classification=="DOMINANCE_EXCHANGE",
    c("policy","season","path","origin","actual_window","v_left","v_right",
      "left_week","right_week","early_minus_late_left","early_minus_late_right",
      "left_branch","right_branch")]
  p27_write(exchanges,out,"h1_exchanges")
  claims <- data.frame(claim=c("H1a","H1b","H2_magnitude","H2_direction","H3"),
    verdict=c(if(any(edges$classification=="DOMINANCE_EXCHANGE") &&
      any(edges$classification=="CONTINUOUS_DISPLACEMENT")) "SUPPORTED_UNDER_M2" else "INCONCLUSIVE_UNDER_M2",
      if(full_h1 && any(season_h1$displacement_event_fraction>.5,na.rm=TRUE))
        "SEASON_DEPENDENT" else "INCONCLUSIVE_UNDER_M2",
      if(any(counts$beta_larger>0)) "CONTEXT_DEPENDENT" else "INCONCLUSIVE_UNDER_M2",
      "CONTEXT_DEPENDENT",
      if(is.null(effects)) "NOT_RUN" else if(!full_h3) "PARTIAL_MODEL_EVIDENCE" else
        if(any(is.finite(effects$timing_change) & abs(effects$timing_change)>1e-6))
          "SUPPORTED_UNDER_M2" else "NOT_SUPPORTED_UNDER_M2"),
    endpoint=inv$parameter_hash[1L],
    source_code_md5=inv$archived_phase27_md5[1L],
    classification_code_md5=unname(tools::md5sum("R/m2_hypotheses.R")),
    current_phase27_md5=unname(tools::md5sum("R/phase27.R")),
    report_code_md5=unname(tools::md5sum("scripts/m2_hypotheses.R")),
    evidence_scope=c(rep(if(full_h1) "MIXED_HISTORICAL_AND_NEW_SIMULATION" else
      inv$evidence_scope[1L],2L),rep(inv$evidence_scope[1L],2L),
      if(is.null(effects)) "NOT_RUN" else if(file.exists(file.path(state_dir,
        "run_provenance.csv"))) "NEW_SIMULATION" else "HISTORICAL_EXPERIMENT"),
    intervention_origin=c("2013 and matched-season","2013 and matched-season",
      "see H2 row-level design","see H2 row-level design","matched-season"),
    actual_window=c("epidemic season through 2019; later calendar years on full-history path",
      "epidemic season through 2019; later calendar years on full-history path",
      "see H2 row-level design","see H2 row-level design","epidemic season"),
    estimand=c("selected reconstructed peak week","parent-grid timing event",
      "absolute finite log-step slope","signed finite log-step slope",
      "matched-origin reconstructed peak time"),stringsAsFactors=FALSE)
  p27_write(claims,out,"claims")
  lines <- c("# M2 hypothesis evidence", "",
    paste("Endpoint parameter hash:",inv$parameter_hash[1L]),
    paste("Archived Phase 27 code MD5:",inv$archived_phase27_md5[1L]),
    paste("Source scope:",inv$evidence_scope[1L]),"",
    "## Coverage",paste(capture.output(print(coverage,row.names=FALSE)),collapse="\n"),"",
    "## H1",sprintf("Saved refined intervals: %d dominance exchanges, %d displacements.",
      sum(edges$classification=="DOMINANCE_EXCHANGE"),
      sum(edges$classification=="CONTINUOUS_DISPLACEMENT")),
    "Original-grid timing-event bounds by policy and path:",
    paste(capture.output(print(bounds,row.names=FALSE)),collapse="\n"),
    "Targeted signed competitor-height reversals:",
    paste(capture.output(print(exchanges,row.names=FALSE)),collapse="\n"),
    "The bounds keep mixed and unresolved parent intervals. Blending is a numerical robustness arm.",
    if(full_h1) paste(capture.output(print(season_h1,row.names=FALSE)),collapse="\n") else
      "Full-season H1 experiment has not run; no full-period H1b verdict is available.",
    if(!is.null(movement_h1)) paste(capture.output(print(movement_h1,row.names=FALSE)),collapse="\n"),
    "Absolute movement may exceed signed net endpoint movement. Equal-season summaries retain no-shift seasons.","",
    "## H2",paste(capture.output(print(counts,row.names=FALSE)),collapse="\n"),
    if(!is.null(season_h2)) paste(capture.output(print(season_h2,row.names=FALSE)),collapse="\n"),
    "Eligible matched-origin magnitude-order exceptions:",
    paste(capture.output(print(exceptions,row.names=FALSE)),collapse="\n"),
    "Magnitude and sign are separate conclusions. Rows are year-step diagnostics, not independent samples.",
    "Matched-origin 2021/22 contains the saved magnitude-order exceptions; 2013-origin 2017 slopes can be positive.","",
    "## H3",if(is.null(state)) "Not run. H3 is not inferred from immediate incidence alone." else
      paste("State-memory scenarios supplied:",nrow(state)),
    if(!is.null(local_h3)) paste(capture.output(print(local_h3,row.names=FALSE)),collapse="\n"),
    if(!is.null(effects)) paste(capture.output(print(effects,row.names=FALSE)),collapse="\n"),
    "Finite donor-state factorial contrasts:",
    if(!is.null(interactions)) paste(capture.output(print(interactions,row.names=FALSE)),collapse="\n") else
      "Not run.",
    if(!full_h3) "Full-season H3 coverage is incomplete; no full-period H3 verdict is available." else "",
    "Initial-I timing effects are conditional fitted-model responses, not observed-history causal attribution.","",
    "The selected M2 endpoint is retained. Prior Phase 26 adequacy failures remain in docs/phase26_interpretation.md.")
  writeLines(lines,file.path(out,"evidence_report.md"))
  grDevices::pdf(file.path(out,"h2_slopes.pdf"))
  q <- h2[h2$reference=="M2" & h2$paired_eligible,]
  if(nrow(q)) {
    plot(q$E_S,q$E_beta,col=ifelse(q$design=="2013_ORIGIN",1,2),
      xlab="E_S",ylab="E_beta",main="M2 finite log-step responses")
    abline(0,1,lty=2)
  } else plot.new()
  grDevices::dev.off()
  grDevices::pdf(file.path(out,"chronological_seasons.pdf"),width=9,height=10)
  old_par <- graphics::par(mfrow=c(3,1),mar=c(4,5,2,1))
  if(!is.null(season_h1)) {
    q <- season_h1[season_h1$policy=="raw",]
    graphics::plot(seq_len(nrow(q)),q$total_absolute_movement,type="b",pch=16,
      xaxt="n",xlab="Season",ylab="Absolute peak movement (weeks)",main="H1 raw")
    graphics::axis(1,at=seq_len(nrow(q)),labels=q$season_id,las=2,cex.axis=.7)
  } else graphics::plot.new()
  if(!is.null(season_h2)) {
    q <- season_h2[season_h2$policy=="raw",]
    graphics::plot(seq_len(nrow(q)),q$selected_E_beta,type="b",pch=16,
      xaxt="n",xlab="Season",ylab="Peak-time slope",main="H2 raw")
    graphics::lines(seq_len(nrow(q)),q$selected_E_S,type="b",pch=1,lty=2)
    graphics::axis(1,at=seq_len(nrow(q)),labels=q$season_id,las=2,cex.axis=.7)
    graphics::legend("topright",c("beta","S"),pch=c(16,1),lty=c(1,2),bty="n")
  } else graphics::plot.new()
  if(!is.null(effects)) {
    q <- effects[effects$policy=="raw",]
    q <- q[order(q$year),]
    graphics::plot(q$year,q$timing_change,type="b",pch=16,xaxt="n",
      xlab="Season",ylab="Low-to-high peak change (weeks)",main="H3 raw")
    graphics::axis(1,at=q$year,labels=sprintf("%d/%02d",q$year,(q$year+1L)%%100L),
      las=2,cex.axis=.7)
    graphics::abline(h=0,lty=2)
  } else graphics::plot.new()
  graphics::par(old_par)
  grDevices::dev.off()
  invisible(claims)
}
