# Targeted Phase 16 validation uses saved coarse evidence and fixed candidate paths.
p16t_spec <- function() data.frame(
  target_id=c("primary","secondary_low","secondary_high"),
  priority=c(1L,2L,2L),profile_id=c("profile-004","profile-009","profile-009"),
  year=4L,axis="beta",
  left_scenario=c("annual-0019","annual-0044","annual-0041"),
  right_scenario=c("annual-0016","annual-0041","annual-0045"),
  stringsAsFactors=FALSE)
p16t_budget <- function() list(refinement_depth=2L,refinement_per_stage=3L,
  neighborhood_scenarios=6L,primary_only_neighborhood=TRUE)
p16t_check_source <- function(coarse,old,c) {
  if(!identical(coarse$settings,c)||!identical(old$settings,c)||
    !identical(coarse$source_hash,old$source_hash)) stop("Coarse and validation provenance mismatch")
  if(!identical(coarse$type,"phase16_evaluation")||!identical(old$type,"phase16_validation"))
    stop("Expected Phase 16 coarse and validation artifacts")
  spec<-p16t_spec()
  for(i in seq_len(nrow(spec))) {
    z<-spec[i,,drop=FALSE]
    p<-coarse$profiles[coarse$profiles$profile_id==z$profile_id,,drop=FALSE]
    if(nrow(p)!=1L||p$amplitude<=0||p$continuous_maxima!=1L||
      p$dense_grid_maxima!=1L||p$sampled_status!="UNIQUE_SAMPLED_MAXIMUM")
      stop("Target must have verified positive-amplitude annual-only forcing")
    q<-coarse$scenarios[coarse$scenarios$scenario_id %in% c(z$left_scenario,z$right_scenario),,drop=FALSE]
    if(nrow(q)!=2L||any(q$profile_id!=z$profile_id)||any(q$status!="VALID")||
      any(q$u!=0)||length(unique(q$v))!=2L) stop("Saved target endpoints unavailable")
    hit<-old$screening[old$screening$year==z$year&old$screening$profile_id==z$profile_id&
      old$screening$axis=="beta"&old$screening$left_scenario==z$left_scenario&
      old$screening$right_scenario==z$right_scenario&
      old$screening$comparison_span=="ADJACENT",,drop=FALSE]
    if(nrow(hit)!=1L) stop("Target not present as a saved adjacent comparison")
  }
  invisible(spec)
}
p16t_evaluation_catalog <- function(root,source_hash,c) {
  dirs<-list.dirs(root,full.names=TRUE,recursive=FALSE)
  rows<-list();objects<-list()
  for(dir in dirs) {
    path<-file.path(dir,"result.rds")
    if(!file.exists(path)||!file.exists(file.path(dir,"STATUS"))) next
    x<-tryCatch(read_phase2_result(path)$result,error=function(e) NULL)
    if(is.null(x)||!identical(x$type,"phase16_evaluation")||
      !identical(x$source_hash,source_hash)||!identical(x$settings,c)) next
    objects[[normalizePath(path)]]<-x
    rows[[length(rows)+1L]]<-data.frame(path=normalizePath(path),
      md5=unname(tools::md5sum(path)),scenario_count=nrow(x$scenarios))
  }
  list(manifest=if(length(rows)) do.call(rbind,rows) else data.frame(),objects=objects)
}
p16t_input_status <- function(plan_status,found,consumed) {
  if(plan_status!="PLANNED") return(plan_status)
  if(!length(found)) return("NO_SAVED_EVALUATION_FOUND")
  if(!any(found %in% consumed)) return("SAVED_EVALUATION_OMITTED_FROM_VALIDATION")
  "EVALUATION_CONSUMED"
}
p16t_old_input_audit <- function(old_dir,coarse,c,artifact_root) {
  old_path<-file.path(normalizePath(old_dir,mustWork=TRUE),"result.rds")
  old_saved<-read_phase2_result(old_path)
  consumed<-normalizePath(old_saved$metadata$input_paths,mustWork=FALSE)
  catalog<-p16t_evaluation_catalog(artifact_root,coarse$source_hash,c)
  rows<-list()
  for(kind in c("refinement","neighborhood")) {
    path<-file.path(old_dir,paste0(kind,"-plan"),"result.rds")
    if(!file.exists(path)) {
      rows[[length(rows)+1L]]<-data.frame(plan_kind=kind,scenario_id="",
        plan_status="PLAN_ARTIFACT_ABSENT",evaluation_path="",evaluation_md5="",
        input_status="PLAN_ARTIFACT_ABSENT")
      next
    }
    plan<-read_phase2_result(path)$result
    for(i in seq_len(nrow(plan$scenarios))) {
      q<-plan$scenarios[i,,drop=FALSE]
      found<-character()
      for(p in names(catalog$objects)) {
        e<-catalog$objects[[p]]$scenarios
        a<-e[e$scenario_id==q$scenario_id& e$profile_id==q$profile_id&
          abs(e$u-q$u)<1e-12&abs(e$v-q$v)<1e-12,,drop=FALSE]
        if(nrow(a)==1L) found<-c(found,p)
      }
      status<-p16t_input_status(q$status,found,consumed)
      rows[[length(rows)+1L]]<-data.frame(plan_kind=kind,scenario_id=q$scenario_id,
        plan_status=q$status,evaluation_path=paste(found,collapse=";"),
        evaluation_md5=paste(vapply(found,function(p) unname(tools::md5sum(p)),character(1)),collapse=";"),
        input_status=status)
    }
  }
  do.call(rbind,rows)
}
p16t_plan_row <- function(r,z,level,left,right) {
  mid<-(left$v+right$v)/2
  data.frame(scenario_id=sprintf("target-%s-l%d",z$target_id,level),
    profile_id=z$profile_id,axis="beta",u=0,v=mid,s=r$accepted$S0,
    a=exp(mid),level=level,status="PLANNED",reason="TARGETED_PATH_MIDPOINT",
    target_id=z$target_id,left_anchor=left$scenario_id,right_anchor=right$scenario_id,
    stringsAsFactors=FALSE)
}
p16t_plan_object <- function(r,rows,profiles,c,kind) {
  list(type="phase16_plan",status=paste0("TARGETED_",toupper(kind),"_PLAN"),settings=c,
    source_hash=r$source_hash,accepted=r$accepted,profiles=profiles,
    scenarios=rows,search_scope=paste0("TARGETED_",toupper(kind),"_ONLY"))
}
p16t_initial_refinement <- function(coarse,c) {
  spec<-p16t_spec();rows<-list()
  for(i in seq_len(nrow(spec))) {
    z<-spec[i,,drop=FALSE]
    l<-coarse$scenarios[coarse$scenarios$scenario_id==z$left_scenario,,drop=FALSE]
    r<-coarse$scenarios[coarse$scenarios$scenario_id==z$right_scenario,,drop=FALSE]
    rows[[i]]<-p16t_plan_row(coarse,z,1L,l,r)
  }
  p16t_plan_object(coarse,do.call(rbind,rows),coarse$profiles,c,"refinement")
}
p16t_initial_neighborhood <- function(coarse,old,c) {
  spec<-p16t_spec()[1,,drop=FALSE]
  screen<-old$screening[old$screening$year==spec$year&
    old$screening$left_scenario==spec$left_scenario&
    old$screening$right_scenario==spec$right_scenario&
    old$screening$comparison_span=="ADJACENT",,drop=FALSE]
  settings<-c;settings$neighborhood$cap<-p16t_budget()$neighborhood_scenarios
  plan<-p16_neighborhood_plan(coarse,screen,settings)
  old_ids<-plan$profiles$profile_id
  plan$profiles$profile_id<-sub("^neighbor-","target-neighbor-",old_ids)
  plan$scenarios$profile_id<-sub("^neighbor-","target-neighbor-",plan$scenarios$profile_id)
  plan$scenarios$scenario_id<-sub("^neighbor-","target-neighbor-",plan$scenarios$scenario_id)
  plan$settings<-c;plan$status<-"TARGETED_PRIMARY_NEIGHBORHOOD_PLAN"
  plan$search_scope<-"PRIMARY_PROFILE_ONLY_PRESPECIFIED_NEIGHBORHOOD"
  plan
}
p16t_stage2 <- function(r,plan1,c,d) {
  rows<-list();spec<-p16t_spec()
  for(i in seq_len(nrow(spec))) {
    z<-spec[i,,drop=FALSE];midid<-plan1$scenarios$scenario_id[plan1$scenarios$target_id==z$target_id]
    if(length(midid)!=1L) stop("Target stage-one request missing")
    mid<-r$scenarios[r$scenarios$scenario_id==midid,,drop=FALSE]
    l<-r$scenarios[r$scenarios$scenario_id==z$left_scenario,,drop=FALSE]
    rr<-r$scenarios[r$scenarios$scenario_id==z$right_scenario,,drop=FALSE]
    if(nrow(mid)!=1L||mid$status!="VALID") {
      q<-p16t_plan_row(r,z,2L,l,rr);q$status<-"DEFERRED";q$reason<-"STAGE_ONE_INPUT_MISSING"
    } else {
      l$axis<-mid$axis<-rr$axis<-"beta"
      a<-p16_pair(r,l,mid,z$year,c,d);b<-p16_pair(r,mid,rr,z$year,c,d)
      choose_left<-is.finite(a$absolute_weekly_difference)&&
        (!is.finite(b$absolute_weekly_difference)||a$absolute_weekly_difference>=b$absolute_weekly_difference)
      if(choose_left) q<-p16t_plan_row(r,z,2L,l,mid) else q<-p16t_plan_row(r,z,2L,mid,rr)
      if(!is.finite(a$absolute_weekly_difference)&&!is.finite(b$absolute_weekly_difference)) {
        q$status<-"DEFERRED";q$reason<-"NO_FINITE_STAGE_ONE_RESPONSE"
      }
    }
    rows[[i]]<-q
  }
  p16t_plan_object(r,do.call(rbind,rows),r$profiles,c,"refinement")
}
p16t_coverage <- function(plan,r,stage) {
  x<-p16_coverage(plan,r)
  if(nrow(x)) x$stage<-stage
  x
}
p16t_endpoint_screen <- function(r,c,d) {
  rows<-list();spec<-p16t_spec()
  for(i in seq_len(nrow(spec))) {
    z<-spec[i,,drop=FALSE]
    l<-r$scenarios[r$scenarios$scenario_id==z$left_scenario,,drop=FALSE]
    rr<-r$scenarios[r$scenarios$scenario_id==z$right_scenario,,drop=FALSE]
    l$axis<-rr$axis<-"beta"
    q<-p16_pair(r,l,rr,z$year,c,d)
    q$comparison_span<-"TARGET_ENDPOINT"
    rows[[i]]<-q
  }
  do.call(rbind,rows)
}
p16t_paths <- function(value) if(is.null(value)) character() else
  trimws(strsplit(value,",",fixed=TRUE)[[1]])
p16t_load_many <- function(value) lapply(p16t_paths(value),p16_load,type="phase16_evaluation")
p16t_serialize_plan <- function(plan,out,root,config_path,inputs,name) {
  path<-file.path(out,name);if(!dir.create(path)) stop("Cannot create targeted child plan")
  meta<-phase2_metadata(root,config_path,inputs,paste0("phase16-targeted-",name))
  meta$phase<-16L;meta$method<-"Targeted saved-data plan; no simulation or fitting"
  write_phase2_result(plan,path,meta)
  p16_csv(plan$scenarios,file.path(path,"scenarios.csv"))
  invisible(path)
}
p16t_claims <- function(checks,coverage,spec,c) {
  complete<-nrow(coverage)>0L&&all(coverage$observed_status=="VALID"&coverage$coordinate_matches)
  rows<-list()
  for(i in seq_len(nrow(spec))) {
    z<-spec[i,,drop=FALSE]
    lo<-if(z$target_id=="secondary_high") 0 else -log(1.2)
    hi<-if(z$target_id=="secondary_high") log(1.2) else 0
    x<-checks[checks$profile_id==z$profile_id&checks$year==z$year&checks$axis=="beta"&
      checks$comparison_span=="ADJACENT"&checks$left_coordinate>=lo-1e-8&
      checks$right_coordinate<=hi+1e-8,,drop=FALSE]
    endpoint<-checks[checks$profile_id==z$profile_id&checks$year==z$year&
      checks$comparison_span=="TARGET_ENDPOINT"&
      checks$left_scenario==z$left_scenario&checks$right_scenario==z$right_scenario,,drop=FALSE]
    coherent<-nrow(x)>0L&&all(x$same_branch&x$forward_reverse_consistent&
      x$state_valid&!x$low_incidence&!x$tie_or_plateau)
    finite<-nrow(endpoint)==1L&&endpoint$validated_status=="SUPPORTED_FINITE_DISPLACEMENT"&&coherent
    abrupt<-any(x$abruptness_supported&x$absolute_weekly_difference>=c$large_weeks)
    switch<-any(x$validated_status=="SUPPORTED_GENUINE_SWITCH_FINITE_HORIZON")
    window<-any(c(x$validated_status,endpoint$validated_status)=="SUPPORTED_WINDOW_SELECTION_EFFECT")
    rows[[i]]<-data.frame(target_id=z$target_id,profile_id=z$profile_id,year=z$year,
      input_status=if(complete) "COMPLETE" else "PARTIAL_EVIDENCE",
      finite_displacement=if(complete) finite else NA,
      abrupt_finite_response=if(complete) abrupt else NA,
      genuine_dominance_switch=if(complete) switch else NA,
      window_selection=if(complete) window else NA,
      unresolved=if(complete) !any(c(finite,abrupt,switch,window)) else TRUE,
      long_run="PHASE17_PENDING",comparison_rows=nrow(x))
  }
  do.call(rbind,rows)
}
p16t_branch_heights <- function(r,screen) {
  ids<-unique(c(screen$left_scenario,screen$right_scenario))
  ids<-ids[!is.na(ids)]
  rows<-list()
  for(id in ids) {
    w<-r$windows[r$windows$scenario==id&r$windows$year==4L,,drop=FALSE]
    p<-r$peaks[r$peaks$scenario==id&r$peaks$year==4L,,drop=FALSE]
    if(nrow(w)!=1L) next
    if(!nrow(p)) {
      rows[[length(rows)+1L]]<-data.frame(scenario=id,year=4L,peak_id=NA_integer_,
        week=NA_real_,peak_height=NA_real_,genuine=FALSE,retained=FALSE,
        annual_winner_position=w$position,annual_winner_height=w$height,
        annual_winner_status=w$status,is_selected=FALSE,
        signed_winner_minus_peak=NA_real_)
    } else for(i in seq_len(nrow(p))) {
      a<-p[i,,drop=FALSE]
      rows[[length(rows)+1L]]<-data.frame(scenario=id,year=4L,peak_id=a$peak_id,
        week=a$week,peak_height=a$height,genuine=a$genuine,retained=a$retained,
        annual_winner_position=w$position,annual_winner_height=w$height,
        annual_winner_status=w$status,is_selected=a$position==w$position,
        signed_winner_minus_peak=w$height-a$height)
    }
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
