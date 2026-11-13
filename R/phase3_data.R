# Temporal contracts and fixed-value provenance; these functions never optimize.
p3_text <- function(x,label) {
  if(!is.character(x) || length(x)!=1L || is.na(x) || !nzchar(x))
    stop("Set ",label,call.=FALSE)
  invisible(x)
}
p3_markers <- function(data,ids=seq_along(data$cases)) {
  data$time_markers[ids,setdiff(names(data$time_markers),"cases"),drop=FALSE]
}
p3_prefix <- function(data,end) {
  integer_scalar(end,"prefix end",2L)
  if(end>length(data$cases)) stop("Prefix exceeds data",call.=FALSE)
  x <- data; ids <- seq_len(end)
  for(n in c("cases","seasonal_week","npieff","processed_week_marker","likelihood_include")) x[[n]] <- data[[n]][ids]
  x$time_markers <- data$time_markers[ids,,drop=FALSE]
  x$observation_order <- ids
  x$provenance$training_indices <- ids
  validate_phase2_data(x)
  x
}
p3_splits <- function(data,spec,complete_years=TRUE) {
  if(is.null(spec)) return(NULL)
  if(!identical(names(spec),c("A","B","C"))) stop("Splits must be ordered A, B, C",call.=FALSE)
  ids <- lapply(names(spec),function(n) {
    b <- spec[[n]]; finite_numeric(b,paste("split",n),2L)
    if(any(b!=floor(b)) || b[1]<1 || b[2]<b[1]) stop("Empty or invalid split ",n,call.=FALSE)
    seq.int(b[1],b[2])
  });names(ids) <- names(spec)
  if(!identical(as.integer(unlist(ids,use.names=FALSE)),seq_len(length(data$cases))))
    stop("Splits must partition all analysis rows without gaps, overlap or reordering",call.=FALSE)
  for(n in names(ids)) {
    ii <- ids[[n]]; w <- data$seasonal_week[ii]
    if(length(ii)<52L) stop("Each split requires at least one 52-week cycle",call.=FALSE)
    if(isTRUE(complete_years) && (w[1]!=1L || tail(w,1)!=52L || length(ii)%%52L!=0L))
      stop("Complete-year splits must start at week 1 and end at week 52",call.=FALSE)
  }
  list(indices=ids,markers=lapply(ids,function(i) p3_markers(data,i)),
       mask=lapply(ids,function(i) data$likelihood_include[i]),
       forecast_mask=list(B=rep(TRUE,length(ids$B)),C=rep(TRUE,length(ids$C))))
}
p3_provenance <- function(cfg,training_end,forecast=TRUE) {
  fixed <- c("omega","delta","npieff","pop","mu","gamma",
             names(cfg$estimation$free)[!unlist(cfg$estimation$free)])
  rows <- lapply(fixed,function(n) {
    s <- cfg$fixed_provenance[[n]]
    if(!is.list(s)) stop("Missing fixed provenance for ",n,call.=FALSE)
    p3_text(s$citation,paste(n,"citation"));p3_text(s$category,paste(n,"category"))
    if(!s$category %in% c("external_independent","training_only","full_data","user_conditional"))
      stop("Unknown fixed-value source category: ",n,call.=FALSE)
    integer_scalar(s$last_used_index,paste(n,"last_used_index"),0L)
    if(s$category=="external_independent" && s$last_used_index!=0L)
      stop("Independent assumptions require last_used_index=0",call.=FALSE)
    if(forecast && s$category=="training_only" && s$last_used_index>training_end)
      stop("Training-only source reaches beyond training boundary: ",n,call.=FALSE)
    if(forecast && is.null(s$known_at_index)) stop("Set fixed provenance known_at_index for ",n,call.=FALSE)
    known <- if(is.null(s$known_at_index)) s$last_used_index else s$known_at_index
    integer_scalar(known,paste(n,"known_at_index"),0L)
    leakage <- s$category %in% c("full_data","user_conditional") || known>training_end || s$last_used_index>training_end
    if(forecast && n=="npieff" && s$category=="full_data")
      stop("Full-period posterior NPI cannot be used in future evaluation; use descriptive only",call.=FALSE)
    if(forecast && cfg$analysis_type=="strict_future" && leakage)
      stop("Strict future evaluation blocks full-data/conditional fixed values: ",n,call.=FALSE)
    data.frame(parameter=n,category=s$category,citation=s$citation,
      last_used_index=s$last_used_index,known_at_index=known,conditional_or_leakage=leakage)
  })
  do.call(rbind,rows)
}
p3_future_input <- function(data,ids,cfg,stage,origin) {
  f <- cfg$future[[stage]]
  if(!is.list(f)) stop("Missing future input contract for ",stage,call.=FALSE)
  finite_numeric(f$values,paste(stage,"future npieff"),length(ids))
  if(any(f$values<0)) stop("Future NPI must be nonnegative",call.=FALSE)
  if(!identical(as.numeric(f$values),as.numeric(data$npieff[ids])))
    stop("Future input must match the explicitly adapted NPI scenario",call.=FALSE)
  p3_text(f$citation,paste(stage,"future input citation"))
  if(!is.character(f$kind) || length(f$kind)!=1L || !f$kind %in% c("known","scenario","origin_generated"))
    stop("Future NPI requires known/scenario; full-data posterior inputs are forbidden",call.=FALSE)
  if(f$kind %in% c("known","origin_generated")) {
    if(f$kind=="origin_generated") p3_text(f$mechanism,paste(stage,"origin-generated mechanism"))
    integer_scalar(f$known_at_index,"future input knowledge index",0L)
    if(f$known_at_index>origin && cfg$analysis_type=="strict_future") stop("Future input not known at forecast origin",call.=FALSE)
    f$conditional <- f$known_at_index>origin
    integer_scalar(f$last_used_index,"future input last_used_index",0L)
    if(f$last_used_index>origin) {
      if(cfg$analysis_type=="strict_future") stop("Future input derived beyond forecast origin",call.=FALSE)
      f$conditional <- TRUE
    }
  } else {
    if(cfg$analysis_type!="conditional_future")
      stop("Unknown future input needs conditional_future with explicit scenario",call.=FALSE)
    p3_text(f$mechanism,paste(stage,"scenario mechanism"))
    integer_scalar(f$known_at_index,"scenario declared index",0L)
    integer_scalar(f$last_used_index,"scenario last_used_index",0L)
    # Retrospective scenarios are kept conditional, never promoted to strict.
    f$conditional <- TRUE
  }
  # There are deliberately no case counts in this interface.
  list(indices=ids,weeks=data$seasonal_week[ids],npieff=f$values,
    markers=p3_markers(data,ids),source=f)
}
p3_load_data <- function(cfg) {
  d <- cfg$data
  for(n in c("cases_file","population_file","region","aggregation","preprocessing_id")) p3_text(d[[n]],paste("data",n))
  if(!is.logical(d$include_npi) || length(d$include_npi)!=1L || is.na(d$include_npi)) stop("Set include_npi",call.=FALSE)
  integer_scalar(d$analysis_end,"analysis_end",2L)
  cases <- utils::read.csv(d$cases_file,stringsAsFactors=FALSE)
  if(!all(c("island","year","week","cases") %in% names(cases)) || anyNA(cases))
    stop("Prepared cases require island/year/week/cases without missing values",call.=FALSE)
  # Crop by explicit row index within the selected island, before adapting NPI length.
  sub <- cases[cases$island==d$region,,drop=FALSE]
  if(nrow(sub)<d$analysis_end) stop("Analysis end exceeds available rows",call.=FALSE)
  x <- adapt_preprocessed_island(sub[seq_len(d$analysis_end),,drop=FALSE],
    utils::read.csv(d$population_file,stringsAsFactors=FALSE),d$region,d$post_npi)
  if(!d$include_npi && length(x$cases)>x$definition$Nnonpi)
    stop("NPI-excluded analysis must end before the original NPI boundary",call.=FALSE)
  x$provenance$preprocessing_id <- d$preprocessing_id
  x$time_markers$original_island_index <- seq_len(length(x$cases))
  x$provenance$input_hashes <- tools::md5sum(c(d$cases_file,d$population_file,d$raw_input_paths))
  x$provenance$case_unit <- "reported_cases_per_week"
  if(!is.null(d$case_unit) && !identical(d$case_unit,x$provenance$case_unit)) stop("Unsupported case units",call.=FALSE)
  if(!is.null(d$preprocessing_manifest)) {
    manifest <- readRDS(d$preprocessing_manifest)
    if(!identical(manifest$preprocessing_id,d$preprocessing_id)) stop("Preprocessing identity mismatch",call.=FALSE)
    if(is.null(manifest$output_md5) || is.null(manifest$mapping_file) || is.null(manifest$mapping_md5))
      stop("Preprocessing manifest lacks output hashes/mapping; regenerate into a new directory",call.=FALSE)
    for(path in c(d$cases_file,d$population_file)) {
      stored <- manifest$output_md5[basename(names(manifest$output_md5))==basename(path)]
      if(length(stored)!=1L || !identical(unname(stored),unname(tools::md5sum(path))))
        stop("Prepared file hash differs from preprocessing manifest: ",path,call.=FALSE)
    }
    # Check the declared current upstream inputs as well as the prepared outputs.
    for(path in d$raw_input_paths) if(basename(path) %in% basename(names(manifest$input_md5))) {
      stored <- manifest$input_md5[basename(names(manifest$input_md5))==basename(path)]
      if(length(stored)!=1L || !identical(unname(stored),unname(tools::md5sum(path))))
        stop("Upstream input differs from preprocessing manifest: ",path,call.=FALSE)
    }
    if(!identical(manifest$mapping_file,"prepared_row_mapping.csv")) stop("Unsupported mapping filename",call.=FALSE)
    mapping_path <- file.path(dirname(d$preprocessing_manifest),manifest$mapping_file)
    if(!identical(manifest$mapping_md5,unname(tools::md5sum(mapping_path)))) stop("Raw mapping hash mismatch",call.=FALSE)
    mapping <- utils::read.csv(mapping_path,stringsAsFactors=FALSE)
    required <- c("raw_row","prefecture","year","raw_week","prepared_week","island","original_island_index")
    if(!all(required %in% names(mapping)) || anyNA(mapping[required])) stop("Invalid raw/prepared row mapping",call.=FALSE)
    if(anyDuplicated(mapping$raw_row) || any(!mapping$raw_week %in% 1:53) ||
       any(mapping$prepared_week!=pmin(mapping$raw_week,52)) ||
       any(mapping$original_island_index!=(mapping$year-2013L)*52L+mapping$prepared_week))
      stop("Inconsistent raw/prepared week mapping",call.=FALSE)
    mapping <- mapping[mapping$island==d$region & mapping$original_island_index<=d$analysis_end,,drop=FALSE]
    if(!setequal(mapping$original_island_index,seq_len(d$analysis_end))) stop("Raw mapping does not cover analysis rows",call.=FALSE)
    x$provenance$raw_mapping <- mapping
    x$provenance$manifest <- manifest
  }
  x
}
p3_check <- function(cfg,data) {
  p3_text(cfg$campaign_id,"campaign_id")
  if(!grepl("^[A-Za-z0-9_-]+$",cfg$campaign_id)) stop("Use a filename-safe campaign_id",call.=FALSE)
  if(!is.character(cfg$analysis_type) || length(cfg$analysis_type)!=1L ||
     !cfg$analysis_type %in% c("strict_future","conditional_future","descriptive")) stop("Set analysis_type",call.=FALSE)
  validate_phase2_data(data)
  if(!is.logical(cfg$require_complete_years) || length(cfg$require_complete_years)!=1L || is.na(cfg$require_complete_years))
    stop("Set require_complete_years explicitly",call.=FALSE)
  split <- p3_splits(data,cfg$splits,cfg$require_complete_years)
  if(!is.null(split) && is.null(data$provenance$manifest))
    stop("Temporal campaigns require data.preprocessing_manifest with prepared hashes and raw row mapping",call.=FALSE)
  if(is.null(split) && cfg$analysis_type!="descriptive") stop("No splits: descriptive only",call.=FALSE)
  if(!is.null(split) && cfg$analysis_type=="descriptive") stop("Use no splits for a descriptive-only campaign",call.=FALSE)
  p3_candidates(cfg$candidates)
  scalar_number(cfg$selection$tie_tolerance,"tie tolerance",0)
  if(!identical(cfg$selection$rule,"smaller_K_then_config_order")) stop("Unknown tie rule",call.=FALSE)
  if(!identical(cfg$peaks$ties,"first_chronological") || !cfg$peaks$incomplete_years %in% c("exclude","include"))
    stop("Set supported peak policy",call.=FALSE)
  ends <- if(is.null(split)) length(data$cases) else max(split$indices$A)
  provenance <- p3_provenance(cfg,ends,forecast=!is.null(split))
  train <- p3_prefix(data,ends)
  if(length(cfg$initial_parameters$free_coordinates) &&
     any(!is.finite(cfg$initial_parameters$free_coordinates) | cfg$initial_parameters$free_coordinates!=0))
    stop("Phase 3 base shape starts must be zero; use declared shape perturbations",call.=FALSE)
  checked_candidates <- cfg$candidates
  if(is.null(split)) {
    choice <- cfg$descriptive_choice
    if(is.null(choice)) stop("Descriptive-only check requires descriptive_choice",call.=FALSE)
    checked_candidates <- data.frame(id="declared_descriptive",K=choice$K,lambda=choice$lambda)
    p3_candidates(checked_candidates)
  }
  for(i in seq_len(nrow(checked_candidates))) {
    model <- p3_model(cfg,checked_candidates[i,,drop=FALSE])
    problem <- prepare_estimation_problem(train,model,cfg$estimation)
    initial <- cfg$initial_parameters;initial$free_coordinates <- rep(0,model$K-1L)
    encode_parameters(initial,problem$map)
    # Check every grid endpoint against guards without evaluating dynamics or RNG.
    for(beta in cfg$starts$beta0_candidates) for(phi in cfg$starts$phi_candidates) {
      initial$beta0 <- beta;initial$phi <- phi
      encode_parameters(initial,problem$map)
    }
  }
  validate_starts_config(cfg$starts)
  if(length(cfg$starts$user_starts)) stop("Phase 3 uses matched beta0/phi grids only; user_starts must be NULL",call.=FALSE)
  if(length(cfg$estimation$ndeps)!=1L || length(cfg$estimation$parscale)!=1L)
    stop("Phase 3 needs scalar ndeps/parscale shared across different active dimensions",call.=FALSE)
  m0 <- p3_model(cfg,constant=TRUE)
  problem0 <- prepare_estimation_problem(train,m0,cfg$estimation)
  initial0 <- cfg$initial_parameters;initial0$free_coordinates <- NULL
  encode_parameters(initial0,problem0$map)
  if(length(cfg$starts$beta0_candidates)*length(cfg$starts$phi_candidates)*(1+cfg$starts$n_perturbed)+
     length(cfg$starts$user_starts)>32L) stop("Phase 3 caps each candidate at 32 starts",call.=FALSE)
  if(!is.null(split)) {
    p3_text(cfg$final_output,"final_output")
    p3_future_input(data,split$indices$B,cfg,"B",max(split$indices$A))
    p3_future_input(data,split$indices$C,cfg,"C",max(split$indices$B))
    p3_information(cfg,max(split$indices$A))
    p3_evaluation_policy(cfg,data,"B",max(split$indices$A),split$indices$B)
    p3_evaluation_policy(cfg,data,"C",max(split$indices$B),split$indices$C)
  }
  list(type="phase3_check",status="CONFIG_CHECKED_NO_FIT",splits=split,provenance=provenance,
    population=data$fixed_parameters,markers=p3_markers(data),mask=data$likelihood_include,
    data_audit=p3_data_audit(data,cfg),raw_mapping=data$provenance$raw_mapping,config=cfg)
}
