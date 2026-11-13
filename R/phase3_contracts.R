# Information timing and matched-comparison contracts; no fitting on source.
p3_source_row <- function(source,name,origin,strict) {
  if(!is.list(source)) stop("Missing information provenance for ",name,call.=FALSE)
  p3_text(source$citation,paste(name,"citation"))
  if(!is.character(source$category) || length(source$category)!=1L ||
     !source$category %in% c("external_independent","training_only","full_data","user_conditional"))
    stop("Invalid information category for ",name,call.=FALSE)
  integer_scalar(source$last_used_index,paste(name,"last_used_index"),0L)
  integer_scalar(source$known_at_index,paste(name,"known_at_index"),0L)
  if(source$category=="external_independent" && source$last_used_index!=0L)
    stop("Independent information must use last_used_index=0: ",name,call.=FALSE)
  conditional <- source$category %in% c("full_data","user_conditional") ||
    source$last_used_index>origin || source$known_at_index>origin
  if(source$category=="training_only" && source$last_used_index>origin)
    stop("Training-only information reaches beyond origin: ",name,call.=FALSE)
  if(strict && conditional) stop("Strict forecast information unavailable at origin: ",name,call.=FALSE)
  data.frame(input=name,category=source$category,citation=source$citation,
    last_used_index=source$last_used_index,known_at_index=source$known_at_index,
    conditional=conditional)
}
p3_information <- function(cfg,origin) {
  strict <- identical(cfg$analysis_type,"strict_future")
  rows <- do.call(rbind,lapply(c("preprocessing","initialization"),function(n)
    p3_source_row(cfg$information[[n]],n,origin,strict)))
  # Conditional nuisance assumptions do not authorize outcome-adaptive starts/scaling.
  if(any(rows$last_used_index>origin))
    stop("Preprocessing/initialization must not use held-out outcomes, including in conditional mode",call.=FALSE)
  rows
}
p3_evaluation_policy <- function(cfg,data,stage,origin,ids) {
  fixed <- p3_provenance(cfg,origin,TRUE)
  information <- p3_information(cfg,origin)
  future <- p3_future_input(data,ids,cfg,stage,origin)
  conditional <- any(fixed$conditional_or_leakage) || any(information$conditional) ||
    isTRUE(future$source$conditional)
  list(class=if(conditional) "CONDITIONAL_FORECAST" else "STRICT_FORECAST",
    origin=origin,stage=stage,fixed=fixed,information=information,future=future$source,
    basis="Recorded source categories and availability indices; citations require user audit")
}
p3_matched <- function(a,b) {
  # Compare only the shared scientific assumptions, never K/lambda or objectives.
  if(is.null(a$fit$representative) || is.null(b$fit$representative)) return(invisible(FALSE))
  x <- a$fit$representative$comparison_context;y <- b$fit$representative$comparison_context
  for(n in c("cases","weeks","npieff","mask","markers","fixed","free"))
    if(!identical(x[[n]],y[[n]])) stop("Unmatched M0/M1 context: ",n,call.=FALSE)
  for(n in c("S0","I0","rho","omega","delta"))
    if(!identical(x$model[[n]],y$model[[n]])) stop("Unmatched nuisance parameter: ",n,call.=FALSE)
  if(!is.null(a$forecast) && !is.null(b$forecast) &&
     (!identical(a$forecast$indices,b$forecast$indices) ||
      !identical(a$forecast$markers,b$forecast$markers))) stop("Unmatched forecast observations",call.=FALSE)
  if(!is.null(a$forecast$scores) && !is.null(b$forecast$scores))
    for(n in c("observed","included"))
      if(!identical(a$forecast$scores$weekly[[n]],b$forecast$scores$weekly[[n]]))
        stop("Unmatched forecast scoring set: ",n,call.=FALSE)
  invisible(TRUE)
}
p3_data_audit <- function(data,cfg) {
  text_or_na <- function(x) if(is.null(x) || !length(x)) NA_character_ else x
  data.frame(region=data$definition$region,rows=length(data$cases),
    likelihood_rows=sum(data$likelihood_include),first_marker=data$processed_week_marker[1],
    last_marker=tail(data$processed_week_marker,1),missing_counts=sum(is.na(data$cases)),
    duplicated_weeks=anyDuplicated(data$processed_week_marker),case_unit="reported_cases_per_week",
    population=data$fixed_parameters$pop,mu=data$fixed_parameters$mu,gamma=data$fixed_parameters$gamma,
    preprocessing_id=cfg$data$preprocessing_id,cases_file=text_or_na(cfg$data$cases_file),
    population_file=text_or_na(cfg$data$population_file),aggregation=cfg$data$aggregation,
    week53="Mapped to 52 upstream; rounded prefecture means then island sums",
    raw_mapping="See prepared_row_mapping.csv when preprocessing manifest is supplied",
    stringsAsFactors=FALSE)
}
