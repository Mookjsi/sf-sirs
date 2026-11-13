# Stage-specific orchestration. No action is called by source().
p3_inputs <- function(cfg) {
  manifest <- cfg$data$preprocessing_manifest
  mapping <- if(is.null(manifest)) NULL else file.path(dirname(manifest),"prepared_row_mapping.csv")
  unique(c(cfg$data$cases_file,cfg$data$population_file,cfg$data$raw_input_paths,manifest,mapping))
}
p3_identity <- function(root,cfg,config_path) {
  files <- unique(c(list.files(file.path(root,"R"),pattern="[.]R$",full.names=TRUE),
    list.files(file.path(root,"scripts"),pattern="[.](R|sh)$",full.names=TRUE),
    file.path(root,"config/phase2_defaults.R"),config_path,p3_inputs(cfg)))
  if(any(!file.exists(files))) stop("Provenance input missing")
  tools::md5sum(normalizePath(files))
}
p3_same_campaign <- function(artifact,cfg,identity) {
  if(!identical(artifact$config,cfg) || !identical(artifact$identity,identity))
    stop("Campaign config/code/input hashes changed; preserve this campaign and create a separate one",call.=FALSE)
  invisible(TRUE)
}
p3_candidate_stage <- function(cfg,data,check,identity,optimizer=stats::optim) {
  if(is.null(check$splits)) stop("Candidates need explicit A/B/C splits")
  ids <- check$splits$indices
  train <- p3_prefix(data,max(ids$A))
  future <- p3_future_input(data,ids$B,cfg,"B",max(ids$A))
  # Only A observations reach the fitting engine; B is scoring-only, C is absent.
  records <- p3_candidate_campaign(train,future,p3_outcomes(data,ids$B),cfg,optimizer)
  baseline <- p3_baseline_record(train,future,p3_outcomes(data,ids$B),cfg,optimizer)
  for(record in records) p3_matched(baseline,record)
  table <- p3_candidate_table(records);selection <- p3_choose(table,cfg$selection)
  list(type="phase3_candidates",status=selection$status,config=cfg,identity=identity,
    check=check,training_indices=ids$A,records=records,table=table,selection=selection,
    baseline=baseline,evaluation_policy=p3_evaluation_policy(cfg,data,"B",max(ids$A),ids$B),
    data_audit=check$data_audit,analysis_type=cfg$analysis_type,information_policy=check$provenance)
}
p3_lock <- function(campaign,cfg,identity) {
  p3_same_campaign(campaign,cfg,identity)
  if(!identical(campaign$type,"phase3_candidates")) stop("Lock requires a candidate artifact")
  choice <- p3_choose(campaign$table,cfg$selection)
  if(is.null(choice$id)) stop("All candidates unavailable: selection cannot be locked")
  if(file.exists(cfg$final_output)) stop("Locked final output must not exist")
  # Absolute path under an existing parent prevents later --out substitution.
  parent <- normalizePath(dirname(cfg$final_output),mustWork=TRUE)
  destination <- file.path(parent,basename(cfg$final_output))
  list(type="phase3_lock",status="SELECTION_LOCKED_NO_C_SCORE",config=cfg,identity=identity,
    choice=choice,candidate_table=campaign$table,final_output=destination,
    splits=campaign$check$splits,information_policy=campaign$information_policy,
    evaluation_policy=campaign$evaluation_policy,data_audit=campaign$data_audit,
    assumptions=list(model=cfg$model,estimation=cfg$estimation,starts=cfg$starts,
      initial_parameters=cfg$initial_parameters,future=cfg$future,information=cfg$information),
    analysis_type=cfg$analysis_type,selection_uses="A fits and B scores only")
}
p3_evaluate_model <- function(data,train_end,ids,cfg,candidate=NULL,constant=FALSE,optimizer=stats::optim) {
  result <- list(model=if(constant) "M0" else "M1",status="FAILED",fit=NULL,forecast=NULL,reason=NULL)
  tryCatch({
    result$evaluation_policy <- p3_evaluation_policy(cfg,data,"C",train_end,ids)
    result$evaluation_class <- result$evaluation_policy$class
    result$information_policy <- p3_provenance(cfg,train_end,TRUE)
    train <- p3_prefix(data,train_end)
    result$fit <- p3_fit(train,cfg,candidate,constant,optimizer)
    if(is.null(result$fit$representative)) {
      result$status <- result$fit$status;result$evaluation_class <- "UNAVAILABLE"
      result$reason <- "No converged valid A+B refit";return(result)
    }
    input <- p3_future_input(data,ids,cfg,"C",train_end)
    result$forecast <- p3_score_forecast(p3_forecast(train,result$fit$representative,input),
      p3_outcomes(data,ids),cfg$peaks)
    if(!is.finite(result$forecast$scores$summary$mean_nb2_nll)) stop("Nonfinite final score")
    result$status <- "AVAILABLE";result
  },error=function(e) {result$status <- "FAILED";result$evaluation_class <- "UNAVAILABLE";result$reason <- conditionMessage(e);result})
}
p3_final <- function(lock,cfg,data,identity,optimizer=stats::optim) {
  p3_same_campaign(lock,cfg,identity)
  if(!identical(lock$type,"phase3_lock")) stop("Final evaluation requires a locked selection")
  ids <- lock$splits$indices;end <- max(ids$B)
  # No feedback, retry, hyperparameter selection or automatic C-dependent refit.
  models <- list(M0=p3_evaluate_model(data,end,ids$C,cfg,constant=TRUE,optimizer=optimizer),
                 M1=p3_evaluate_model(data,end,ids$C,cfg,lock$choice$candidate,optimizer=optimizer))
  p3_matched(models$M0,models$M1)
  list(type="phase3_final",status=if(all(vapply(models,function(x) identical(x$status,"AVAILABLE"),logical(1))))
    "FINAL_EVALUATION_STORED" else "FINAL_EVALUATION_HAS_FAILURES",config=cfg,identity=identity,
    lock=lock,models=models,splits=lock$splits,analysis_type=cfg$analysis_type,
    observed=data$cases,markers=p3_markers(data),information_policy=lock$information_policy,
    data_audit=p3_data_audit(data,cfg),raw_mapping=data$provenance$raw_mapping,
    evaluation_class=if(all(vapply(models,function(x) identical(x$status,"AVAILABLE"),logical(1))))
      models$M1$evaluation_class else "UNAVAILABLE",
    interpretation="One fixed final destination per lock; no C-driven model selection")
}
p3_descriptive <- function(cfg,data,identity,lock=NULL,optimizer=stats::optim) {
  if(is.null(lock)) {
    if(!is.null(cfg$splits)) stop("A split campaign's descriptive fit requires its lock")
    declared <- cfg$descriptive_choice
    if(is.null(declared)) stop("Set descriptive_choice or provide a lock")
    candidate <- data.frame(id="declared_descriptive",K=declared$K,lambda=declared$lambda)
    p3_candidates(candidate)
  } else {
    p3_same_campaign(lock,cfg,identity)
    if(!identical(lock$type,"phase3_lock")) stop("Descriptive stage requires lock artifact")
    candidate <- lock$choice$candidate
  }
  provenance <- p3_provenance(cfg,length(data$cases),forecast=FALSE)
  fits <- lapply(c("M0","M1"),function(n) tryCatch({
    fit <- p3_fit(data,cfg,candidate,constant=n=="M0",optimizer=optimizer)
    if(is.null(fit$representative)) return(list(model=n,status=fit$status,fit=fit,reason="No converged valid descriptive fit"))
    e <- fit$representative$evaluation;ids <- which(data$likelihood_include)
    list(model=n,status="AVAILABLE",fit=fit,reason="",indices=ids,
      scores=p3_metrics(data$cases[ids],e$states$expected_cases[ids],e$parameters$phi),
      peaks=p3_peaks(p3_markers(data,ids),data$cases[ids],e$states$expected_cases[ids],cfg$peaks))
  },error=function(e) list(model=n,status="FAILED",fit=NULL,reason=conditionMessage(e))))
  names(fits) <- c("M0","M1")
  p3_matched(fits$M0,fits$M1)
  list(type="phase3_descriptive",status="DESCRIPTIVE_ONLY_NOT_HELD_OUT",config=cfg,identity=identity,
    models=fits,candidate=candidate,markers=p3_markers(data),observed=data$cases,
    likelihood_include=data$likelihood_include,information_policy=provenance,
    reference=p3_reference(cfg$reference,data),analysis_type="descriptive",
    data_audit=p3_data_audit(data,cfg),raw_mapping=data$provenance$raw_mapping,
    splits=if(is.null(lock)) NULL else lock$splits)
}
p3_metadata <- function(root,config_path,cfg,paths,run_id) {
  md <- phase2_metadata(root,config_path,unique(c(p3_inputs(cfg),paths)),run_id)
  md$phase <- 3L;md$code_hashes <- p3_identity(root,cfg,config_path)
  md$source_preprocessing <- cfg$data$preprocessing_id
  md$campaign_id <- cfg$campaign_id
  md$observation_model <- "NB2 mean=reported expected incidence, size=phi"
  md
}
