phase22_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L || !args[1L] %in% c("inspect","fit","score") || args[2L]!="--out")
    stop("Usage: Rscript scripts/phase22.R inspect|fit|score --out NEW_DIR")
  root <- normalizePath(file.path(dirname(sub("^--file=","",
    grep("^--file=",commandArgs(),value=TRUE)[1L])),".."))
  setwd(root)
  sys.source("R/load_phase6.R",globalenv())
  load_spline_sirs_phase6(root,globalenv())
  for(f in c("R/result_storage.R","R/accepted_extension.R","R/phase21.R","R/phase22.R"))
    sys.source(f,globalenv())
  env <- new.env(parent=baseenv());sys.source("config/phase22.R",env)
  cfg <- env$phase22_config()
  inputs <- p22_inputs(cfg)
  out <- reserve_phase2_output(args[3L])
  input_paths <- unlist(cfg[c("reference","replay","panel","plan","fit","evaluation")])
  if(args[1L]=="score") input_paths <- c(input_paths,joint_fit=cfg$joint_fit)
  metadata <- phase2_metadata(root,"config/phase22.R",input_paths)
  metadata$config <- cfg
  if(args[1L]=="inspect") {
    result <- list(type="phase22_inspect",status="VERIFIED",hashes=inputs$hashes,
      replay=inputs$replay$checks,candidate=inputs$candidate,
      basis_dimension=inputs$basis$dimension,calendar=inputs$calendar)
  } else if(args[1L]=="fit") {
    primary <- p22_fit(inputs,p22_problem(inputs),cfg)
    biology <- p22_fit(inputs,p22_problem(inputs,free_biology=TRUE),cfg,
      seed=cfg$seed+1L)
    result <- list(type="phase22_fit",status="ATTEMPTS_SAVED",
      primary=primary,biology_sensitivity=biology,hashes=inputs$hashes)
  } else {
    if(!file.exists(cfg$joint_fit)) stop("Run the Phase 22 fit stage before score")
    joint_saved <- read_phase2_result(cfg$joint_fit)
    joint_fit <- joint_saved$result
    if(!identical(joint_fit$type,"phase22_fit") ||
       !identical(joint_fit$status,"ATTEMPTS_SAVED") ||
       is.null(joint_fit$primary$attempts)) stop("Saved Phase 22 joint fit is incomplete")
    p <- p22_problem(inputs)
    hist <- inputs$ref$record$retained_coordinates[p$problem$map$order]
    frozen <- p22_score(hist,inputs$frozen$theta,inputs,p)
    historical <- p22_score(hist,rep(0,inputs$basis$dimension),inputs,p)
    primary_attempt <- joint_fit$primary$attempts[[joint_fit$primary$best]]
    joint <- primary_attempt$score
    joint_row <- data.frame(reference="JOINT",historical_loglik=joint$historical_loglik,
      extension_loglik=joint$extension_loglik)
    score_table <- rbind(data.frame(reference=c("HIST","FROZEN"),
        historical_loglik=c(historical$historical_loglik,frozen$historical_loglik),
        extension_loglik=c(historical$extension_loglik,frozen$extension_loglik)),joint_row)
    result <- list(type="phase22_score",status="SCORED",hashes=inputs$hashes,
      joint_fit_hash=unname(tools::md5sum(cfg$joint_fit)),
      scores=score_table,
      joint_diagnostics=data.frame(status=primary_attempt$status,
        optimizer_convergence=primary_attempt$convergence,
        objective=joint$objective,loglik=joint$loglik,
        seasonal_penalty=joint$seasonal_penalty,post_penalty=joint$post_penalty,
        replay_difference=primary_attempt$replay_difference,
        competitive_endpoints=length(joint_fit$primary$competitive)),
      frozen_adequacy=inputs$saved$evaluation$result$adequacy)
    utils::write.csv(result$scores,file.path(out,"period_scores.csv"),row.names=FALSE)
  }
  write_phase2_result(result,out,metadata)
  invisible(result)
}
if(sys.nframe()==0L) phase22_main()
