phase27_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  stages <- c("replay","mechanism","sensitivity","annual","fourier","validation","report")
  if(length(args)<3L || !args[1L] %in% stages || args[2L]!="--out" ||
     (length(args)>3L && (length(args)-3L)%%2L))
    stop("Usage: Rscript scripts/phase27.R STAGE --out NEW_DIR [--replay PATH] [--mechanism PATH] [--sensitivity PATH] [--annual PATH] [--fourier PATH] [--validation PATH]")
  root <- normalizePath(file.path(dirname(sub("^--file=","",
    grep("^--file=",commandArgs(),value=TRUE)[1L])),".."))
  setwd(root)
  sys.source("R/load_phase6.R",globalenv());load_spline_sirs_phase6(root,globalenv())
  for(f in c("R/phase5_reporting.R","R/phase10_peak_identity.R",
    "R/result_storage.R","R/accepted_extension.R","R/phase21.R","R/phase22.R",
    "R/phase22_analysis.R","R/phase24.R","R/phase25.R","R/phase26.R",
    "R/phase16_targeted.R","R/phase27.R"))
    sys.source(f,globalenv())
  for(f in c("config/phase22.R","config/phase24.R","config/phase16.R",
    "config/phase25.R","config/phase26.R","config/phase27.R"))
    sys.source(f,globalenv())
  c22 <- phase22_config();c25 <- phase25_config()
  c26 <- phase26_config();c27 <- phase27_config()
  if(args[1L] %in% c("fourier","validation")) {
    c26$adaptive_gradient <- TRUE
    c26$refine_descent <- TRUE
    c26$refine_gradient_step <- c27$refine_gradient_step
    c26$refine_maxit <- c27$refine_maxit
    c26$descent_repair_max_steps <- c27$descent_repair_max_steps
  }
  rest <- args[-(1:3)]
  keys <- if(length(rest)) rest[seq.int(1L,length(rest),by=2L)] else character()
  vals <- if(length(rest)) rest[seq.int(2L,length(rest),by=2L)] else character()
  if(any(!keys %in% paste0("--",stages)) || anyDuplicated(keys)) stop("Invalid options")
  opt <- as.list(setNames(vals,sub("^--","",keys)))
  stage <- args[1L];out <- reserve_phase2_output(args[3L])
  log_path <- file.path(out,"execution.log")
  log_event <- function(event,detail="") {
    line <- paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S %z"),
      stage,event,detail,sep=" | ")
    cat(line,"\n",file=log_path,append=TRUE,sep="")
    message(line)
  }
  old_progress <- getOption("phase27.progress")
  options(phase27.progress=log_event)
  on.exit(options(phase27.progress=old_progress),add=TRUE)
  writeLines(c("RUNNING",stage),file.path(out,"STATUS"))
  log_event("START",paste("output",out))
  tryCatch({
  log_event("INPUTS_START")
  result <- {
    if(stage!="replay") {
      if(is.null(opt$replay)) stop("Empirical stage requires --replay")
      gate <- read_phase2_result(opt$replay)
      if(!identical(gate$metadata$phase,27L) ||
         !identical(gate$metadata$stage,"replay") ||
         !identical(gate$result$status,"REPLAY_PASSED")) stop("Phase 27 replay gate missing")
    }
    inputs <- p22_inputs(c22)
    ref <- p27_reference(inputs,c22,c26,c27)
    if(stage=="replay") p27_replay(ref,inputs,c26,c27,out) else
    if(stage=="mechanism") p27_mechanism(ref,inputs,c27,c25,out) else
    if(stage=="sensitivity") p27_sensitivity(ref,inputs,c27,c22,out) else
    if(stage=="annual") p27_annual(ref,inputs,c22,c27,out) else
    if(stage=="fourier") p27_fit_fourier(ref,inputs,c26,c27,out) else
    if(stage=="validation") p27_validation(ref,inputs,c26,c27,c22,out) else {
      needed <- c("mechanism","sensitivity","annual","fourier","validation")
      paths <- unlist(opt[needed])
      if(length(paths)!=length(needed) || anyNA(paths)) stop("Report requires all empirical stage paths")
      paths <- setNames(paths,needed)
      p27_report(ref,c27,c22,paths,out)
    }
  }
  log_event("COMPUTATION_COMPLETED")
  sources <- c(file.path("config",paste0("phase",c(22,24,25,26,27),".R")),
    unlist(c22[c("reference","replay","panel","plan","fit","evaluation","joint_fit")]),
    c27$fit,c27$evaluation,c27$phase25,c27$joint_sensitivity,c27$joint_matched,
    c27$joint_annual,c27$joint_fourier,c27$joint_validation,
    c27$mechanism_grid,c22$annual_profiles,c22$annual_scenarios,
    c22$annual_windows,unlist(opt))
  metadata <- phase2_metadata(root,"config/phase27.R",sources)
  metadata$phase <- 27L;metadata$stage <- stage
  metadata$reference_status <- c27$reference_status
  metadata$endpoint_parameter_md5 <- p27_parameter_hash(ref$fit$endpoint)
  metadata$configuration <- c27
  p27_write(data.frame(path=names(metadata$input_hashes),
    md5=unname(metadata$input_hashes)),out,"provenance")
  dput(list(stage=stage,reference_status=c27$reference_status,
    endpoint_fit=c27$fit,endpoint_fit_md5=unname(tools::md5sum(c27$fit)),
    endpoint_parameter_hash=p27_parameter_hash(ref$fit$endpoint),
    code_commit=metadata$code_commit,config=c27,code_hashes=metadata$code_hashes,
    input_hashes=metadata$input_hashes),file=file.path(out,"manifest.R"))
  write_phase2_result(result,out,metadata)
  log_event("STORAGE_COMPLETED")
  invisible(result)
  },error=function(e) {
    log_event("FAILED",conditionMessage(e))
    if(!identical(readLines(file.path(out,"STATUS"),n=1L),"STORAGE_FAILED"))
      writeLines(c("RUN_FAILED",conditionMessage(e)),file.path(out,"STATUS"))
    stop(e)
  })
}
if(sys.nframe()==0L) phase27_main()
