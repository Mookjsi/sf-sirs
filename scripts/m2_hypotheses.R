m2_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  stages <- c("inventory","plan","summarize-h1","h1-full","summarize-h2","state-memory","report","bundle")
  if(length(args)<3L || !args[1L] %in% stages || args[2L]!="--out" ||
     (length(args)-3L)%%2L) stop("Usage: Rscript scripts/m2_hypotheses.R STAGE --out NEW_DIR [--replay DIR --mechanism DIR --sensitivity DIR] [--h1 DIR --h2 DIR --inventory DIR --state DIR]")
  script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1L])
  root <- normalizePath(file.path(dirname(script),".."));setwd(root)
  sys.source("R/load_phase6.R",globalenv());load_spline_sirs_phase6(root,globalenv())
  for(f in c("R/phase5_reporting.R","R/phase10_peak_identity.R","R/result_storage.R",
    "R/accepted_extension.R","R/phase21.R","R/phase22.R","R/phase22_analysis.R",
    "R/phase24.R","R/phase25.R","R/phase26.R","R/phase16_targeted.R",
    "R/phase27.R","R/m2_hypotheses.R")) sys.source(f,globalenv())
  for(f in c("config/phase22.R","config/phase24.R","config/phase16.R",
    "config/phase25.R","config/phase26.R","config/phase27.R")) sys.source(f,globalenv())
  rest <- args[-(1:3)];keys <- rest[seq.int(1L,length(rest),by=2L)]
  vals <- rest[seq.int(2L,length(rest),by=2L)]
  allowed <- c("--replay","--mechanism","--sensitivity","--h1","--h2","--inventory","--state","--report","--resume","--reuse-cells")
  if(any(!keys %in% allowed) || anyDuplicated(keys)) stop("Invalid options")
  opt <- as.list(setNames(vals,sub("^--","",keys)))
  stage <- args[1L]
  if(!is.null(opt$resume) && !opt$resume %in% c("yes","no")) stop("--resume must be yes or no")
  if(identical(opt$resume,"yes")) {
    out <- normalizePath(args[3L],mustWork=TRUE)
    if(!identical(readLines(file.path(out,"STATUS"),n=1L),"RUN_FAILED"))
      stop("Only a failed output directory may be resumed")
  } else out <- reserve_phase2_output(args[3L])
  log <- file.path(out,"execution.log")
  event <- function(label,detail="") {
    line <- paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S %z"),stage,label,detail,sep=" | ")
    cat(line,"\n",file=log,append=TRUE);message(line)
  }
  writeLines("RUNNING",file.path(out,"STATUS"));event("START")
  tryCatch({
    if(!is.null(opt[["reuse-cells"]])) {
      source <- opt[["reuse-cells"]]
      if(stage!="state-memory" || !identical(readLines(file.path(source,"STATUS"),n=1L),
        "STORAGE_COMPLETED")) stop("--reuse-cells requires a completed state-memory source")
      cells <- list.files(file.path(source,"cells"),pattern="^h3-.*[.]rds$",
        full.names=TRUE)
      if(!length(cells)) stop("No state-memory cells to reuse")
      dir.create(file.path(out,"cells"),showWarnings=FALSE)
      if(!all(file.copy(cells,file.path(out,"cells")))) stop("Could not copy reusable cells")
      event("COPIED_NUMERICAL_CELLS",paste(length(cells),"from",source))
    }
    if(stage=="bundle") {
      if(is.null(opt$report) || !identical(readLines(file.path(opt$report,"STATUS"),n=1L),"STORAGE_COMPLETED"))
        stop("Bundle requires a completed --report directory")
      files <- sort(list.files(opt$report,full.names=TRUE,recursive=TRUE))
      hashes <- tools::md5sum(files[file.info(files)$isdir==FALSE])
      write.csv(data.frame(file=substring(names(hashes),nchar(normalizePath(opt$report))+2L),
        md5=unname(hashes)),
        file.path(out,"bundle_manifest.csv"),row.names=FALSE)
      target <- file.path(out,"m2_results.tar.gz")
      args <- c("--sort=name","--mtime=@0","--owner=0","--group=0","--numeric-owner",
        "-czf",target,"-C",normalizePath(opt$report),".")
      if(system2("tar",args)!=0L) stop("Could not create results bundle")
    } else if(stage=="report") {
      required <- c("inventory","h1","h2")
      if(any(!required %in% names(opt))) stop("Report requires --inventory, --h1 and --h2")
      files <- list(inventory=c("existing_experiment_inventory"),
        h1=c("classification_audit","h1_event_summary"),
        h2=c("h2_magnitude_summary","h2_direction_summary","h2_step_audit"))
      for(k in names(files)) {
        if(!identical(readLines(file.path(opt[[k]],"STATUS"),n=1L),"STORAGE_COMPLETED"))
          stop("Incomplete analysis stage: ",k)
        for(f in files[[k]]) file.copy(file.path(opt[[k]],paste0(f,".csv")),out)
      }
      for(f in c("h1_season_summary","h1_equal_season_summary",
        "h1_movement_by_mechanism","season_coverage"))
        if(file.exists(file.path(opt$h1,paste0(f,".csv"))))
          file.copy(file.path(opt$h1,paste0(f,".csv")),out)
      if(file.exists(file.path(opt$h2,"h2_season_summary.csv")))
        file.copy(file.path(opt$h2,"h2_season_summary.csv"),out)
      reference <- utils::read.csv(file.path(opt$inventory,"existing_experiment_inventory.csv"))
      reference <- reference[reference$selected,]
      for(k in c("h1","h2")) {
        check <- utils::read.csv(file.path(opt[[k]],"existing_experiment_inventory.csv"))
        check <- check[check$selected,]
        if(!identical(check$parameter_hash,reference$parameter_hash) ||
           !identical(check$fit_md5,reference$fit_md5))
          stop("Report inputs use different M2 endpoints")
      }
      if(!is.null(opt$state) && !identical(readLines(file.path(opt$state,"STATUS"),n=1L),"STORAGE_COMPLETED"))
        stop("Incomplete state-memory stage")
      if(!is.null(opt$state)) {
        check <- utils::read.csv(file.path(opt$state,"existing_experiment_inventory.csv"))
        check <- check[check$selected,]
        if(!identical(check$parameter_hash,reference$parameter_hash))
          stop("State-memory output uses a different M2 endpoint")
      }
      m2_report(out,out,opt$state)
      for(k in c("inventory","h1","h2","state")) if(!is.null(opt[[k]])) {
        destination <- file.path(out,"sources",k)
        dir.create(destination,recursive=TRUE,showWarnings=FALSE)
        source_files <- list.files(opt[[k]],pattern="[.]csv$",full.names=TRUE)
        if(length(source_files) && !all(file.copy(source_files,destination)))
          stop("Could not copy report source tables: ",k)
      }
    } else {
      needed <- c("replay","mechanism","sensitivity")
      if(any(!needed %in% names(opt))) stop("Stage requires --replay, --mechanism and --sensitivity")
      evidence <- m2_inventory(opt[needed],out)
      if(stage=="plan") m2_plan(evidence,out)
      if(stage=="summarize-h1") m2_h1(evidence,out)
      if(stage=="h1-full") m2_h1_full(evidence,out,progress=event)
      if(stage=="summarize-h2") m2_h2(evidence,out)
      if(stage=="state-memory") m2_state_memory(evidence,out,progress=event)
    }
    writeLines("STORAGE_COMPLETED",file.path(out,"STATUS"));event("STORAGE_COMPLETED")
  },error=function(e) {
    writeLines(c("RUN_FAILED",conditionMessage(e)),file.path(out,"STATUS"))
    event("FAILED",conditionMessage(e));stop(e)
  })
}
if(sys.nframe()==0L) m2_main()
