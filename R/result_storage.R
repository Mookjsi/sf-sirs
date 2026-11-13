# Provenance and guarded storage; a computed fit is not labelled saved until writes finish.
reserve_phase2_output <- function(path) {
  if(length(path)!=1L || !nzchar(path) || file.exists(path)) stop("Use an explicit new output directory",call.=FALSE)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  if(!dir.create(path)) stop("Cannot reserve output directory",call.=FALSE)
  normalizePath(path)
}
phase2_metadata <- function(root,config_path,input_paths=character(),run_id=NULL) {
  if(is.null(run_id)) run_id <- paste0(format(Sys.time(),"%Y%m%dT%H%M%S%z"),"-",Sys.getpid())
  commit <- tryCatch(suppressWarnings(system2("git",c("-C",shQuote(root),"rev-parse","HEAD"),
    stdout=TRUE,stderr=FALSE)),error=function(e) character())
  if(length(commit)!=1L || !grepl("^[0-9a-f]{40}$",commit)) commit <- "UNAVAILABLE: no verified Git commit"
  code <- c(list.files(file.path(root,"R"),pattern="[.]R$",full.names=TRUE),
             list.files(file.path(root,"scripts"),pattern="[.]R$",full.names=TRUE),
             list.files(file.path(root,"config"),pattern="^phase2.*[.]R$",full.names=TRUE),config_path)
  list(run_id=run_id,code_commit=commit,hash_algorithm="MD5 file integrity, not a commit identifier",
       code_hashes=tools::md5sum(code),config_path=normalizePath(config_path),
       input_paths=input_paths,input_hashes=if(length(input_paths)) tools::md5sum(input_paths) else character(),
       R=R.version.string,session=capture.output(utils::sessionInfo()),RNG=RNGkind(),
       method="penalized maximum likelihood; not Bayesian posterior or MAP")
}
validate_phase2_metadata <- function(metadata) {
  required <- c("run_id","code_commit","hash_algorithm","code_hashes","config_path","input_paths",
                "input_hashes","R","session","RNG","method")
  if(!is.list(metadata) || !all(required %in% names(metadata)))
    stop("Missing required provenance fields for result storage",call.=FALSE)
  for(n in c("run_id","code_commit","config_path","R","method"))
    if(!is.character(metadata[[n]]) || length(metadata[[n]])!=1L || is.na(metadata[[n]]) || !nzchar(metadata[[n]]))
      stop("Invalid provenance field ",n,call.=FALSE)
  invisible(metadata)
}
write_phase2_result <- function(result,out,metadata) {
  validate_phase2_metadata(metadata)
  if(!dir.exists(out) || file.exists(file.path(out,"result.rds"))) stop("Output is not reserved or already contains a result",call.=FALSE)
  temporary <- file.path(out,"result.rds.partial")
  if(file.exists(temporary)) stop("Partial output exists; preserve it and choose another directory",call.=FALSE)
  tryCatch({
    saveRDS(list(schema_version=1L,metadata=metadata,result=result),temporary)
    if(!file.rename(temporary,file.path(out,"result.rds"))) stop("Atomic result rename failed")
    capture.output(dput(metadata),file=file.path(out,"metadata.txt"))
    writeLines(c("STORAGE_COMPLETED",paste("COMPUTATION_STATUS",if(is.null(result$status)) result$type else result$status)),
               file.path(out,"STATUS"))
  },error=function(e) {
    try(writeLines(c("STORAGE_FAILED",conditionMessage(e)),file.path(out,"STATUS")),silent=TRUE)
    stop("Storage failed; computed result is not certified as saved: ",conditionMessage(e),call.=FALSE)
  })
  invisible(file.path(out,"result.rds"))
}
read_phase2_result <- function(path) {
  if(!file.exists(file.path(dirname(path),"STATUS")) ||
     !identical(readLines(file.path(dirname(path),"STATUS"),n=1L),"STORAGE_COMPLETED"))
    stop("Artifact storage completion marker is absent",call.=FALSE)
  x <- readRDS(path)
  if(!identical(x$schema_version,1L)) stop("Unsupported result schema",call.=FALSE)
  validate_phase2_metadata(x$metadata)
  x
}
