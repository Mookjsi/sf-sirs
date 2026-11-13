args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3L || !args[1] %in% c("evaluate","report") || args[2]!="--out")
  stop("Usage: Rscript scripts/phase29.R evaluate --out NEW_DIR | report --out NEW_DIR --input EVALUATION_DIR")
root <- normalizePath(file.path(dirname(sub("^--file=","",
  grep("^--file=",commandArgs(),value=TRUE)[1])),".."))
setwd(root)
for(f in c("R/validation.R","R/nb2.R","R/phase29.R","config/phase29.R"))
  sys.source(f,globalenv())
cfg <- phase29_config();out <- args[3]
if(file.exists(out)) stop("Output directory already exists: ",out)
if(args[1]=="report") {
  if(length(args)!=5L || args[4]!="--input") stop("report requires --input EVALUATION_DIR")
  source_dir <- args[5]
  required <- c("input_manifest","fitting_metrics","timing_shift",
    "peak_timing_by_season","peak_timing_summary","comparison_correction",
    "comparison_seasonal_function")
  paths <- file.path(source_dir,paste0(required,".csv"))
  if(any(!file.exists(paths))) stop("Incomplete evaluation directory")
  dir.create(out,recursive=TRUE)
  file.copy(paths,out)
  if(dir.exists(file.path(source_dir,"figures"))) {
    dir.create(file.path(out,"figures"))
    file.copy(list.files(file.path(source_dir,"figures"),full.names=TRUE),file.path(out,"figures"))
  }
  result <- setNames(lapply(paths,utils::read.csv),required)
} else {
  if(length(args)!=3L) stop("evaluate accepts only --out NEW_DIR")
  missing <- unlist(cfg)[!file.exists(unlist(cfg))]
  if(length(missing)) stop("Missing input(s): ",paste(missing,collapse=", "))
  dir.create(out,recursive=TRUE)
  models <- p29_weekly(cfg)
  result <- p29_evaluate(models)
  for(name in names(result)) utils::write.csv(result[[name]],file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
  manifest <- data.frame(source=names(cfg),path=unlist(cfg),md5=unname(tools::md5sum(unlist(cfg))),
    scoring="rows 2 onward; calendar years; seasons week 27 through next week 26; earliest maximum; raw boundary")
  utils::write.csv(manifest,file.path(out,"input_manifest.csv"),row.names=FALSE)
  p29_figures(models,result,out)
}
a <- result$comparison_correction;b <- result$comparison_seasonal_function
lines <- c("# Phase 29: five internal models","",
  "Retrospective reconstruction of reported cases and epidemic peak timing.","",
  "## Correction structure", "",capture.output(print(a,row.names=FALSE)),"",
  "## Seasonal function with M2 correction","",capture.output(print(b,row.names=FALSE)),"",
  "MAE and RMSE assess waveform reconstruction. The shift and seasonal peak errors assess timing.",
  "See fitting_metrics.csv for WAPE, NB2 log likelihood and later/full periods; timing_shift.csv for endpoint errors.",
  "These metrics support a simulation baseline choice. They do not establish causal validity, forecasting performance, or complete simulator validation.")
writeLines(lines,file.path(out,"phase29_report.md"))
