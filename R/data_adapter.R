# Accept original preprocessed island tables, not unprocessed surveillance rows.
original_region_definition <- function(region) {
  if(!is.character(region) || length(region)!=1L || is.na(region) ||
     !region %in% c("Hokkaido","Honshu","Shikoku","Kyushu","Ryukyu"))
    stop("Specify exactly one original island explicitly",call.=FALSE)
  n <- switch(region,Hokkaido=208L,Ryukyu=312L,364L)
  list(region=region,variant=if(region=="Shikoku") "sirs_npi_tweak" else "sirs_npi",
       Nnonpi=n,npi_starts= switch(region,Hokkaido="2017 week 1",Ryukyu="2019 week 1","2020 week 1"),
       independent_npi_normal_prior=(region!="Shikoku"),
       source_commit="475c7dcdadc12987033d332b101f025625074301")
}
original_npi_effect <- function(n,Nnonpi,post_npi=NULL) {
  integer_scalar(n,"n"); integer_scalar(Nnonpi,"Nnonpi")
  free <- max(0L,n-Nnonpi)
  if(free && is.null(post_npi)) stop("Post-baseline NPI multipliers are unset; supply them explicitly",call.=FALSE)
  if(!free && !is.null(post_npi) && length(post_npi)) stop("No free NPI weeks in this input",call.=FALSE)
  if(free) {
    finite_numeric(post_npi,"post_npi",free)
    if(any(post_npi<0)) stop("NPI multipliers must be nonnegative",call.=FALSE)
  }
  c(rep(1,min(n,Nnonpi)),if(free) post_npi else numeric())
}
adapt_preprocessed_island <- function(cases,populations,region,post_npi=NULL) {
  definition <- original_region_definition(region)
  if(!is.data.frame(cases) || !all(c("island","year","week","cases") %in% names(cases)))
    stop("Cases require original output columns: island, year, week, cases",call.=FALSE)
  if(!is.data.frame(populations) || !all(c("island","pop") %in% names(populations)))
    stop("Population requires original output columns: island, pop",call.=FALSE)
  if(anyNA(cases[,c("island","year","week","cases")]) || anyNA(populations[,c("island","pop")]))
    stop("Missing source values are not imputed",call.=FALSE)
  d <- cases[as.character(cases$island)==region,,drop=FALSE]
  p <- populations[as.character(populations$island)==region,,drop=FALSE]
  if(!nrow(d) || nrow(p)!=1L) stop("Region is absent or population is not a unique island denominator",call.=FALSE)
  finite_numeric(d$year,"year"); finite_numeric(d$cases,"cases")
  w <- seasonal_weeks(d$week)
  if(any(d$year!=floor(d$year)) || any(d$cases<0 | d$cases!=floor(d$cases)) ||
     any(d$cases>.Machine$integer.max)) stop("Invalid year or integer case counts",call.=FALSE)
  key <- 52*d$year + w
  if(anyDuplicated(key) || any(diff(key)!=1)) stop("Input must already be unique, consecutive and chronologically ordered; no sorting/aggregation is performed",call.=FALSE)
  if(d$year[1]!=2013L || w[1]!=1L)
    stop("Original Nnonpi indexing requires an initial state at 2013 week 1; cropped series need a separate explicit model decision",call.=FALSE)
  scalar_number(p$pop,"population",0,strict_lower=TRUE)
  n <- nrow(d)
  structure(list(observation_order=seq_len(n),time_markers=d,
       processed_week_marker=paste(d$year,w,sep="-W"),seasonal_week=w,
       cases=as.numeric(d$cases),fixed_parameters=list(pop=as.numeric(p$pop),mu=1/80/52,gamma=1),
       npieff=original_npi_effect(n,definition$Nnonpi,post_npi),
       likelihood_include=c(FALSE,rep(TRUE,n-1L)),definition=definition,
       provenance=list(input_level="already aggregated island/week",
         week53="Original loader maps 53 to 52, then rounds mean counts; not reimplemented here",
         original_raw_week53_markers="Unavailable in this prepared schema; time_markers preserves all supplied columns")),
    class="spline_sirs_data")
}
read_preprocessed_island <- function(cases_file,population_file,region,post_npi=NULL) {
  adapt_preprocessed_island(utils::read.csv(cases_file,stringsAsFactors=FALSE),
    utils::read.csv(population_file,stringsAsFactors=FALSE),region,post_npi)
}
