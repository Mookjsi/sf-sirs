# Fixed-origin forecast from the saved Phase 26 M2 endpoint.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=8L || !identical(args[c(3L,5L,7L)],
  c('--reference','--hash','--population')))
  stop('Usage: Rscript scripts/forcast/forecast.R FIT_RESULT_RDS OUTPUT_DIR --reference REFERENCE_CSV --hash REFERENCE_HASH_TXT --population POPULATION_CSV')
fit_path <- args[1L]; out <- args[2L]
reference_path <- args[4L];hash_path <- args[6L];pop_path <- args[8L]
for(path in c(fit_path,reference_path,hash_path,pop_path))
  if(!file.exists(path)) stop('Missing required input: ',path)
if(!file.exists(fit_path)) stop('Missing M2 fit: ',fit_path)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source('R/validation.R')
source('R/sirs_core.R')
x <- readRDS(fit_path)$result
if(x$model!='M2' || x$status!='FIT_SAVED') stop('Expected completed Phase 26 M2 fit')
z <- x$score$weekly
stopifnot(nrow(z)==573L,tail(z$year,1)==2023L,tail(z$week,1)==52L,
  tail(z$observed,1)==120,abs(tail(z$expected,1)-115.1777)<1e-3)
p <- x$score$parameters
fixed <- list(pop=sum(tail(z[,c('S','I','R')],1)),mu=p$mu,gamma=p$gamma,
  omega=p$omega,delta=p$delta,rho=p$rho)
beta <- vapply(1:52,function(w) z$seasonal_beta[which(z$seasonal_week==w)[1L]],numeric(1))
year <- c(rep(2024L,52),rep(2025L,52),rep(2026L,53));week <- c(1:52,1:52,1:53)
seasonal <- pmin(week,52L)
multiplier <- rep(tail(z$future_q,1),length(week))
state <- unlist(tail(z[,c('S','I','R')],1),use.names=TRUE)
sim <- simulate_weekly_sirs(beta,state,fixed,c(52L,seasonal),c(1,multiplier),initial_mass=sum(state))
forecast <- data.frame(origin='2023-W52',target=sprintf('%d-W%02d',year,week),year=year,week=week,
  lead=seq_along(week),model='M2',mean=sim$expected_cases[-1L],
  lower80=NA_real_,upper80=NA_real_,lower95=NA_real_,upper95=NA_real_,
  future_scenario='2023-W52 terminal post-2019 multiplier held fixed')
reference <- read.csv(reference_path)
if(nrow(reference)!=66L || !identical(forecast$target[1:66],reference$target) ||
   max(abs(forecast$mean[1:66]-reference$mean))>1e-6)
  stop('First 66 point forecasts failed reference reproduction')
digest <- strsplit(system2('sha256sum',fit_path,stdout=TRUE),' ')[[1L]][1L]
if(!identical(readLines(hash_path,warn=FALSE),paste('fit_sha256',digest)))
  stop('M2 fit SHA-256 differs from reference forecast')
write.csv(forecast,file.path(out,'forecast_weekly.csv'),row.names=FALSE,na='')
population <- read.csv(pop_path)
fixed$pop <- sum(population$pop[population$island=='Honshu'])
fit_denominator <- simulate_weekly_sirs(beta,state,fixed,c(52L,seasonal),
  c(1,multiplier),initial_mass=sum(state))
sensitivity <- forecast
sensitivity$mean <- fit_denominator$expected_cases[-1L]
sensitivity$population_denominator <- fixed$pop
sensitivity$scenario <- 'fixed fitting population denominator'
write.csv(sensitivity,file.path(out,'forecast_fit_denominator_sensitivity.csv'),row.names=FALSE,na='')
manifest_path <- file.path(out,'source_manifest.csv')
if(file.exists(manifest_path)) {
  manifest <- read.csv(manifest_path,stringsAsFactors=FALSE)
  manifest <- rbind(manifest,data.frame(path=normalizePath(fit_path),kind='Phase 26 M2 fit through 2023-W52',sha256=digest))
  write.csv(manifest,manifest_path,row.names=FALSE,na='')
  writeLines(paste('fit_sha256',digest),file.path(out,'forecast_fit_hash.txt'))
}

cat('Fixed forecast saved:',nrow(forecast),'weeks; endpoint expected',tail(z$expected,1),'\n')
