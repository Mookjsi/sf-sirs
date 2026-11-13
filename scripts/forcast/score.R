# Score complete Honshu weeks after the fixed forecast exists.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Usage: Rscript scripts/forcast/score.R OUTPUT_DIR')
out <- args[1L]
f <- read.csv(file.path(out,'forecast_weekly.csv'),stringsAsFactors=FALSE)
o <- read.csv(file.path(out,'honshu_external_cases.csv'),stringsAsFactors=FALSE)
n <- read.csv(file.path(out,'naive_2023.csv'))
h <- o[o$island=='Honshu',]
if(anyDuplicated(paste(h$year,h$week,h$prefecture))) stop('Duplicate observation key')
keys <- unique(h[,c('year','week')])
totals <- lapply(seq_len(nrow(keys)),function(i) {
  z <- h[h$year==keys$year[i] & h$week==keys$week[i],]
  complete <- nrow(z)==34L && length(unique(z$prefecture))==34L && all(!is.na(z$cases))
  data.frame(year=keys$year[i],week=keys$week[i],
    cases_observed=if(complete) sum(z$cases) else NA_real_,
    source_version=if(length(unique(z$source_version))==1L) z$source_version[1L] else 'mixed',
    observation_status=if(complete) 'complete' else 'excluded: missing Honshu report')
})
a <- do.call(rbind,totals)
a <- merge(f,a,by=c('year','week'),all.x=TRUE,sort=FALSE)
a <- merge(a,n,by='week',all.x=TRUE,sort=FALSE)
names(a)[names(a)=='cases'] <- 'cases_naive'
a <- a[order(a$lead),]
a$surveillance_regime <- ifelse(a$year<2025L | (a$year==2025L & a$week<=14L),
  'pre-2025-W15','post-2025-W15')
a$source_version[is.na(a$source_version)] <- 'unavailable'
a$observation_status[is.na(a$observation_status)] <- 'not yet observed'
a$absolute_error_M2 <- abs(a$mean-a$cases_observed)
a$absolute_error_naive <- abs(a$cases_naive-a$cases_observed)
a$raw_count_error_note <- ifelse(a$surveillance_regime=='post-2025-W15',
  'Includes surveillance network change','')
write.csv(a,file.path(out,'scored_weekly.csv'),row.names=FALSE,na='')
periods <- list('2024-W01-W52'=a$year==2024,
  '2025-W01-W14'=a$year==2025 & a$week<=14,
  '2025-W15-W52'=a$year==2025 & a$week>=15 & a$week<=52,
  '2026-W01-W44'=a$year==2026 & a$week<=44)
score_period <- function(label,selection) {
  z <- a[selection & !is.na(a$cases_observed),]
  if(!nrow(z)) return(NULL)
  peak <- function(values) paste(z$target[which(values==max(values))],collapse=';')
  first_peak <- function(values) z$week[which(values==max(values))[1L]]
  mae_m2 <- mean(z$absolute_error_M2);mae_naive <- mean(z$absolute_error_naive)
  data.frame(period=label,source_version=paste(unique(z$source_version),collapse=';'),
    surveillance_regime=unique(z$surveillance_regime)[1L],weeks=nrow(z),
    mae_M2=mae_m2,rmse_M2=sqrt(mean((z$mean-z$cases_observed)^2)),
    cumulative_bias_M2=sum(z$mean-z$cases_observed),
    mae_naive=mae_naive,rmse_naive=sqrt(mean((z$cases_naive-z$cases_observed)^2)),
    cumulative_bias_naive=sum(z$cases_naive-z$cases_observed),
    improvement=if(mae_naive==0) NA_real_ else 1-mae_m2/mae_naive,
    observed_peak_week=first_peak(z$cases_observed),M2_peak_week=first_peak(z$mean),
    naive_peak_week=first_peak(z$cases_naive),
    observed_peak_candidates=peak(z$cases_observed),M2_peak_candidates=peak(z$mean),
    naive_peak_candidates=peak(z$cases_naive),
    M2_peak_week_error=first_peak(z$mean)-first_peak(z$cases_observed),
    observed_peak_height=max(z$cases_observed),M2_peak_height=max(z$mean),
    peak_height_error=max(z$mean)-max(z$cases_observed),
    observed_total=sum(z$cases_observed),M2_total=sum(z$mean),naive_total=sum(z$cases_naive),
    peak_scope=if(label=='2024-W01-W52') 'complete calendar year' else 'maximum within observed period',
    interpretation=if(unique(z$surveillance_regime)=='post-2025-W15')
      'Raw count errors include surveillance network change' else 'Raw count errors')
}
metrics <- do.call(rbind,Map(score_period,names(periods),periods))
write.csv(metrics,file.path(out,'metrics.csv'),row.names=FALSE,na='')
for(yr in c(2024L,2025L,2026L)) {
  z <- a[a$year==yr & a$week<=if(yr==2026L) 44L else 52L,]
  png(file.path(out,sprintf('forecast_%d.png',yr)),width=1200,height=650)
  plot(z$week,z$cases_observed,type='l',lwd=2,xlab=paste(yr,'epidemiological week'),
    ylab='Honshu RSV reports',ylim=range(c(z$cases_observed,z$mean,z$cases_naive),na.rm=TRUE),col='black')
  lines(z$week,z$mean,col='#2166ac',lwd=2)
  lines(z$week,z$cases_naive,col='#b2182b',lwd=2)
  if(yr==2025L) abline(v=14.5,lty=2,col='gray40')
  legend('topright',c('Observed','M2 fixed forecast','2023 seasonal naive'),
    col=c('black','#2166ac','#b2182b'),lwd=2,bty='n')
  dev.off()
}
print(metrics)
