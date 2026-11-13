# Plug-in observation scores, not posterior predictive densities.
p3_circular_distance <- function(a,b) {
  a <- seasonal_weeks(a);b <- seasonal_weeks(b)
  if(length(a)!=length(b)) stop("Peak vector lengths differ",call.=FALSE)
  d <- abs(a-b);pmin(d,52L-d)
}
p3_metrics <- function(y,expected,phi) {
  finite_numeric(y,"evaluation cases");finite_numeric(expected,"expected cases",length(y))
  ll <- nb2_loglik(y,expected,phi,rep(TRUE,length(y)))
  list(summary=data.frame(n=length(y),mean_nb2_nll=-mean(ll$per_week),
    MAE=mean(abs(y-expected)),RMSE=sqrt(mean((y-expected)^2)),bias=mean(expected-y),
    observed_total=sum(y),expected_total=sum(expected)),
    weekly=data.frame(observed=y,expected=expected,nb2_logp=ll$per_week,included=TRUE),
    interpretation="Plug-in NB2, fitted parameters fixed; no parameter uncertainty integration; bias=expected-observed")
}
p3_peaks <- function(markers,y,expected,policy=list(incomplete_years="exclude",ties="first_chronological")) {
  if(!all(c("year","week") %in% names(markers))) stop("Peak markers need year/week",call.=FALSE)
  finite_numeric(y,"peak observations",nrow(markers));finite_numeric(expected,"peak expected",nrow(markers))
  seasonal_weeks(markers$week)
  if(!identical(policy$ties,"first_chronological") || !policy$incomplete_years %in% c("include","exclude")) stop("Unknown peak policy")
  rows <- lapply(unique(markers$year),function(yr) {
    ids <- which(markers$year==yr);w <- markers$week[ids]
    complete <- identical(as.integer(w),seq_len(52L))
    eligible <- complete || policy$incomplete_years=="include"
    oi <- which(y[ids]==max(y[ids]));ei <- which(expected[ids]==max(expected[ids]))
    ow <- w[oi[1]];ew <- w[ei[1]]
    data.frame(year=yr,complete_year=complete,eligible=eligible,
      observed_peak_week=if(eligible) ow else NA_integer_,
      expected_peak_week=if(eligible) ew else NA_integer_,
      circular_error=if(eligible) p3_circular_distance(ow,ew) else NA_integer_,
      signed_error=if(eligible) ew-ow else NA_integer_,
      absolute_error=if(eligible) abs(ew-ow) else NA_integer_,
      signed_circular_error=if(eligible) ((ew-ow+26) %% 52)-26 else NA_integer_,
      evaluated_weeks=length(ids),excluded_weeks=52L-length(unique(w)),
      observed_ties=length(oi),expected_ties=length(ei),
      reason=if(eligible) "first chronological maximum" else "incomplete year excluded")
  })
  do.call(rbind,rows)
}
