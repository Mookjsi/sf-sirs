# Phase 28 compares the selected M2 reconstruction with the original Park main model.
p28_write <- function(x,out,name) {
  utils::write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
  invisible(x)
}
p28_panel <- function(inputs,ref,cfg) {
  clock <- p26_clock(inputs)
  z <- ref$score$weekly
  if(nrow(clock)!=nrow(z) || any(clock$year!=z$year | clock$week!=z$week))
    stop("M2 panel and replay rows differ")
  if(anyDuplicated(paste(clock$year,clock$week)) ||
     sum(clock$year==2020L & clock$week==53L)!=1L ||
     !identical(z$in_likelihood,c(FALSE,rep(TRUE,nrow(z)-1L))))
    stop("Chronology, week 53, or likelihood mask differs")
  original <- utils::read.csv(cfg$original_cases)
  original <- original[original$island==cfg$region,,drop=FALSE]
  population <- utils::read.csv(cfg$original_population)
  population <- population$pop[population$island==cfg$region]
  if(length(population)!=1L || population!=inputs$replay$data$fixed_parameters$pop)
    stop("Original and M2 population denominators differ")
  key <- paste(clock$year,clock$week)
  source_key <- paste(original$year,original$week)
  matched <- match(key,source_key)
  alignment <- data.frame(row=clock$row,year=clock$year,week=clock$week,
    seasonal_week=clock$seasonal_week,observed=z$observed,
    original_prepared=original$cases[matched],
    included=z$in_likelihood,
    source_match=!is.na(matched),stringsAsFactors=FALSE)
  exclusions <- data.frame(row=clock$row[!z$in_likelihood],
    reason="initial state; no likelihood contribution")
  exclusions <- rbind(exclusions,data.frame(
    row=clock$row[clock$year==2020L & clock$week==53L],
    reason="present in common M2 panel; original Park preprocessing merged this week into week 52"))
  ordinary <- !(clock$year==2020L & clock$week %in% c(52L,53L)) &
    !is.na(matched)
  if(any(z$observed[ordinary]!=original$cases[matched[ordinary]]))
    stop("Material case-count disagreement with original prepared panel")
  if(any(is.na(matched) & clock$week!=53L))
    stop("Original prepared panel lacks a non-week-53 row")
  alignment$adaptation <- ifelse(clock$year==2020L &
      clock$week %in% c(52L,53L),
    "Original Park preprocessing merged week 53 into week 52; Phase 28 retains both M2 rows with the original 52-week beta lookup",
    "none")
  if(any(clock$seasonal_week<1L | clock$seasonal_week>52L) ||
     clock$seasonal_week[clock$week==53L]!=52L)
    stop("Unverified seasonal-week mapping")
  n_nonpi <- which(clock$year>=2020L)[1L]-1L
  if(!is.finite(n_nonpi) || n_nonpi<1L || n_nonpi>=nrow(clock))
    stop("Cannot derive pre-intervention boundary")
  list(alignment=alignment,exclusions=exclusions,pop=population,
    N=nrow(clock),Nnonpi=n_nonpi,clock=clock,observed=z$observed,
    mask=z$in_likelihood,source_rows=nrow(original))
}
p28_standata <- function(panel) list(N=panel$N,Nnonpi=panel$Nnonpi,
  week=as.integer(panel$clock$seasonal_week),cases=as.integer(panel$observed),
  mu=1/80/52,pop=panel$pop,gamma=1)
p28_reference_recurrence <- function(beta,initial,fixed,week,multiplier) {
  n <- length(week)
  x <- matrix(0,n,5L,dimnames=list(NULL,c("S","I","R","incidence","expected_cases")))
  x[1L,1:3] <- initial
  delta <- fixed$delta;mu <- fixed$mu;gamma <- fixed$gamma
  for(i in 2:n) {
    foi <- beta[week[i]]*(x[i-1L,"I"]+fixed$omega)/fixed$pop*multiplier[i]
    sout <- (1-exp(-(foi+mu)))*x[i-1L,"S"]
    stoi <- foi/(foi+mu)*sout
    iout <- (1-exp(-(gamma+mu)))*x[i-1L,"I"]
    itor <- gamma/(gamma+mu)*iout
    rout <- (1-exp(-(delta+mu)))*x[i-1L,"R"]
    rtos <- delta/(delta+mu)*rout
    x[i,"S"] <- x[i-1L,"S"]+rtos-sout+mu*fixed$pop
    x[i,"I"] <- x[i-1L,"I"]+stoi-iout
    x[i,"R"] <- x[i-1L,"R"]+itor-rout
    x[i,"incidence"] <- stoi
    x[i,"expected_cases"] <- stoi*fixed$rho
  }
  x
}
p28_verify_equations <- function(panel,ref,inputs,tolerance=1e-7) {
  beta <- vapply(1:52,function(w) ref$score$weekly$seasonal_beta[
    which(panel$clock$seasonal_week==w)[1L]],numeric(1))
  mult <- ref$score$weekly$effective_beta/ref$score$weekly$seasonal_beta
  initial <- unlist(ref$score$states[1L,c("S","I","R")],use.names=TRUE)
  p <- ref$score$parameters
  fixed <- list(pop=panel$pop,mu=1/80/52,gamma=1,
    omega=ref_or(p$omega,inputs$ref$parameters$omega),
    delta=ref_or(p$delta,inputs$ref$parameters$delta),rho=p$rho)
  if(!is.finite(fixed$omega) || !is.finite(fixed$delta))
    stop("Missing fixed nuisance coordinates")
  actual <- simulate_weekly_sirs(beta,initial,fixed,
    panel$clock$seasonal_week,mult)
  expected <- p28_reference_recurrence(beta,initial,fixed,
    panel$clock$seasonal_week,mult)
  fields <- colnames(expected)
  state_error <- max(abs(as.matrix(actual[,fields])-expected))
  ll_a <- nb2_loglik(panel$observed,actual$expected_cases,p$phi,panel$mask)$total
  ll_b <- sum(stats::dnbinom(panel$observed[panel$mask],
    mu=expected[panel$mask,"expected_cases"],size=p$phi,log=TRUE))
  data.frame(check=c("SIRS and reported cases","NB2 included likelihood"),
    error=c(state_error,abs(ll_a-ll_b)),tolerance=tolerance,
    passed=c(state_error<=tolerance,abs(ll_a-ll_b)<=tolerance))
}
p28_curve_metrics <- function(observed,fitted,mask,clock,model,summary) {
  periods <- list(full=rep(TRUE,length(observed)),
    `2013-2019`=clock$year<=2019L,
    `2015-2017`=clock$year>=2015L & clock$year<=2017L,
    `2020-2023`=clock$year>=2020L & clock$year<=2023L)
  for(y in 2013:2022) {
    ix <- (clock$year==y & clock$week>=27L) |
      (clock$year==y+1L & clock$week<=26L)
    if(sum(ix)>=52L) periods[[paste0("season_",y)]] <- ix
  }
  do.call(rbind,lapply(names(periods),function(k) {
    ix <- periods[[k]] & mask
    a <- observed[ix];b <- fitted[ix]
    data.frame(model=model,summary=summary,period=k,rows=length(a),
      mae=mean(abs(a-b)),rmse=sqrt(mean((a-b)^2)),
      wape=sum(abs(a-b))/sum(a),signed_burden_error=sum(b-a))
  }))
}
p28_timing <- function(observed,fitted,clock,model,summary) {
  rows <- lapply(2013:2022,function(y) {
    ix <- which((clock$year==y & clock$week>=27L) |
      (clock$year==y+1L & clock$week<=26L))
    complete <- length(ix)==if(y==2020L) 53L else 52L
    if(!complete) return(data.frame(model=model,summary=summary,season=y,
      complete=FALSE,observed_peak=NA_integer_,fitted_peak=NA_integer_,
      observed_ties=NA_character_,fitted_ties=NA_character_,
      observed_cog=NA_real_,fitted_cog=NA_real_,peak_error=NA_real_,
      cog_error=NA_real_,boundary_peak=NA))
    a <- observed[ix];b <- fitted[ix]
    ap <- which(a==max(a));bp <- which(b==max(b))
    ac <- sum(seq_along(a)*a)/sum(a);bc <- sum(seq_along(b)*b)/sum(b)
    data.frame(model=model,summary=summary,season=y,complete=TRUE,
      observed_peak=ap[1L],fitted_peak=bp[1L],
      observed_ties=paste(ap,collapse=";"),fitted_ties=paste(bp,collapse=";"),
      observed_cog=ac,fitted_cog=bc,peak_error=bp[1L]-ap[1L],
      cog_error=bc-ac,boundary_peak=any(bp %in% c(1L,length(ix))))
  })
  do.call(rbind,rows)
}
p28_inventory <- function(ref,panel,park_fit=NULL) {
  names <- ref$fit$parameter_table$parameter[seq_along(ref$fit$endpoint)]
  seasonal <- sum(grepl("^log_beta0$|^z[0-9]+$",names))
  correction <- sum(grepl("^future_q_|^level_|^sin_|^cos_",names))
  initial <- sum(grepl("^log_S_over_R$|^log_I_over_R$",names))
  observation <- sum(grepl("^log_phi$|^logit_rho$",names))
  if(sum(c(seasonal,correction,initial,observation))!=length(names))
    stop("Unclassified M2 free coordinates")
  park_correction <- panel$N-panel$Nnonpi
  data.frame(model=c("M2","PARK"),
    free_coordinates=c(length(names),52L+park_correction+8L),
    fixed_parameter_count=c(4L,2L),
    seasonal_coordinates=c(seasonal,52L),
    time_varying_correction_coordinates=c(correction,park_correction),
    initial_state_coordinates=c(initial,2L),
    observation_coordinates=c(observation,2L),
    other_epidemiological_coordinates=c(0L,2L),
    estimated_prior_or_penalty_hyperparameters=c(0L,2L),
    note=c("M2 free coordinates counted from saved endpoint; anchored level excluded; fixed mu, gamma, omega, delta",
      "PARK includes 52 beta, post-intervention npi, and S0, I0, rho, sigma, phi, tau, omega, sigma_npi; fixed mu and gamma"))
}
p28_clock_time <- function() as.numeric(proc.time()[[1L]]+proc.time()[[2L]])
p28_timing_row <- function(stage,model,start,prep,compile,fit,post) {
  data.frame(stage=stage,model=model,cpu_model=tryCatch(readLines("/proc/cpuinfo",
    warn=FALSE)[grep("model name",readLines("/proc/cpuinfo",warn=FALSE))[1L]],
    error=function(e) NA_character_),
    allocated_cores=if(model=="PARK") 4L else 1L,
    actual_thread_setting=if(model=="PARK") "four RStan chains on four cores" else
      "single R optimizer process",
    preparation_seconds=prep,compilation_seconds=compile,
    fit_seconds=fit,postprocess_seconds=post,
    total_seconds=as.numeric(difftime(Sys.time(),start,units="secs")),
    process_tree_cpu_seconds=NA_real_,peak_memory_bytes=NA_real_)
}
