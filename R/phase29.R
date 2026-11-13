# Phase 29 scores saved jointly fitted trajectories without optimization.
p29_read <- function(path) {
  if (!file.exists(path)) stop("Missing input: ",path)
  readRDS(path)
}
p29_weekly <- function(cfg) {
  m0 <- p29_read(cfg$m0_fit)
  if (m0$result$type!="phase22_fit" || m0$result$primary$status!="ENDPOINTS_RECORDED")
    stop("Wrong Phase 22 joint fit")
  base <- m0$result$primary$attempts[[m0$result$primary$best]]$score
  w0 <- utils::read.csv(cfg$m0_weekly)
  w0 <- w0[w0$reference=="JOINT",]
  if (nrow(w0)!=length(base$observed) ||
      max(abs(w0$expected-base$states$expected_cases))>1e-6)
    stop("Phase 22 trajectory disagrees with saved fit")
  w0 <- data.frame(row=w0$row,year=w0$year,week=w0$week,
    observed=w0$observed,expected=w0$expected,
    in_likelihood=w0$likelihood_included)
  result <- list(SPL_M0=list(weekly=w0,phi=base$parameters$phi))
  for (key in c("m1","m2")) {
    fit <- p29_read(cfg[[paste0(key,"_fit")]])
    ev <- p29_read(cfg[[paste0(key,"_evaluation")]])
    model <- toupper(key)
    if (fit$metadata$phase!=26L || fit$metadata$stage!=paste0(key,"-fit") ||
        ev$metadata$phase!=26L || ev$metadata$stage!=paste0(key,"-evaluate") ||
        fit$result$model!=model || ev$result$model!=model ||
        !identical(ev$result$status,"EVALUATED"))
      stop("Phase 26 fit/evaluation identity mismatch: ",key)
    if (max(abs(fit$result$score$weekly$expected-ev$result$score$weekly$expected))>1e-7 ||
        !identical(fit$result$penalty$id,"primary")) stop("Phase 26 trajectory mismatch: ",key)
    result[[paste0("SPL_",model)]] <- list(weekly=ev$result$score$weekly,
      phi=ev$result$score$parameters$phi)
  }
  f <- p29_read(cfg$fourier_result)
  if (f$metadata$phase!=27L || f$metadata$stage!="fourier" ||
      f$result$status!="FOURIER_SAVED") stop("Wrong Phase 27 artifact")
  for (family in c("A","AS")) {
    key <- paste0("fourier_",family)
    fit <- p29_read(cfg[[paste0(key,"_fit")]])
    w <- utils::read.csv(cfg[[paste0(key,"_weekly")]])
    registered <- f$result$fit_objects[[paste0("raw_",family)]]
    if (is.null(registered) || !identical(fit$endpoint,registered$endpoint) ||
        fit$model!="M2" ||
        fit$family!=switch(family,A="FOURIER_1",AS="FOURIER_2") ||
        fit$penalty$id!="primary" ||
        nrow(w)!=nrow(fit$score$weekly) ||
        max(abs(w$expected-fit$score$weekly$expected))>1e-7)
      stop("Fourier fit/weekly identity mismatch: ",family)
    result[[paste0(family,"_M2")]] <- list(weekly=fit$score$weekly,
      phi=fit$score$parameters$phi)
  }
  result
}
p29_align <- function(models) {
  ref <- models[["SPL_M0"]]$weekly
  for (name in names(models)) {
    w <- models[[name]]$weekly
    if (nrow(w)!=nrow(ref) || any(w$row!=ref$row | w$year!=ref$year |
        w$week!=ref$week | w$observed!=ref$observed |
        w$in_likelihood!=ref$in_likelihood) ||
        !identical(as.logical(w$in_likelihood),c(FALSE,rep(TRUE,nrow(w)-1L))))
      stop("Common observation/calendar/mask mismatch: ",name)
  }
  if (sum(ref$year==2020L & ref$week==53L)!=1L) stop("Week 53 mismatch")
  invisible(TRUE)
}
p29_metric <- function(observed,expected,phi,mask) {
  a <- observed[mask];b <- expected[mask]
  if (!length(a)) stop("Empty scoring period")
  data.frame(n=length(a),mae=mean(abs(a-b)),rmse=sqrt(mean((a-b)^2)),
    wape=if(sum(a)==0) NA_real_ else sum(abs(a-b))/sum(a),
    wape_reason=if(sum(a)==0) "zero observed denominator" else "",
    nb2_loglik=nb2_loglik(observed,expected,phi,mask)$total)
}
p29_peak <- function(w,season,model) {
  ix <- which((w$year==season & w$week>=27L) |
    (w$year==season+1L & w$week<=26L))
  complete <- length(ix)==52L && all(w$in_likelihood[ix])
  if (!complete) return(data.frame(model=model,season=season,complete=FALSE,
    observed_peak=NA_integer_,fitted_peak=NA_integer_,signed_error=NA_integer_,
    absolute_error=NA_integer_,observed_ties="",fitted_ties="",
    observed_boundary=NA,fitted_boundary=NA))
  a <- which(w$observed[ix]==max(w$observed[ix]))
  b <- which(w$expected[ix]==max(w$expected[ix]))
  data.frame(model=model,season=season,complete=TRUE,observed_peak=a[1],
    fitted_peak=b[1],signed_error=b[1]-a[1],absolute_error=abs(b[1]-a[1]),
    observed_ties=paste(a,collapse=";"),fitted_ties=paste(b,collapse=";"),
    observed_boundary=any(a %in% c(1L,length(ix))),
    fitted_boundary=any(b %in% c(1L,length(ix))))
}
p29_shift <- function(peaks) {
  x <- peaks[match(c(2015L,2017L),peaks$season),]
  if (anyNA(x$observed_peak)) stop("Transition season unavailable")
  observed_change <- diff(x$observed_peak)
  fitted_change <- diff(x$fitted_peak)
  data.frame(model=x$model[1],observed_peak_2015=x$observed_peak[1],
    observed_peak_2017=x$observed_peak[2],fitted_peak_2015=x$fitted_peak[1],
    fitted_peak_2017=x$fitted_peak[2],endpoint_error_2015=x$signed_error[1],
    endpoint_error_2017=x$signed_error[2],observed_change=observed_change,
    fitted_change=fitted_change,shift_error=abs(fitted_change-observed_change),
    fitted_advance=-fitted_change)
}
p29_evaluate <- function(models) {
  p29_align(models)
  periods <- list(`2013-2019`=2013:2019,`2015-2017`=2015:2017,
    `2020-2023`=2020:2023,`2013-2023`=2013:2023)
  metrics <- do.call(rbind,lapply(names(models),function(name) {
    x <- models[[name]];w <- x$weekly
    do.call(rbind,lapply(names(periods),function(period) {
      z <- p29_metric(w$observed,w$expected,x$phi,
        as.logical(w$in_likelihood) & w$year %in% periods[[period]])
      data.frame(model=name,period=period,z)
    }))
  }))
  peaks <- do.call(rbind,lapply(names(models),function(name)
    do.call(rbind,lapply(2012:2019,function(y) p29_peak(models[[name]]$weekly,y,name)))))
  shifts <- do.call(rbind,lapply(split(peaks,peaks$model),p29_shift))
  summary <- do.call(rbind,lapply(split(peaks,peaks$model),function(x)
    data.frame(model=x$model[1],mean_absolute_peak_error=mean(x$absolute_error[x$complete]),
      maximum_absolute_peak_error=max(x$absolute_error[x$complete]),
      evaluated_seasons=sum(x$complete))))
  table <- function(labels) {
    do.call(rbind,lapply(labels,function(name) {
      m <- metrics[metrics$model==name,];s <- shifts[shifts$model==name,]
      p <- summary[summary$model==name,]
      h <- m[m$period=="2013-2019",];t <- m[m$period=="2015-2017",]
      data.frame(model=name,mae_2013_2019=h$mae,rmse_2013_2019=h$rmse,
        mae_2015_2017=t$mae,rmse_2015_2017=t$rmse,
        fitted_peak_advance=s$fitted_advance,shift_error=s$shift_error,
        historical_mean_absolute_peak_error=p$mean_absolute_peak_error,
        historical_maximum_absolute_peak_error=p$maximum_absolute_peak_error)
    }))
  }
  list(fitting_metrics=metrics,timing_shift=shifts,peak_timing_by_season=peaks,
    peak_timing_summary=summary,
    comparison_correction=table(c("SPL_M0","SPL_M1","SPL_M2")),
    comparison_seasonal_function=table(c("A_M2","AS_M2","SPL_M2")))
}
p29_figures <- function(models,result,out) {
  dir.create(file.path(out,"figures"))
  cols <- c(SPL_M0="#777777",SPL_M1="#A67928",SPL_M2="#1B6EAB",
    A_M2="#B44B49",AS_M2="#7758A6")
  draw <- function(name,fun,width=10,height=7) {
    for (ext in c("pdf","png")) {
      path <- file.path(out,"figures",paste0(name,".",ext))
      if(ext=="pdf") grDevices::pdf(path,width=width,height=height) else
        grDevices::png(path,width=width*120,height=height*120,res=120)
      fun();grDevices::dev.off()
    }
  }
  draw("case_reconstruction",function() {
    old <- graphics::par(mfrow=c(5,2),mar=c(2.5,3,2,1));on.exit(graphics::par(old))
    for(name in names(models)) for(period in c("2013-2019","2015-2017")) {
      w <- models[[name]]$weekly;years <- if(period=="2013-2019") 2013:2019 else 2015:2017
      ix <- which(w$year %in% years);ylim <- range(c(w$observed,w$expected),finite=TRUE)
      graphics::plot(ix,w$observed[ix],type="l",col="black",ylim=ylim,
        xlab="Week index",ylab="Reported cases",main=paste(name,period))
      graphics::lines(ix,w$expected[ix],col=cols[name],lwd=1.5)
    }
  })
  draw("timing_shift",function() {
    s <- result$timing_shift;old <- graphics::par(mar=c(5,7,3,2));on.exit(graphics::par(old))
    labels <- c("Observed",s$model);y <- rev(seq_along(labels));
    graphics::plot(NA,xlim=c(1,52),ylim=range(y)+c(-.5,.5),yaxt="n",
      xlab="Peak week within epidemic season",ylab="",main="2015 to 2017 peak timing")
    graphics::axis(2,at=y,labels=labels,las=1)
    for(i in seq_along(labels)) {
      x <- if(i==1) c(s$observed_peak_2015[1],s$observed_peak_2017[1]) else
        c(s$fitted_peak_2015[i-1],s$fitted_peak_2017[i-1])
      col <- if(i==1) "black" else cols[labels[i]]
      graphics::segments(x[1],y[i],x[2],y[i],col=col,lwd=2)
      graphics::points(x,y=rep(y[i],2),pch=c(16,17),col=col)
      graphics::text(max(x)+2,y[i],labels=paste0("change ",x[2]-x[1]),pos=4,cex=.7)
    }
    graphics::legend("topright",c("2015","2017"),pch=c(16,17),bty="n")
  })
  draw("epidemic_timing",function() {
    x <- result$peak_timing_by_season;old <- graphics::par(mfrow=c(5,1),mar=c(3,4,2,1));on.exit(graphics::par(old))
    ymax <- max(x$absolute_error,na.rm=TRUE)
    for(name in names(models)) {
      z <- x[x$model==name,];graphics::plot(z$season,z$absolute_error,type="b",
        col=cols[name],pch=16,ylim=c(0,max(1,ymax)),xlim=c(2012,2019),
        xlab="Season start year",ylab="Absolute error",main=name)
      special <- z$complete & (grepl(";",z$observed_ties) |
        grepl(";",z$fitted_ties) | z$observed_boundary | z$fitted_boundary)
      graphics::points(z$season[special],z$absolute_error[special],pch=1,cex=1.6)
      if(any(!z$complete)) graphics::text(z$season[!z$complete],0,"NA",pos=3)
    }
    graphics::mtext("Open circles: tied or boundary peak",side=1,outer=TRUE,line=-1)
  })
}
