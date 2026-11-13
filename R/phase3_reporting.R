# Tables and base-R graphics are an explicit post-fit action, never a load side effect.
p3_report <- function(out,run_id,candidates=NULL,final=NULL,descriptive=NULL) {
  artifacts <- Filter(Negate(is.null),list(candidates,final,descriptive))
  if(!length(artifacts)) stop("Report needs at least one completed Phase 3 artifact")
  cfg <- artifacts[[1]]$config
  for(a in artifacts) if(!identical(a$config,cfg) || !identical(a$identity,artifacts[[1]]$identity))
    stop("Do not combine different campaigns/configurations/inputs")
  for(a in artifacts) if(!is.null(a$models)) p3_matched(a$models$M0,a$models$M1)
  if(!is.null(candidates$baseline)) for(r in candidates$records) p3_matched(candidates$baseline,r)
  prefix <- paste0(run_id,"-")
  csv <- function(x,name) utils::write.csv(x,file.path(out,paste0(prefix,name,".csv")),row.names=FALSE)
  pdf <- function(name,fun) {
    grDevices::pdf(file.path(out,paste0(prefix,name,".pdf")),width=10,height=6)
    on.exit(grDevices::dev.off(),add=TRUE);fun()
  }
  base <- if(!is.null(final)) final else if(!is.null(descriptive)) descriptive else candidates$check
  splits <- base$splits
  marks <- if(!is.null(final)) final$markers else if(!is.null(descriptive)) descriptive$markers else candidates$check$markers
  context <- paste(cfg$data$region,paste(marks$year[1],marks$week[1],sep="-W"),
    "to",paste(tail(marks$year,1),tail(marks$week,1),sep="-W"),"|",run_id)
  summary <- if(is.null(splits)) data.frame(split="full descriptive",first_index=1L,last_index=nrow(marks),n=nrow(marks)) else
    do.call(rbind,lapply(names(splits$indices),function(n) {
      i <- splits$indices[[n]];data.frame(split=n,first_index=min(i),last_index=max(i),n=length(i))
    }))
  summary$first_marker <- paste(marks$year[summary$first_index],marks$week[summary$first_index],sep="-W")
  summary$last_marker <- paste(marks$year[summary$last_index],marks$week[summary$last_index],sep="-W")
  summary$region <- cfg$data$region;summary$aggregation <- cfg$data$aggregation
  summary$population <- if(!is.null(candidates)) candidates$check$population$pop else {
    available <- artifacts[[1]]$models
    f <- Filter(function(x) !is.null(x$fit$representative),available)
    if(length(f)) f[[1]]$fit$representative$comparison_context$fixed$pop else NA_real_
  }
  summary$first_global_observation_excluded <- TRUE
  csv(summary,"data-splits")
  rows <- list();statuses <- list()
  for(stage in c("final","descriptive")) {
    a <- if(stage=="final") final else descriptive
    if(is.null(a)) next
    for(n in names(a$models)) {
      x <- a$models[[n]];r <- x$fit$representative
      statuses[[length(statuses)+1L]] <- data.frame(stage=stage,model=n,status=x$status,
        reason=if(is.null(x$reason)) "" else x$reason,
        comparison="Matched M0/M1 assumptions; primary comparison requires both AVAILABLE; local optimizer endpoints only")
      if(!is.null(r)) {
        t <- r$parameter_record;t$model <- n;t$stage <- stage
        rows[[length(rows)+1L]] <- t
      }
    }
  }
  if(length(rows)) csv(do.call(rbind,rows),"model-parameters") else csv(data.frame(status="No available fitted parameter records"),"model-parameters")
  if(length(statuses)) csv(do.call(rbind,statuses),"model-comparability") else csv(data.frame(status="No final/descriptive fits supplied"),"model-comparability")
  provenance <- artifacts[[1]]$information_policy
  if(!is.null(provenance)) csv(provenance,"fixed-value-provenance")
  if(!is.null(candidates)) {
    t <- candidates$table;csv(t,"candidate-B-scores")
    # Starts and endpoints are saved in RDS; export key failure/objective diagnostics too.
    starts <- list()
    for(rec in candidates$records) if(length(rec$fit$runs)) for(i in seq_along(rec$fit$runs)) {
      r <- rec$fit$runs[[i]];e <- r$evaluation
      starts[[length(starts)+1L]] <- data.frame(candidate=rec$candidate$id,start=i,status=r$status,
        objective=if(is.null(e)) NA_real_ else e$objective,
        loglik=if(is.null(e)) NA_real_ else e$loglik,penalty=if(is.null(e)) NA_real_ else e$penalty,
        message=if(is.null(r$message)) "" else paste(r$message,collapse="; "),
        boundary=paste(r$boundary$coordinates,collapse="; "))
    }
    if(length(starts)) csv(do.call(rbind,starts),"candidate-start-diagnostics")
    pdf("candidate-B-scores",function() {
      good <- t$status=="AVAILABLE" & is.finite(t$mean_nb2_nll)
      if(any(good)) {
        lim <- range(t$mean_nb2_nll[good]);pad <- max(diff(lim)*.15,.01)
        graphics::plot(seq_len(nrow(t))[good],t$mean_nb2_nll[good],xlim=c(.5,nrow(t)+.5),
          ylim=lim+c(-pad,pad),xaxt="n",pch=16,xlab="Explicit K/lambda candidate",ylab="B mean plug-in NB2 NLL",
          main=paste("Candidate validation |",candidates$evaluation_policy$class),sub=context)
        if(any(!good)) {
          graphics::points(which(!good),rep(lim[1]-pad/2,sum(!good)),pch=4,col="red")
          graphics::text(which(!good),rep(lim[1]-pad/2,sum(!good)),labels=t$status[!good],pos=3,cex=.6,col="red")
        }
        graphics::axis(1,seq_len(nrow(t)),paste0(t$id," [K=",t$K,", L=",t$lambda,"]"),cex.axis=.7)
      } else {graphics::plot.new();graphics::title(main="All candidate scores unavailable",sub=context)}
    })
  } else csv(data.frame(status="UNAVAILABLE: no candidate artifact supplied"),"candidate-B-scores")
  final_metrics <- list();descriptive_metrics <- list();peak_rows <- list();series <- list()
  for(stage in c("final","descriptive")) {
    a <- if(stage=="final") final else descriptive
    if(is.null(a)) next
    for(n in names(a$models)) {
      x <- a$models[[n]]
      if(!identical(x$status,"AVAILABLE")) next
      label <- if(stage=="final") {
        if(is.null(x$evaluation_class)) "UNAVAILABLE_CLASSIFICATION" else x$evaluation_class
      } else "descriptive"
      scores <- if(stage=="final") x$forecast$scores else x$scores
      metrics <- scores$summary;metrics$model <- n;metrics$analysis <- label
      if(stage=="final") metrics$evaluation_class <- if(is.null(x$evaluation_class)) "UNAVAILABLE" else x$evaluation_class
      if(stage=="final") final_metrics[[n]] <- metrics
      if(stage=="descriptive") descriptive_metrics[[n]] <- metrics
      peaks <- if(stage=="final") x$forecast$peaks else x$peaks
      peaks$model <- n;peaks$analysis <- label;peak_rows[[length(peak_rows)+1L]] <- peaks
      e <- x$fit$representative$evaluation
      expected <- if(stage=="final") c(e$states$expected_cases,x$forecast$states$expected_cases) else e$states$expected_cases
      t <- cbind(index=seq_len(nrow(a$markers)),a$markers,observed=a$observed,expected=expected)
      t$model <- n;t$analysis <- label
      t$scored <- if(stage=="final") t$index %in% x$forecast$indices else a$likelihood_include
      series[[length(series)+1L]] <- t
    }
  }
  if(length(final_metrics)) csv(do.call(rbind,final_metrics),"final-C-metrics") else
    csv(data.frame(status="UNAVAILABLE: no successful final evaluation supplied; not held-out evidence"),"final-C-metrics")
  if(length(descriptive_metrics)) csv(do.call(rbind,descriptive_metrics),"descriptive-metrics")
  if(length(series)) {
    all <- do.call(rbind,series);csv(all,"cases-expected-series")
    for(stage in unique(all$analysis)) pdf(paste0(stage,"-cases"),function() {
      d <- all[all$analysis==stage,,drop=FALSE]
      graphics::plot(d$index,d$observed,pch=16,cex=.35,col="grey60",xlab="Original prepared observation index (52 weeks/year)",
        ylab="Reported cases/week",main=paste("Cases and expected cases |",stage),sub=context,
        ylim=range(c(d$observed,d$expected)))
      for(n in unique(d$model)) {s <- d[d$model==n,];graphics::lines(s$index,s$expected,col=if(n=="M0") "blue" else "red",lwd=1.5)}
      if(!is.null(splits)) graphics::abline(v=c(max(splits$indices$A),max(splits$indices$B))+.5,lty=2)
      models <- unique(d$model);colors <- c(M0="blue",M1="red")
      graphics::legend("topright",c("observed",paste(models,"expected")),col=c("grey60",colors[models]),
        lty=c(NA,rep(1,length(models))),pch=c(16,rep(NA,length(models))),bty="n")
    })
  }
  if(length(peak_rows)) {
    peaks <- do.call(rbind,peak_rows);csv(peaks,"yearly-peaks")
    pdf("yearly-peaks",function() {
      p <- peaks[peaks$eligible,,drop=FALSE]
      if(!nrow(p)) {graphics::plot.new();graphics::title("No eligible complete-year peaks",sub=context);return(invisible(NULL))}
      graphics::plot(p$year,p$observed_peak_week,pch=16,ylim=c(1,52),xlab="Prepared year",ylab="Peak seasonal week (1..52)",
        main="Observed vs expected case peaks; first chronological ties",sub=context)
      groups <- interaction(p$model,p$analysis,drop=TRUE)
      for(i in seq_along(levels(groups))) {
        d <- p[groups==levels(groups)[i],];graphics::points(d$year,d$expected_peak_week,pch=1+i,col=1+i)
      }
      graphics::legend("topright",c("observed",levels(groups)),pch=c(16,seq_along(levels(groups))+1),col=c(1,seq_along(levels(groups))+1),bty="n",cex=.7)
    })
  } else csv(data.frame(status="No available case-peak comparison"),"yearly-peaks")
  chosen <- if(!is.null(final) && identical(final$models$M1$status,"AVAILABLE")) final$models$M1 else
    if(!is.null(descriptive) && identical(descriptive$models$M1$status,"AVAILABLE")) descriptive$models$M1 else NULL
  if(is.null(chosen) && !is.null(candidates)) {
    records <- Filter(function(r) identical(r$candidate$id,candidates$selection$id) && identical(r$status,"AVAILABLE"),candidates$records)
    if(length(records)) chosen <- records[[1]]
  }
  if(!is.null(chosen)) {
    e <- chosen$fit$representative$evaluation
    beta <- data.frame(week=1:52,beta=e$beta_by_week,shape=e$beta_by_week/e$parameters$beta0)
    csv(beta,"selected-beta-shape")
    pdf("selected-beta-shape",function() {
      graphics::par(mfrow=c(2,1),mar=c(4,4,3,1))
      label <- if(!is.null(final)) final$evaluation_class else if(!is.null(descriptive)) "descriptive" else "A fit selected on B"
      graphics::plot(beta$week,beta$beta,type="l",xlab="Week w (coordinate w-1)",ylab="beta (week^-1)",main=paste("M1 transmission |",label))
      graphics::plot(beta$week,beta$shape,type="l",xlab="Prepared seasonal week",ylab="Shape (weekly mean 1)",main=context)
    })
  }
  ref <- if(is.null(descriptive)) list(status="UNAVAILABLE",reason="No descriptive artifact supplied") else descriptive$reference
  csv(data.frame(status=ref$status,reason=ref$reason,forecast_eligible=FALSE),"Mref-status")
  if(!is.null(ref$audit)) csv(ref$audit,"Mref-audit")
  if(identical(ref$status,"AVAILABLE_DESCRIPTIVE_ONLY")) {
    ids <- which(descriptive$likelihood_include)
    # No NB2 score: pointwise medians and no coherent joint plug-in phi.
    d <- data.frame(model="Mref",analysis="descriptive_pointwise_posterior_medians",
      MAE=mean(abs(descriptive$observed[ids]-ref$expected[ids])),
      RMSE=sqrt(mean((descriptive$observed[ids]-ref$expected[ids])^2)),
      observed_total=sum(descriptive$observed[ids]),expected_total=sum(ref$expected[ids]),
      mean_nb2_nll=NA_real_,reason="No joint plug-in parameters; Bayesian inference differs")
    csv(d,"Mref-descriptive")
    csv(cbind(index=seq_along(ref$expected),descriptive$markers,observed=descriptive$observed,
      expected=ref$expected,analysis="descriptive_pointwise_posterior_medians"),"Mref-cases-expected")
    csv(p3_peaks(descriptive$markers[ids,,drop=FALSE],descriptive$observed[ids],ref$expected[ids],cfg$peaks),"Mref-yearly-peaks")
    csv(data.frame(week=1:52,beta=ref$beta,definition="pointwise posterior medians; no mean-one normalization"),"Mref-beta")
  }
  p3_report_outputs(out,run_id,candidates,final,descriptive)
  list(type="phase3_report",status="REPORT_WRITTEN_NO_FITTING",config=cfg,identity=artifacts[[1]]$identity,
    files=list.files(out,recursive=TRUE),run_id=run_id,interpretation="Matched M0/M1 forecasts and Mref descriptive reference are separate; no uncertainty integration")
}
