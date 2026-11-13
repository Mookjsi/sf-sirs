# Stable staged outputs from saved evaluations only; no fitting or propagation.
p3_bind_rows <- function(rows) {
  rows <- Filter(Negate(is.null),rows)
  if(!length(rows)) return(data.frame(status="UNAVAILABLE",reason="No applicable saved records"))
  columns <- unique(unlist(lapply(rows,names),use.names=FALSE))
  do.call(rbind,lapply(rows,function(x) {
    for(n in setdiff(columns,names(x))) x[[n]] <- rep(NA,nrow(x))
    x[columns]
  }))
}
p3_report_outputs <- function(out,run_id,candidates=NULL,final=NULL,descriptive=NULL) {
  for(n in c("tables","figures","diagnostics")) dir.create(file.path(out,n))
  write <- function(x,name,folder="tables") utils::write.csv(x,file.path(out,folder,paste0(name,".csv")),row.names=FALSE)
  figure <- function(name,fun) {
    grDevices::pdf(file.path(out,"figures",paste0(name,".pdf")),width=10,height=6)
    on.exit(grDevices::dev.off(),add=TRUE);fun()
  }
  artifacts <- Filter(Negate(is.null),list(candidates,final,descriptive));cfg <- artifacts[[1]]$config
  copy <- function(old,new,folder="tables") {
    path <- file.path(out,paste0(run_id,"-",old,".csv"))
    if(file.exists(path)) {
      file.copy(path,file.path(out,folder,paste0(new,".csv")),overwrite=FALSE)
    } else write(data.frame(status="UNAVAILABLE",reason=paste("No",old,"records")),new,folder)
  }
  copy("data-splits","split_definition");copy("yearly-peaks","seasonal_peak_metrics")
  copy("model-parameters","parameter_estimates");copy("fixed-value-provenance","fixed_value_provenance")
  audits <- Filter(Negate(is.null),lapply(artifacts,function(a) a$data_audit))
  write(if(length(audits)) audits[[1]] else data.frame(status="UNAVAILABLE",reason="Legacy artifact has no data audit"),"data_audit")
  maps <- Filter(Negate(is.null),lapply(artifacts,function(a)
    if(is.null(a$raw_mapping)) a$check$raw_mapping else a$raw_mapping))
  write(if(length(maps)) maps[[1]] else data.frame(status="UNAVAILABLE",reason="No raw/prepared mapping in stored artifact"),"prepared_row_mapping")
  selected <- if(!is.null(final)) final$lock$choice else if(!is.null(candidates)) candidates$selection else NULL
  if(is.null(selected)) {
    write(data.frame(status="DESCRIPTIVE_ONLY",reason="No temporal selection supplied"),"selected_model")
  } else {
    decision <- if(is.null(selected$candidate)) data.frame(status=selected$status) else selected$candidate
    decision$decision_status <- selected$status;decision$tie_rule <- cfg$selection$rule
    decision$tie_tolerance <- cfg$selection$tie_tolerance;decision$locked <- !is.null(final)
    decision$seed <- cfg$starts$seed
    decision$interpretation <- "B plug-in NB2 NLL only; penalized A objectives are not ranking scores"
    write(decision,"selected_model")
  }
  contexts <- list()
  if(!is.null(candidates)) {
    if(!is.null(candidates$baseline)) contexts[[length(contexts)+1L]] <- list(
      stage="A_to_B",id="M0",model="M0",K=NA_integer_,lambda=0,record=candidates$baseline)
    for(r in candidates$records) contexts[[length(contexts)+1L]] <- list(stage="A_to_B",
      id=r$candidate$id,model="M1",K=r$candidate$K,lambda=r$candidate$lambda,record=r)
  }
  peak_rows <- list()
  legacy_peaks <- file.path(out,paste0(run_id,"-yearly-peaks.csv"))
  if(file.exists(legacy_peaks)) peak_rows[[1L]] <- utils::read.csv(legacy_peaks,stringsAsFactors=FALSE)
  for(context in contexts) {
    r <- context$record
    if(!identical(r$status,"AVAILABLE")) next
    endpoint <- r$fit$representative;ctx <- endpoint$comparison_context;ids <- which(ctx$mask)
    for(stage in c("A","B")) {
      t <- if(stage=="B") r$forecast$peaks else p3_peaks(ctx$markers[ids,,drop=FALSE],
        ctx$cases[ids],endpoint$evaluation$states$expected_cases[ids],cfg$peaks)
      t$model <- context$model;t$candidate <- context$id;t$analysis <- stage
      peak_rows[[length(peak_rows)+1L]] <- t
    }
  }
  peak_table <- p3_bind_rows(peak_rows)
  write(peak_table,"seasonal_peak_metrics")
  figure("seasonal_peak_comparison",function() {
    if(!"eligible" %in% names(peak_table)) {
      graphics::plot.new();graphics::title("No available seasonal peaks");return(invisible(NULL))
    }
    p <- peak_table[!is.na(peak_table$eligible) & peak_table$eligible,,drop=FALSE]
    if(!nrow(p)) {graphics::plot.new();graphics::title("No complete eligible seasons");return(invisible(NULL))}
    if(!"candidate" %in% names(p)) p$candidate <- "selected"
    p$candidate[is.na(p$candidate)] <- "selected"
    groups <- interaction(p$model,p$analysis,p$candidate,drop=TRUE)
    graphics::plot(p$year,p$observed_peak_week,pch=16,ylim=c(1,52),xlab="Prepared year",ylab="Peak week",
      main=paste(cfg$data$region,"reported-case seasonal peaks"),sub="First chronological ties; incomplete-season exclusions in table")
    for(i in seq_along(levels(groups))) {
      d <- p[groups==levels(groups)[i],];graphics::points(d$year,d$expected_peak_week,col=i+1,pch=i+1)
    }
    graphics::legend("topright",c("observed",levels(groups)),col=c(1,seq_along(levels(groups))+1),
      pch=c(16,seq_along(levels(groups))+1),bty="n",cex=.6)
  })
  for(stage in c("final","descriptive")) {
    artifact <- if(stage=="final") final else descriptive
    if(is.null(artifact)) next
    for(n in names(artifact$models)) {
      r <- artifact$models[[n]];endpoint <- r$fit$representative
      contexts[[length(contexts)+1L]] <- list(stage=stage,id=n,model=n,
        K=if(is.null(endpoint)) NA_integer_ else endpoint$basis$K,
        lambda=if(is.null(endpoint)) NA_real_ else endpoint$lambda,record=r)
    }
  }
  endpoints <- list();coordinates <- list();candidate_rows <- list()
  value <- function(x) if(is.null(x) || !length(x)) NA_real_ else x
  text <- function(x) if(is.null(x)) "" else paste(x,collapse="; ")
  for(context in contexts) {
    record <- context$record;runs <- record$fit$runs
    if(!length(runs)) runs <- list(list(status=record$status,message=record$reason))
    for(i in seq_along(runs)) {
      r <- runs[[i]];e <- r$evaluation
      row <- data.frame(stage=context$stage,model=context$model,candidate=context$id,K=context$K,
        lambda=context$lambda,start=i,start_id=text(r$start_label),status=r$status,
        convergence=value(r$convergence),reason=text(c(r$message,e$failure_reason,record$reason)),
        objective=value(e$objective),loglik=value(e$loglik),penalty=value(e$penalty),
        B_mean_nb2_nll=if(context$stage=="A_to_B") value(record$forecast$scores$summary$mean_nb2_nll) else NA_real_,
        evaluation_class=if(context$stage=="A_to_B" && identical(record$status,"AVAILABLE")) {
          if(is.null(candidates$evaluation_policy$class)) "UNAVAILABLE" else candidates$evaluation_policy$class
        } else if(context$stage=="final" && identical(record$status,"AVAILABLE") &&
                  !is.null(record$evaluation_class)) record$evaluation_class else "UNAVAILABLE",
        boundary=text(r$boundary$coordinates),warnings=text(r$warnings),
        evaluations=value(r$failure_diagnostics$evaluations),invalid_evaluations=value(r$failure_diagnostics$invalid),
        failure_reasons=paste(capture.output(dput(r$failure_diagnostics$reasons)),collapse=" "))
      endpoints[[length(endpoints)+1L]] <- row
      if(context$stage=="A_to_B") candidate_rows[[length(candidate_rows)+1L]] <- row
      if(length(r$map$order)) {
        q <- function(x) if(is.null(x) || length(x)!=length(r$map$order)) rep(NA_real_,length(r$map$order)) else unname(x)
        coordinates[[length(coordinates)+1L]] <- data.frame(stage=context$stage,model=context$model,
          candidate=context$id,start=i,status=r$status,coordinate=r$map$order,
          initial=q(r$initial_coordinates),estimate=q(r$coordinates),lower=unname(r$map$lower),upper=unname(r$map$upper),
          guard_hit=r$map$order %in% r$boundary$coordinates)
      }
    }
  }
  write(p3_bind_rows(candidate_rows),"candidate_results");write(p3_bind_rows(coordinates),"parameter_diagnostics")
  diagnostic_table <- p3_bind_rows(endpoints);write(diagnostic_table,"optimizer_endpoints","diagnostics")
  figure("convergence_summary",function() {
    if(!"model" %in% names(diagnostic_table)) {
      graphics::plot.new();graphics::title("No saved optimizer endpoints");return(invisible(NULL))
    }
    counts <- table(diagnostic_table$status,diagnostic_table$stage)
    graphics::barplot(counts,beside=TRUE,legend.text=rownames(counts),
      ylab="Recorded starts",xlab="Fit stage",main=paste(cfg$data$region,"optimizer endpoint statuses"))
  })
  metrics <- list()
  if(!is.null(final)) for(n in names(final$models)) {
    r <- final$models[[n]]
    t <- if(identical(r$status,"AVAILABLE")) r$forecast$scores$summary else data.frame(reason=text(r$reason))
    t$model <- n;t$status <- r$status;t$evaluation_class <- if(is.null(r$evaluation_class)) "UNAVAILABLE" else r$evaluation_class
    t$score_definition <- "Plug-in NB2; no parameter uncertainty integration";metrics[[length(metrics)+1L]] <- t
  }
  write(p3_bind_rows(metrics),"test_metrics")
  ref <- if(is.null(descriptive)) NULL else descriptive$reference
  write(data.frame(category="Mref descriptive/reference only",status=if(is.null(ref)) "UNAVAILABLE" else ref$status,
    reason=if(is.null(ref)) "No descriptive artifact supplied" else ref$reason,
    diagnostic_status=if(is.null(ref$diagnostic_status)) "UNVERIFIED" else ref$diagnostic_status,
    forecast_eligible=FALSE,matched_primary_baseline=FALSE),"reference_comparability")
  if(!is.null(ref$audit)) write(ref$audit,"reference_audit")
  if(!is.null(ref$chain_diagnostics)) write(ref$chain_diagnostics,"reference_chains","diagnostics")
  if(!is.null(ref$expected_summary)) write(ref$expected_summary,"reference_expected_posterior_intervals")
  if(!is.null(ref$beta_summary)) write(ref$beta_summary,"reference_beta_posterior_intervals")
  policies <- Filter(Negate(is.null),list(candidates$evaluation_policy,
    final$models$M0$evaluation_policy,final$models$M1$evaluation_policy))
  write(p3_bind_rows(lapply(policies,function(x) data.frame(stage=x$stage,origin=x$origin,
    evaluation_class=x$class,definition=x$basis))),"evaluation_policy")
  write(p3_bind_rows(lapply(policies,function(x) {
    t <- x$information;t$stage <- x$stage;t$origin <- x$origin;t
  })),"information_timing")
  series_plot <- function(name,rows,label,boundary=NULL) {
    table <- p3_bind_rows(rows);write(table,paste0(name,"_series"))
    figure(name,function() {
      if(!"expected" %in% names(table) || !nrow(table)) {
        graphics::plot.new();graphics::title(main=paste(label,"UNAVAILABLE"));return(invisible(NULL))
      }
      first <- table[table$model==table$model[1],]
      graphics::plot(first$index,first$observed,pch=16,cex=.4,col="grey50",xlab="Prepared observation index",
        ylab="Reported cases/week",ylim=range(c(table$observed,table$expected)),
        main=paste(cfg$data$region,label),sub=paste(first$year[1],first$week[1],"through",tail(first$year,1),tail(first$week,1)))
      colors <- c(M0="blue",M1="red",Mref="darkgreen")
      for(n in unique(table$model)) {s <- table[table$model==n,];graphics::lines(s$index,s$expected,col=colors[[n]],lwd=1.5)}
      if(length(boundary)) graphics::abline(v=boundary+.5,lty=2)
      graphics::legend("topright",c("observed",unique(table$model)),col=c("grey50",colors[unique(table$model)]),
        pch=c(16,rep(NA,length(unique(table$model)))),lty=c(NA,rep(1,length(unique(table$model)))),bty="n")
    })
  }
  rows_for <- function(record,model,stage,marks,observed,ids) {
    if(!identical(record$status,"AVAILABLE")) return(NULL)
    expected <- if(stage=="A") record$fit$representative$evaluation$states$expected_cases else record$forecast$states$expected_cases
    data.frame(index=ids,marks[ids,,drop=FALSE],observed=observed,expected=expected,model=model,check.names=FALSE)
  }
  if(!is.null(candidates)) {
    records <- Filter(function(r) identical(r$candidate$id,candidates$selection$id),candidates$records)
    records <- list(M0=candidates$baseline,M1=if(length(records)) records[[1]] else NULL)
    for(stage in c("A","B")) {
      ids <- candidates$check$splits$indices[[stage]]
      rows <- lapply(names(records),function(n) {
        r <- records[[n]];if(is.null(r)) return(NULL)
        observed <- if(stage=="A") r$fit$representative$comparison_context$cases else r$forecast$scores$weekly$observed
        if(is.null(observed)) return(NULL)
        rows_for(r,n,stage,candidates$check$markers,observed,ids)
      })
      label <- if(stage=="A") "A fitted (initial row unscored)" else paste("B prediction",candidates$evaluation_policy$class)
      series_plot(paste0(stage,"_cases"),rows,label,max(candidates$check$splits$indices$A))
    }
  }
  if(!is.null(final)) {
    ids <- final$splits$indices$C
    rows <- lapply(names(final$models),function(n) rows_for(final$models[[n]],n,"C",final$markers,final$observed[ids],ids))
    series_plot("C_cases",rows,paste("C future evaluation",final$evaluation_class),max(final$splits$indices$B))
  }
  for(name in c("candidate-B-scores","selected-beta-shape","yearly-peaks")) {
    file <- file.path(out,paste0(run_id,"-",name,".pdf"))
    if(file.exists(file)) file.copy(file,file.path(out,"figures",paste0(name,".pdf")))
  }
  invisible(NULL)
}
