phase28_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)<3L || !args[1L] %in% c("prepare","verify","fit","compare") ||
     (length(args)-1L)%%2L)
    stop("Usage: Rscript scripts/phase28.R STAGE --out NEW_DIR [--prepare PATH] [--verify PATH] [--model M2|PARK] [--m2 PATH] [--park PATH]")
  root <- normalizePath(file.path(dirname(sub("^--file=","",
    grep("^--file=",commandArgs(),value=TRUE)[1L])),".."))
  setwd(root)
  sys.source("R/load_phase6.R",globalenv());load_spline_sirs_phase6(root,globalenv())
  for(f in c("R/phase5_reporting.R","R/phase10_peak_identity.R","R/result_storage.R",
    "R/accepted_extension.R","R/phase21.R","R/phase22.R","R/phase22_analysis.R",
    "R/phase24.R","R/phase25.R","R/phase26.R","R/phase27.R","R/phase28.R"))
    sys.source(f,globalenv())
  for(f in c("config/phase22.R","config/phase24.R","config/phase16.R",
    "config/phase25.R","config/phase26.R","config/phase27.R","config/phase28.R"))
    sys.source(f,globalenv())
  c22 <- phase22_config();c26 <- phase26_config();c27 <- phase27_config()
  cfg <- phase28_config()
  rest <- args[-1L]
  keys <- if(length(rest)) rest[seq(1L,length(rest),by=2L)] else character()
  vals <- if(length(rest)) rest[seq(2L,length(rest),by=2L)] else character()
  if(anyDuplicated(keys) || any(!keys %in% c("--out","--prepare","--verify","--model","--m2","--park")))
    stop("Invalid options")
  opt <- as.list(setNames(vals,substring(keys,3L)))
  if(is.null(opt$out) || !nzchar(opt$out)) stop("Missing --out NEW_DIR")
  stage <- args[1L]
  if(stage=="fit" && (is.null(opt$model) || !opt$model %in% c("M2","PARK")))
    stop("fit requires --model M2 or PARK")
  out <- reserve_phase2_output(opt$out)
  log_path <- file.path(out,"execution.log")
  log_event <- function(event,detail="") {
    line <- paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S %z"),
      stage,if(stage=="fit") opt$model else "",event,detail,sep=" | ")
    cat(line,"\n",file=log_path,append=TRUE,sep="")
    message(line)
  }
  writeLines(c("RUNNING",stage),file.path(out,"STATUS"))
  start <- Sys.time()
  log_event("START",paste("output",out))
  read_stage <- function(path,want) {
    if(is.null(path)) stop("Missing --",want)
    x <- read_phase2_result(path)
    if(!identical(x$metadata$phase,28L) || !identical(x$metadata$stage,want))
      stop("Wrong Phase 28 stage: ",path)
    x$result
  }
  tryCatch({
    log_event("PREPARATION_START")
    inputs <- p22_inputs(c22)
    ref <- p27_reference(inputs,c22,c26,c27)
    panel <- p28_panel(inputs,ref,cfg)
    log_event("PREPARATION_COMPLETED",
      sprintf("rows=%d pre_intervention=%d",panel$N,panel$Nnonpi))
    if(stage=="prepare") {
      log_event("ALIGNMENT_EXPORT_START")
      p28_write(panel$alignment,out,"data_alignment")
      p28_write(panel$exclusions,out,"data_exclusions")
      result <- list(status="PREPARED",N=panel$N,Nnonpi=panel$Nnonpi,
        source_rows=panel$source_rows)
      log_event("ALIGNMENT_EXPORT_COMPLETED")
    } else if(stage=="verify") {
      log_event("VERIFICATION_START")
      prepared <- read_stage(opt$prepare,"prepare")
      if(prepared$N!=panel$N || prepared$Nnonpi!=panel$Nnonpi)
        stop("Prepared panel dimensions changed")
      p28_write(ref$checks,out,"m2_replay_checks")
      equations <- p28_verify_equations(panel,ref,inputs)
      p28_write(equations,out,"model_equation_checks")
      if(!all(ref$checks$passed) || !all(equations$passed))
        stop("M2 replay or common equations failed")
      result <- list(status="VERIFIED",checks=ref$checks,equations=equations)
      log_event("VERIFICATION_COMPLETED")
    } else if(stage=="fit") {
      log_event("FIT_INPUT_CHECK_START")
      prepared <- read_stage(opt$prepare,"prepare")
      if(prepared$N!=panel$N || prepared$Nnonpi!=panel$Nnonpi)
        stop("Prepared panel dimensions changed")
      read_stage(opt$verify,"verify")
      prep_end <- Sys.time()
      log_event("FIT_INPUT_CHECK_COMPLETED")
      if(opt$model=="M2") {
        penalty <- ref$fit$penalty
        log_event("OPTIMIZATION_START","three Phase 26 starts")
        fit <- p26_fit(inputs,c22,c26,"M2",penalty,out)
        fit_end <- Sys.time()
        log_event("OPTIMIZATION_COMPLETED",
          sprintf("selected_start=%d objective=%.9f",fit$best,fit$score$objective))
        log_event("POSTPROCESSING_START")
        score <- p26_score(fit$endpoint,inputs,p22_problem(inputs),c26,"M2",penalty)
        saveRDS(fit,file.path(out,"fit.rds"))
        p28_write(fit$parameter_table,out,"parameter_inventory_detail")
        p28_write(fit$numerical,out,"numerical_diagnostics")
        p28_write(score$weekly,out,"weekly_fit")
        result <- list(status="FIT_SAVED",model="M2",endpoint=fit$endpoint,
          penalty=penalty,objective=score$objective,
          fit_seconds=as.numeric(difftime(fit_end,prep_end,units="secs")))
        compile_seconds <- 0
        log_event("POSTPROCESSING_COMPLETED")
      } else {
        if(!requireNamespace("rstan",quietly=TRUE)) stop("PARK fit requires rstan")
        log_event("COMPILATION_START",cfg$park_stan)
        compiled <- rstan::stan_model(cfg$park_stan)
        compile_end <- Sys.time()
        log_event("COMPILATION_COMPLETED",
          sprintf("seconds=%.1f",as.numeric(difftime(compile_end,prep_end,units="secs"))))
        standata <- p28_standata(panel)
        init <- lapply(seq_len(cfg$park_chains),function(chain) {
          set.seed(cfg$park_seed+chain)
          list(S0=.08,I0=.0005,rho=.02,
            beta=stats::rnorm(52L,8,.01),sigma=.5,phi=10,
            tau=200,omega=100,
            npi=rep(1,standata$N-standata$Nnonpi),sigma_npi=.05)
        })
        log_event("SAMPLING_START",sprintf("chains=%d cores=%d iterations=%d",
          cfg$park_chains,cfg$park_cores,cfg$park_iter))
        fit <- rstan::sampling(compiled,data=standata,seed=cfg$park_seed,
          chains=cfg$park_chains,cores=cfg$park_cores,iter=cfg$park_iter,
          init=init,
          control=list(max_treedepth=cfg$park_treedepth,
            adapt_delta=cfg$park_adapt_delta))
        fit_end <- Sys.time()
        log_event("SAMPLING_COMPLETED",
          sprintf("seconds=%.1f",as.numeric(difftime(fit_end,compile_end,units="secs"))))
        log_event("POSTPROCESSING_START")
        draws <- rstan::extract(fit,pars=c("C","beta","npi","S","I","R","rho","phi"),
          permuted=TRUE,inc_warmup=FALSE)
        expected <- draws$C
        if(ncol(expected)!=panel$N || any(!is.finite(expected)))
          stop("PARK posterior expected-case draws invalid")
        curves <- data.frame(panel$clock,observed=panel$observed,
          expected_mean=colMeans(expected),
          expected_median=apply(expected,2L,stats::median))
        p28_write(curves,out,"weekly_fit")
        saveRDS(fit,file.path(out,"fit.rds"))
        saveRDS(draws,file.path(out,"posterior_draws.rds"))
        if(!requireNamespace("posterior",quietly=TRUE))
          stop("PARK diagnostics require posterior")
        pars <- c("S0","I0","rho","beta","sigma","phi","tau","omega",
          "npi","sigma_npi")
        diagnostic_draws <- posterior::as_draws_array(as.array(fit,pars=pars))
        diagnostics <- posterior::summarise_draws(diagnostic_draws,
          rhat=posterior::rhat,ess_bulk=posterior::ess_bulk,
          ess_tail=posterior::ess_tail)
        p28_write(diagnostics,out,"numerical_diagnostics")
        sampler <- rstan::get_sampler_params(fit,inc_warmup=FALSE)
        hmc <- do.call(rbind,lapply(seq_along(sampler),function(i) {
          x <- sampler[[i]]
          energy <- x[,"energy__"]
          data.frame(chain=i,divergences=sum(x[,"divergent__"]),
            treedepth_hits=sum(x[,"treedepth__"]>=cfg$park_treedepth),
            bfmi=mean(diff(energy)^2)/stats::var(energy))
        }))
        p28_write(hmc,out,"hmc_diagnostics")
        if(!all(c("rhat","ess_bulk","ess_tail") %in% names(diagnostics)) ||
           !nrow(diagnostics)) stop("PARK diagnostics columns missing")
        diagnostic_pass <- all(is.finite(diagnostics$rhat) &
          diagnostics$rhat<=cfg$rhat_max &
          diagnostics$ess_bulk>=cfg$ess_min &
          diagnostics$ess_tail>=cfg$ess_min) &&
          all(hmc$divergences==0L & hmc$treedepth_hits==0L &
            hmc$bfmi>=cfg$bfmi_min)
        result <- list(status="FIT_SAVED",model="PARK",draws=nrow(expected),
          fit_seconds=as.numeric(difftime(fit_end,compile_end,units="secs")),
          diagnostic_pass=diagnostic_pass)
        compile_seconds <- as.numeric(difftime(compile_end,prep_end,units="secs"))
        log_event("POSTPROCESSING_COMPLETED",
          paste("diagnostic_pass",isTRUE(diagnostic_pass)))
      }
      timing <- p28_timing_row("fit",opt$model,start,
        as.numeric(difftime(prep_end,start,units="secs")),compile_seconds,
        result$fit_seconds,as.numeric(difftime(Sys.time(),fit_end,units="secs")))
      p28_write(timing,out,"timing_metrics")
      log_event("TIMING_SAVED")
    } else {
      log_event("COMPARISON_START")
      read_stage(opt$prepare,"prepare")
      read_stage(opt$verify,"verify")
      m2 <- read_stage(opt$m2,"fit");park <- read_stage(opt$park,"fit")
      if(m2$model!="M2" || park$model!="PARK") stop("Model arm mismatch")
      park_diag <- utils::read.csv(file.path(dirname(opt$park),
        "numerical_diagnostics.csv"),check.names=FALSE)
      names(park_diag) <- sub("^posterior::","",names(park_diag))
      park_hmc <- utils::read.csv(file.path(dirname(opt$park),"hmc_diagnostics.csv"))
      required_diag <- c("rhat","ess_bulk","ess_tail")
      if(!all(required_diag %in% names(park_diag)) || !nrow(park_diag) ||
         !nrow(park_hmc)) stop("PARK diagnostics missing")
      park$diagnostic_pass <- isTRUE(all(is.finite(park_diag$rhat) &
        park_diag$rhat<=cfg$rhat_max & park_diag$ess_bulk>=cfg$ess_min &
        park_diag$ess_tail>=cfg$ess_min)) &&
        isTRUE(all(park_hmc$divergences==0L &
          park_hmc$treedepth_hits==0L & park_hmc$bfmi>=cfg$bfmi_min))
      p28_write(park_diag,out,"numerical_diagnostics")
      mz <- utils::read.csv(file.path(dirname(opt$m2),"weekly_fit.csv"))
      pz <- utils::read.csv(file.path(dirname(opt$park),"weekly_fit.csv"))
      if(nrow(mz)!=panel$N || nrow(pz)!=panel$N ||
         any(mz$observed!=panel$observed | pz$observed!=panel$observed))
        stop("Fit panels differ")
      curves <- list(M2=mz$expected,PARK_median=pz$expected_median,
        PARK_mean=pz$expected_mean)
      metrics <- do.call(rbind,lapply(names(curves),function(k)
        p28_curve_metrics(panel$observed,curves[[k]],panel$mask,
          panel$clock,if(k=="M2") "M2" else "PARK",k)))
      p28_write(metrics,out,"fit_metrics")
      timing <- do.call(rbind,lapply(names(curves),function(k)
        p28_timing(panel$observed,curves[[k]],panel$clock,
          if(k=="M2") "M2" else "PARK",k)))
      p28_write(timing,out,"peak_metrics")
      movement <- do.call(rbind,lapply(split(timing,
        interaction(timing$model,timing$summary,drop=TRUE)),function(x) {
        a <- x[x$season %in% c(2015L,2017L),]
        a <- a[order(a$season),]
        if(nrow(a)!=2L || !all(a$complete)) return(NULL)
        peak_observed <- diff(a$observed_peak)
        peak_fitted <- diff(a$fitted_peak)
        cog_observed <- diff(a$observed_cog)
        cog_fitted <- diff(a$fitted_cog)
        data.frame(model=a$model[1L],summary=a$summary[1L],
          observed_peak_movement=peak_observed,fitted_peak_movement=peak_fitted,
          peak_movement_error=peak_fitted-peak_observed,
          observed_cog_movement=cog_observed,fitted_cog_movement=cog_fitted,
          cog_movement_error=cog_fitted-cog_observed,
          endpoint_peak_pass=all(abs(a$peak_error)<=c26$endpoint_peak_tolerance),
          peak_movement_pass=abs(peak_fitted-peak_observed)<=c26$peak_change_tolerance,
          cog_movement_pass=abs(cog_fitted-cog_observed)<=c26$cog_change_tolerance)
      }))
      p28_write(movement,out,"movement_metrics")
      inventory <- p28_inventory(ref,panel)
      p28_write(inventory,out,"parameter_inventory")
      runtime <- rbind(utils::read.csv(file.path(dirname(opt$m2),"timing_metrics.csv")),
        utils::read.csv(file.path(dirname(opt$park),"timing_metrics.csv")))
      p28_write(runtime,out,"runtime_metrics")
      comparison <- data.frame(dimension=c("shared states","seasonality",
        "changes across years","separable components","confounding",
        "uncertainty"),
        M2=c("S, I, R and expected reported cases",
          "shared seasonal spline","historical level and shape; later spline",
          "shared shape, historical level, historical shape and later correction",
          "transmission remains confounded with reporting and initial states",
          "point estimate and numerical diagnostics"),
        PARK=c("S, I, R and expected reported cases",
          "52 weekly beta coordinates with cyclic random-walk prior",
          "later weekly npi multiplier","shared beta and later multiplier",
          "transmission remains confounded with reporting and initial states",
          "joint posterior draws and HMC diagnostics"))
      p28_write(comparison,out,"interpretability_comparison")
      grDevices::pdf(file.path(out,"observed_fitted.pdf"),width=10,height=7)
      graphics::plot(panel$observed,type="l",col="black",xlab="Common panel row",
        ylab="Reported cases")
      graphics::lines(mz$expected,col="steelblue")
      graphics::lines(pz$expected_median,col="firebrick")
      detail <- which(panel$clock$year>=2015L & panel$clock$year<=2017L)
      graphics::plot(detail,panel$observed[detail],type="l",col="black",
        xlab="Common panel row",ylab="Reported cases")
      graphics::lines(detail,mz$expected[detail],col="steelblue")
      graphics::lines(detail,pz$expected_median[detail],col="firebrick")
      grDevices::dev.off()
      grDevices::pdf(file.path(out,"peak_cog_errors.pdf"),width=8,height=5)
      complete <- timing$complete & timing$season %in% 2015:2017
      graphics::plot(timing$peak_error[complete],timing$cog_error[complete],
        xlab="Peak-week error",ylab="Center-of-gravity error")
      grDevices::dev.off()
      grDevices::pdf(file.path(out,"accuracy_runtime_coordinates.pdf"),width=9,height=6)
      graphics::par(mfrow=c(1,3))
      full <- metrics[metrics$period=="full" &
        metrics$summary %in% c("M2","PARK_median"),]
      graphics::barplot(full$mae,names.arg=full$summary,las=2,main="MAE")
      graphics::barplot(runtime$total_seconds,names.arg=runtime$model,
        main="Wall seconds")
      graphics::barplot(inventory$free_coordinates,names.arg=inventory$model,
        main="Free coordinates")
      grDevices::dev.off()
      park_draws <- readRDS(file.path(dirname(opt$park),"posterior_draws.rds"))
      park_beta <- colMeans(park_draws$beta)
      park_npi <- c(rep(1,panel$Nnonpi),colMeans(park_draws$npi))
      m2_beta <- mz$seasonal_beta
      grDevices::pdf(file.path(out,"transmission_components.pdf"),width=10,height=8)
      graphics::par(mfrow=c(2,3))
      graphics::plot(1:52,m2_beta[match(1:52,mz$seasonal_week)],type="l",
        main="M2 common seasonal beta",xlab="Week",ylab="beta")
      graphics::plot(mz$log_level,type="l",main="M2 historical log level",
        xlab="Row",ylab="log multiplier")
      graphics::lines(mz$log_shape,col="steelblue")
      graphics::lines(log(mz$future_q),col="firebrick")
      graphics::plot(mz$effective_beta,type="l",main="M2 effective beta",
        xlab="Row",ylab="beta")
      graphics::plot(1:52,park_beta,type="l",main="PARK seasonal beta",
        xlab="Week",ylab="beta")
      graphics::plot(park_npi,type="l",main="PARK later multiplier",
        xlab="Row",ylab="multiplier")
      graphics::plot(park_beta[panel$clock$seasonal_week]*park_npi,
        type="l",main="PARK effective beta",xlab="Row",ylab="beta")
      grDevices::dev.off()
      p28_write(data.frame(panel$clock,park_seasonal_beta=park_beta[panel$clock$seasonal_week],
        park_npi=park_npi,park_effective_beta=park_beta[panel$clock$seasonal_week]*park_npi),
        out,"park_transmission_components")
      pick <- function(model,summary) metrics[metrics$period=="full" &
        metrics$model==model & metrics$summary==summary,]
      mm <- pick("M2","M2");pm <- pick("PARK","PARK_median")
      pmean <- pick("PARK","PARK_mean")
      writeLines(c("# Phase 28 focused comparison","",
        "Comparator: original Park main Honshu model, not its extensions.",
        "All accuracy numbers are retrospective on the common panel. Runtime is one timed fit per model.","",
        "## Fitting accuracy","",
        sprintf("Full-period MAE: M2 %.4f; PARK posterior median %.4f.",mm$mae,pm$mae),
        sprintf("Full-period RMSE: M2 %.4f; PARK posterior mean %.4f.",mm$rmse,pmean$rmse),
        "See fit_metrics.csv, peak_metrics.csv, and movement_metrics.csv for periods and endpoint errors.","",
        "## Computational cost","",
        paste("Measured total wall seconds:",paste(runtime$model,
          round(runtime$total_seconds,1),collapse="; ")),
        "The M2 timing depends on earlier fitted nuisance coordinates and produces a point estimate; PARK also produces posterior uncertainty.","",
        "## Model complexity","",
        paste("Independent coordinates:",paste(inventory$model,
          inventory$free_coordinates,collapse="; ")),
        "These are counts of free coordinates, not effective degrees of freedom.","",
        "## Interpretability","",
        "M2 explicitly separates shared seasonal shape, historical level and shape changes, and a later spline correction. PARK separates a shared weekly seasonal beta from later weekly multipliers. Neither decomposition establishes causal identification.","",
        sprintf("PARK diagnostic gate: %s.",if(isTRUE(park$diagnostic_pass)) "passed" else "failed or unavailable"),
        "Do not infer forecasting performance or statistical significance from these fitted differences."),
        file.path(out,"phase28_report.md"))
      result <- list(status="COMPARISON_SAVED",metrics=metrics,
        inventory=inventory,runtime=runtime)
      log_event("COMPARISON_COMPLETED")
    }
    log_event("PROVENANCE_START")
    sources <- c(cfg$park_stan,cfg$park_driver,cfg$park_preprocessing,
      cfg$park_license,cfg$original_cases,cfg$original_population,
      cfg$m2_fit,cfg$m2_evaluation,unlist(opt[c("prepare","verify","m2","park")]))
    metadata <- phase2_metadata(root,"config/phase28.R",sources)
    metadata$phase <- 28L;metadata$stage <- stage
    metadata$park_commit <- cfg$park_commit;metadata$park_license <- "MIT"
    metadata$config <- cfg
    dput(list(stage=stage,code_commit=metadata$code_commit,
      code_hashes=metadata$code_hashes,input_hashes=metadata$input_hashes,
      park_commit=cfg$park_commit,license="MIT"),file=file.path(out,"manifest.R"))
    write_phase2_result(result,out,metadata)
    log_event("STORAGE_COMPLETED")
  },error=function(e) {
    log_event("FAILED",conditionMessage(e))
    if(!identical(readLines(file.path(out,"STATUS"),n=1L),"STORAGE_FAILED"))
      writeLines(c("RUN_FAILED",conditionMessage(e)),file.path(out,"STATUS"))
    stop(e)
  })
}
if(sys.nframe()==0L) phase28_main()
