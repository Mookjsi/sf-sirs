phase26_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  stages <- c("replay","m1-fit","m1-evaluate","m2-fit","m2-evaluate","compare","stress","audit")
  if(length(args)<3L || !args[1L] %in% stages || args[2L]!="--out")
    stop("Usage: Rscript scripts/phase26.R STAGE --out NEW_DIR [--replay PATH] [--fit PATH] [--m1 PATH] [--m2 PATH] [--penalty primary|light|strong]")
  root <- normalizePath(file.path(dirname(sub("^--file=","",
    grep("^--file=",commandArgs(),value=TRUE)[1L])),".."))
  setwd(root)
  sys.source("R/load_phase6.R",globalenv());load_spline_sirs_phase6(root,globalenv())
  for(f in c("R/phase5_reporting.R","R/phase10_peak_identity.R","R/result_storage.R",
             "R/accepted_extension.R","R/phase21.R","R/phase22.R","R/phase22_analysis.R",
             "R/phase24.R","R/phase26.R")) sys.source(f,globalenv())
  source("config/phase22.R",local=TRUE);source("config/phase26.R",local=TRUE)
  cfg <- phase26_config();source_cfg <- phase22_config()
  kv <- args[-(1:3)]
  option_keys <- kv[2L*seq_len(length(kv)%/%2L)-1L]
  if(length(kv)%%2L || any(!option_keys %in%
    c("--replay","--fit","--m1","--m2","--penalty"))) stop("Invalid stage options")
  opt <- if(length(kv)) setNames(as.list(kv[2L*seq_len(length(kv)%/%2L)]),sub("^--","",option_keys)) else list()
  stage <- args[1L];out <- reserve_phase2_output(args[3L])
  writeLines(c("RUNNING",stage),file.path(out,"STATUS"))
  inputs <- p22_inputs(source_cfg)
  paths <- c(unlist(source_cfg[c("reference","replay","panel","plan","fit","evaluation","joint_fit")]),
    cfg$phase24,cfg$phase25,cfg$phase24_weekly)
  read_stage <- function(path,want) {
    if(is.null(path) || !file.exists(path)) stop("Missing ",want," artifact")
    x <- read_phase2_result(path)
    if(!identical(x$metadata$phase,26L) || !identical(x$metadata$stage,want) ||
       !identical(x$result$config_version,cfg$version)) stop("Wrong Phase 26 stage or version: ",path)
    if(anyNA(tools::md5sum(x$metadata$input_paths)) ||
       !identical(unname(tools::md5sum(x$metadata$input_paths)),unname(x$metadata$input_hashes)))
      stop("Stage provenance changed: ",path)
    paths <<- c(paths,path)
    x$result
  }
  result <- tryCatch({
    if(stage=="replay") {
      z <- p26_replay(inputs,source_cfg,cfg)
      utils::write.csv(z$checks,file.path(out,"baseline_replay.csv"),row.names=FALSE)
      if(!z$passed) stop("Baseline replay gate failed: ",
        paste(z$checks$check[!z$checks$passed],collapse=", "))
      list(status="REPLAY_PASSED",config_version=cfg$version,checks=z$checks,
        source_hashes=z$source_hashes)
    } else {
      read_stage(opt[["replay"]],"replay")
      penalty_id <- if(is.null(opt[["penalty"]])) "primary" else opt[["penalty"]]
      penalty <- cfg$penalty_grid[cfg$penalty_grid$id==penalty_id,]
      if(nrow(penalty)!=1L) stop("Unknown preregistered penalty")
      model <- if(grepl("m2",stage)) "M2" else "M1"
      if(stage=="audit") {
        fit_path <- opt[["fit"]]
        if(is.null(fit_path)) stop("Audit requires --fit exact saved candidate fit")
        x <- read_phase2_result(fit_path)
        if(!identical(x$metadata$phase,26L) || !x$metadata$stage %in% c("m1-fit","m2-fit") ||
           !identical(x$result$config_version,cfg$version)) stop("Invalid saved candidate")
        paths <- c(paths,fit_path)
        z <- p26_candidate_audits(x$result,inputs,cfg,out)
        z$config_version <- cfg$version
        z
      } else if(stage=="m1-fit" || stage=="m2-fit") {
        if(model=="M2") {
          m1 <- read_stage(opt[["m1"]],"m1-evaluate")
          if(isTRUE(m1$transition$waveform_adequate)) stop("M2 condition not met: M1 waveform adequate")
          baseline_transition <- p26_transition(m1$baseline_metrics,cfg)
          if(abs(m1$transition$peak_change_error)>=abs(baseline_transition$peak_change_error))
            stop("M2 condition not met: M1 did not improve peak-change timing")
        }
        fit <- p26_fit(inputs,source_cfg,cfg,model,penalty,out)
        fit$config_version <- cfg$version
        fit
      } else if(stage=="m1-evaluate" || stage=="m2-evaluate") {
        fit <- read_stage(opt[["fit"]],paste0(tolower(model),"-fit"))
        ev <- p26_evaluate(fit,inputs,source_cfg,cfg,out)
        ev$config_version <- cfg$version
        ev
      } else if(stage=="stress") {
        clock <- p26_clock(inputs)
        mask <- c(FALSE,clock$year[-1L]<=2017L)
        fit <- p26_fit(inputs,source_cfg,cfg,"M1",penalty,out,mask,
          projected_seasons=2018:2019,holdout=TRUE)
        score <- p26_score(fit$endpoint,inputs,p22_problem(inputs),cfg,"M1",penalty,
          projected_seasons=2018:2019)
        z <- score$weekly
        h <- z$year %in% 2018:2019
        holdout <- data.frame(period="2018-2019",projection="zero",
          rows=sum(h),nb2_loglik=sum(z$nb2_loglik[h]),
          mae=mean(abs(z$observed[h]-z$expected[h])),
          rmse=sqrt(mean((z$observed[h]-z$expected[h])^2)))
        utils::write.csv(holdout,file.path(out,"holdout_metrics.csv"),row.names=FALSE)
        utils::write.csv(z,file.path(out,"weekly_fit.csv"),row.names=FALSE)
        fit$holdout_metrics <- holdout;fit$config_version <- cfg$version;fit
      } else {
        m1 <- read_stage(opt[["m1"]],"m1-evaluate")
        m2 <- if(!is.null(opt[["m2"]])) read_stage(opt[["m2"]],"m2-evaluate") else NULL
        baseline_transition <- p26_transition(m1$baseline_metrics,cfg)
        if(is.null(m2) && !isTRUE(m1$transition$waveform_adequate) &&
           abs(m1$transition$peak_change_error)<abs(baseline_transition$peak_change_error))
          stop("M2 condition met; fit and evaluate M2 before final comparison")
        models <- c("M0","M1",if(!is.null(m2)) "M2")
        scores <- list(m1$baseline,m1$score,if(!is.null(m2)) m2$score)
        fits <- list(NULL,m1,m2)
        comparison <- do.call(rbind,lapply(seq_along(models),function(i) {
          s <- scores[[i]];p <- if(i==1L) 0 else s$level_penalty+s$shape_penalty
          data.frame(model=models[i],loglik=s$loglik,
            penalized_objective=s$objective,added_parameters=if(i==1L) 0L else
              length(cfg$historical_seasons)-1L+if(i==3L) 2L*length(cfg$historical_seasons) else 0L,
            historical_loglik=s$historical_loglik,future_loglik=s$future_loglik,
            correction_penalty=p)
        }))
        utils::write.csv(comparison,file.path(out,"model_comparison.csv"),row.names=FALSE)
        periods <- do.call(rbind,lapply(seq_along(models),function(i)
          p26_periods(scores[[i]],models[i])))
        utils::write.csv(periods,file.path(out,"period_fit_comparison.csv"),row.names=FALSE)
        metrics <- do.call(rbind,lapply(seq_along(models),function(i)
          p26_season_metrics(scores[[i]],inputs,cfg,models[i])))
        utils::write.csv(metrics,file.path(out,"season_metrics.csv"),row.names=FALSE)
        years <- do.call(rbind,lapply(seq_along(models),function(i)
          p26_year_metrics(scores[[i]],models[i])))
        utils::write.csv(years,file.path(out,"year_metrics.csv"),row.names=FALSE)
        transitions <- do.call(rbind,lapply(seq_along(models),function(i)
          p26_transition(metrics[metrics$model==models[i],],cfg)))
        utils::write.csv(transitions,file.path(out,"transition_metrics.csv"),row.names=FALSE)
        selected <- if(!is.null(m2) && isTRUE(m2$transition$waveform_adequate)) "M2" else
          if(isTRUE(m1$transition$waveform_adequate)) "M1" else "NONE"
        criteria <- rbind(transform(m1$transition,model="M1"),
          if(!is.null(m2)) transform(m2$transition,model="M2"))
        report <- c("# Phase 26 empirical comparison","",
          paste("Selected candidate:",selected),"",
          "## Prespecified criteria", "",
          capture.output(print(criteria,row.names=FALSE)),"",
          "## Model scores", "",
          capture.output(print(comparison,row.names=FALSE)),"",
          "## Period scores", "",
          capture.output(print(periods,row.names=FALSE)),"",
          if(selected=="NONE") "Neither candidate met the prespecified waveform standard." else
            paste("The",selected,"candidate met the prespecified waveform standard."),
          "Selection requires timing, center of gravity and weekly fit criteria; likelihood alone does not select a model.",
          "Year-specific corrections are latent effective transmission, not identified NPI or causal effects.",
          "Changes in initial states, reporting and spline shape remain potential confounders.",
          "The through-2017 stress fit is a separate generalization check; it does not test pandemic forecasting.",
          "Phase 25 mechanism claims remain conditional on their original JOINT reference.",
          "The mechanism, S-versus-beta and annual-versus-semiannual audits must be rerun under the exact selected endpoint before transferring conclusions.",
          "", "## Input hashes", "",
          paste(names(tools::md5sum(paths)),tools::md5sum(paths),sep=": "))
        writeLines(report,file.path(out,"phase26_report.md"))
        list(status="COMPARED",config_version=cfg$version,selected=selected,
          comparison=comparison,periods=periods,metrics=metrics,years=years,
          transitions=transitions)
      }
    }
  },error=function(e) {
    writeLines(c("FAILED",conditionMessage(e)),file.path(out,"STATUS"))
    stop(e)
  })
  metadata <- phase2_metadata(root,"config/phase26.R",paths)
  metadata$phase <- 26L;metadata$stage <- stage;metadata$config <- cfg
  write_phase2_result(result,out,metadata)
  invisible(result)
}
if(sys.nframe()==0L) phase26_main()
