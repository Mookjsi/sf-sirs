# Read-only adapter for the inspected reproduction pipeline's exported medians.
# No Stan fit is loaded and no Bayesian engine or diagnostics are executed.
p3_reference <- function(spec,data) {
  diagnostic_status <- "UNVERIFIED";quality <- NULL;chains <- NULL
  unavailable <- function(reason,audit=NULL,paths=character()) list(status="UNAVAILABLE",reason=reason,
    audit=audit,paths=paths,forecast_eligible=FALSE,diagnostic_status=diagnostic_status,
    quality=quality,chain_diagnostics=chains)
  if(is.null(spec$attempt_dir) || is.null(spec$data_dir)) return(unavailable("Mref paths not configured"))
  paths <- c(manifest=file.path(spec$attempt_dir,"manifest.rds"),
    standata=file.path(spec$attempt_dir,"standata.rds"),
    C=file.path(spec$attempt_dir,"extract/C_summary.csv"),
    beta=file.path(spec$attempt_dir,"extract/beta_summary.csv"),
    sampling_status=file.path(spec$attempt_dir,"COMPLETE"),
    extraction_status=file.path(spec$attempt_dir,"extract/STATUS"),
    cases=file.path(spec$data_dir,"japan_rsv_island.csv"),
    population=file.path(spec$data_dir,"japan_pop_island.csv"),
    quality=file.path(spec$attempt_dir,"diagnostics/quality_checks.rds"),
    chains=file.path(spec$attempt_dir,"diagnostics/chain_diagnostics.csv"))
  if(any(!file.exists(paths))) return(unavailable(paste("Missing Mref files:",
    paste(names(paths)[!file.exists(paths)],collapse=", ")),paths=paths[file.exists(paths)]))
  tryCatch({
    if(!identical(readLines(paths[["sampling_status"]],n=1L),"SAMPLING_COMPLETED_NOT_CONVERGENCE_CERTIFIED") ||
       !identical(readLines(paths[["extraction_status"]],n=1L),"COMPLETED")) stop("Mref stage completion markers unsupported or absent")
    m <- readRDS(paths[["manifest"]]);s <- readRDS(paths[["standata"]])
    quality <- readRDS(paths[["quality"]]);chains <- utils::read.csv(paths[["chains"]])
    flags <- c("draws_finite","dimensions_ok","at_least_four_chains","rank_rhat_pass",
      "ess_bulk_pass","ess_tail_pass","divergences_pass","treedepth_pass","bfmi_pass")
    if(!all(flags %in% names(quality)) ||
       !all(vapply(quality[flags],function(x) is.logical(x) && length(x)==1L && !is.na(x),logical(1))))
      stop("Unsupported Mref diagnostic structure")
    diagnostic_status <- if(all(unlist(quality[flags]))) "PASS_DECLARED_DIAGNOSTICS" else "REVIEW_REQUIRED"
    if(!identical(quality$overall,diagnostic_status)) stop("Inconsistent Mref diagnostic status")
    if(!all(c("chain","draws","divergences","depth_hits","ebfmi") %in% names(chains)) || anyNA(chains))
      stop("Unsupported Mref chain diagnostics")
    md <- spec$metadata
    required <- c("preprocessing_id","estimated_parameters","fixed_parameters","priors","inference",
      "npi","expected_definition","beta_definition","representation","case_unit","beta_unit","training_scope")
    for(n in required) p3_text(md[[n]],paste("Mref metadata",n))
    if(!identical(md$representation,"pointwise_posterior_median") ||
       !identical(md$case_unit,"reported_cases_per_week") || !identical(md$beta_unit,"week^-1"))
      stop("Mref representation or units unsupported")
    p3_text(m$commit,"Mref commit");p3_text(m$model_sha256,"Mref model SHA256")
    if(!grepl("^[a-f0-9]{40}$",m$commit) || !grepl("^[a-f0-9]{64}$",m$model_sha256)) stop("Malformed Mref source identity")
    if(!is.list(m$settings) || !all(c("chains","iter","warmup","seed","adapt_delta","max_treedepth") %in% names(m$settings)))
      stop("Mref MCMC settings incomplete")
    if(nrow(chains)!=m$settings$chains || anyDuplicated(chains$chain) ||
       any(chains$draws<=0 | chains$divergences<0 | chains$depth_hits<0 | !is.finite(chains$ebfmi)))
      stop("Inconsistent Mref chain diagnostic dimensions/values")
    if(!identical(quality$bfmi_pass,all(chains$ebfmi>=0.3)) ||
       !identical(quality$divergences_pass,all(chains$divergences==0)) ||
       !identical(quality$treedepth_pass,all(chains$depth_hits==0)))
      stop("Mref quality flags disagree with chain diagnostics")
    if(nrow(m$config)!=1L || !identical(as.character(m$config$region),data$definition$region)) stop("Mref region mismatch")
    if(!identical(as.character(m$config$variant),"seasonal") ||
       !identical(as.character(m$config$model),data$definition$variant)) stop("Mref model variant unsupported for region")
    cases <- utils::read.csv(paths[["cases"]],stringsAsFactors=FALSE)
    cases <- cases[cases$island==data$definition$region,,drop=FALSE]
    pops <- utils::read.csv(paths[["population"]],stringsAsFactors=FALSE)
    pop <- pops$pop[pops$island==data$definition$region]
    C <- utils::read.csv(paths[["C"]]);beta <- utils::read.csv(paths[["beta"]])
    if(!all(c("index","median","year","week") %in% names(C)) || !all(c("index","median") %in% names(beta)))
      stop("Unsupported Mref CSV schema")
    if(any(C$index!=floor(C$index)) || any(beta$index!=floor(beta$index)) ||
       nrow(C)!=nrow(cases) || !identical(as.integer(C$index),seq_len(nrow(C))) ||
       !identical(as.integer(beta$index),seq_len(52L))) stop("Mref period/dimension mismatch")
    n <- length(data$cases)
    if(anyNA(cases) || anyNA(pops) || anyDuplicated(paste(cases$year,cases$week)) ||
       any(diff(52*cases$year+cases$week)!=1) || cases$year[1]!=2013 || cases$week[1]!=1)
      stop("Invalid Mref calendar or missing source values")
    for(summary in list(C,beta)) if(all(c("lower","upper") %in% names(summary))) {
      finite_numeric(summary$lower,"Mref lower interval",nrow(summary))
      finite_numeric(summary$upper,"Mref upper interval",nrow(summary))
      if(any(summary$lower>summary$median | summary$median>summary$upper | summary$lower<0))
        stop("Invalid Mref posterior interval ordering")
    }
    if(nrow(C)<n || !identical(as.integer(C$year),as.integer(cases$year)) ||
       !identical(as.integer(C$week),as.integer(cases$week)) ||
       !identical(as.integer(cases$year[seq_len(n)]),as.integer(data$time_markers$year)) ||
       !identical(as.integer(cases$week[seq_len(n)]),as.integer(data$seasonal_week))) stop("Mref period mismatch")
    if(!isTRUE(all.equal(as.numeric(cases$cases[seq_len(n)]),data$cases,tolerance=0)) ||
       !isTRUE(all.equal(as.numeric(s$cases),as.numeric(cases$cases),tolerance=0)) ||
       !isTRUE(all.equal(as.numeric(s$week),as.numeric(cases$week),tolerance=0))) stop("Mref preprocessing/case mismatch")
    if(length(pop)!=1L || pop!=data$fixed_parameters$pop || s$pop!=pop ||
       s$mu!=data$fixed_parameters$mu || s$gamma!=data$fixed_parameters$gamma || s$N!=nrow(C) ||
       s$Nnonpi!=data$definition$Nnonpi) stop("Mref population/rates/NPI indexing mismatch")
    finite_numeric(C$median,"Mref expected medians",nrow(C));finite_numeric(beta$median,"Mref beta medians",52L)
    if(C$median[1]!=0) stop("Mref initial expected-case sentinel differs from original Stan C[1]=0")
    if(any(C$median<0) || any(beta$median<=0)) stop("Invalid Mref expected/beta values")
    audit <- data.frame(field=c("model","commit","model_sha256","region","training_rows","training_period",
      "population","mu","gamma","Nnonpi","likelihood_mask","MCMC_settings","convergence",required),
      value=c(as.character(m$config$model),m$commit,m$model_sha256,as.character(m$config$region),as.character(nrow(C)),
        paste(paste(C$year[1],C$week[1],sep="-W"),paste(tail(C$year,1),tail(C$week,1),sep="-W"),sep=" through "),
        as.character(s$pop),as.character(s$mu),as.character(s$gamma),as.character(s$Nnonpi),"Original rows 2..N",
        paste(capture.output(dput(m$settings)),collapse=" "),diagnostic_status,unlist(md[required],use.names=FALSE)))
    mismatch <- md$preprocessing_id!=data$provenance$preprocessing_id
    if(mismatch) return(unavailable("Mref preprocessing identity differs",audit,paths))
    list(status="AVAILABLE_DESCRIPTIVE_ONLY",reason="Bayesian priors/inference and nuisance/NPI estimation differ from penalized M0/M1",
      audit=audit,manifest=m,standata=s,metadata=md,paths=paths,
      expected=C$median[seq_len(n)],beta=beta$median,markers=p3_markers(data),
      training_rows=seq_len(nrow(C)),forecast_eligible=FALSE,
      diagnostic_status=diagnostic_status,quality=quality,chain_diagnostics=chains,
      expected_summary=C[seq_len(n),,drop=FALSE],beta_summary=beta,
      uncertainty="Pointwise latent expected-case / beta posterior intervals, not observation predictive intervals",
      interpretation="Pointwise posterior medians, not a joint plug-in fit or posterior predictive score; convergence not certified by this adapter")
  },error=function(e) unavailable(conditionMessage(e),paths=paths))
}
