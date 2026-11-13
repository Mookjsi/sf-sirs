# Phase 2 configuration and structural checks; no model evaluation occurs on source().
validate_estimation_config <- function(config) {
  if(!is.list(config)) stop("estimation_config must be a list",call.=FALSE)
  if(!identical(sort(names(config$free)),sort(c("S0","I0","rho"))) ||
     !all(vapply(config$free,function(x) is.logical(x) && length(x)==1L && !is.na(x),logical(1))))
    stop("free must explicitly contain logical S0, I0 and rho flags",call.=FALSE)
  if(!is.character(config$method) || length(config$method)!=1L || is.na(config$method) ||
     !config$method %in% c("L-BFGS-B","BFGS","Nelder-Mead")) stop("Set a supported optimizer method",call.=FALSE)
  if(!is.null(config$joint)) {
    if(!identical(sort(names(config$joint)),sort(c("omega","delta","mu","gamma"))) ||
       !all(vapply(config$joint,function(x) is.logical(x) && length(x)==1L && !is.na(x),logical(1))))
      stop("joint must explicitly declare omega/delta/mu/gamma free flags",call.=FALSE)
  }
  allowed <- if(config$method=="L-BFGS-B") c("maxit","trace","factr","pgtol","lmm") else c("maxit","trace","reltol","abstol","REPORT")
  if(!is.list(config$control)) stop("optim control must be a named list",call.=FALSE)
  control_names <- names(config$control)
  if(length(config$control) && (is.null(control_names) || anyNA(control_names) ||
     any(!nzchar(control_names)) || anyDuplicated(control_names)))
    stop("optim control names must be unique and nonempty",call.=FALSE)
  if(any(!control_names %in% allowed)) stop("Unknown or method-incompatible optim control",call.=FALSE)
  for(n in names(config$control)) scalar_number(config$control[[n]],paste("control",n),0)
  for(n in intersect(names(config$control),c("maxit","trace","REPORT","lmm")))
    integer_scalar(config$control[[n]],paste("control",n),if(n %in% c("REPORT","lmm")) 1L else 0L)
  scalar_number(config$failure_value,"failure_value",0,strict_lower=TRUE)
  scalar_number(config$boundary_tolerance,"boundary_tolerance",0)
  integer_scalar(config$failure_examples,"failure_examples",0L)
  for(n in c("log_positive","spline","logit")) {
    b <- config$bounds[[n]]; finite_numeric(b,paste("bounds",n),2L)
    if(b[1]>=b[2]) stop("Bounds must be ordered lower < upper",call.=FALSE)
  }
  finite_numeric(config$ndeps,"ndeps"); finite_numeric(config$parscale,"parscale")
  if(any(config$ndeps<=0) || any(config$parscale<=0)) stop("ndeps and parscale must be positive",call.=FALSE)
  invisible(config)
}
validate_phase2_data <- function(data) {
  if(!inherits(data,"spline_sirs_data")) stop("Phase 1 adapted data are required",call.=FALSE)
  n <- length(data$cases)
  if(n<2L) stop("At least one likelihood observation after the initial row is required",call.=FALSE)
  finite_numeric(data$cases,"cases",n)
  if(any(data$cases<0 | data$cases!=floor(data$cases))) stop("Invalid case counts",call.=FALSE)
  w <- seasonal_weeks(data$seasonal_week)
  if(length(w)!=n || !identical(data$observation_order,seq_len(n)) ||
     !identical(data$likelihood_include,c(FALSE,rep(TRUE,n-1L)))) stop("Phase 1 ordering/mask contract mismatch",call.=FALSE)
  finite_numeric(data$npieff,"npieff",n)
  boundary <- data$definition$Nnonpi
  tail <- if(n>boundary) data$npieff[seq.int(boundary+1L,n)] else NULL
  if(!identical(as.numeric(data$npieff),as.numeric(original_npi_effect(n,boundary,tail))))
    stop("Original NPI prefix/tail contract mismatch",call.=FALSE)
  scalar_number(data$fixed_parameters$pop,"pop",0,strict_lower=TRUE)
  scalar_number(data$fixed_parameters$mu,"mu",0,strict_lower=TRUE)
  scalar_number(data$fixed_parameters$gamma,"gamma",0)
  if(nrow(data$time_markers)!=n) stop("Observation markers are missing or misaligned",call.=FALSE)
  invisible(data)
}
expand_control_vector <- function(x,n,label) {
  if(length(x)==1L) x <- rep(x,n)
  finite_numeric(x,label,n)
  if(any(x<=0)) stop(label," must be positive",call.=FALSE)
  as.numeric(x)
}
prepare_estimation_problem <- function(data,model_config,estimation_config) {
  validate_phase2_data(data); validate_estimation_config(estimation_config)
  integer_scalar(model_config$K,"K",4L); scalar_number(model_config$lambda,"lambda",0)
  if(!is.null(model_config$transmission) && !model_config$transmission %in% c("constant","spline","fourier"))
    stop("Unknown transmission representation",call.=FALSE)
  if(identical(model_config$transmission,"constant") && model_config$lambda!=0)
    stop("M0 requires zero penalty",call.=FALSE)
  fourier <- identical(model_config$transmission,"fourier")
  if(fourier && model_config$lambda!=0) stop("Fourier fits require zero penalty",call.=FALSE)
  basis <- if(fourier) new_fourier_basis(model_config$harmonics) else new_cyclic_basis(model_config$K)
  joint <- estimation_config$joint
  for(n in c("omega","delta")) if(!isTRUE(joint[[n]])) scalar_number(model_config[[n]],n,0)
  if(!is.null(joint)) for(n in c("mu","gamma")) model_config[[n]] <- data$fixed_parameters[[n]]
  for(n in c("S0","I0","rho")) if(!estimation_config$free[[n]])
    scalar_number(model_config[[n]],paste("fixed",n),0,1)
  required_sources <- c(c("omega","delta")[!vapply(c("omega","delta"),function(n) isTRUE(joint[[n]]),logical(1))],
    "npieff",names(estimation_config$free)[!unlist(estimation_config$free)])
  for(n in required_sources) {
    s <- model_config$fixed_sources[[n]]
    if(!is.character(s) || length(s)!=1L || is.na(s) || !nzchar(s))
      stop("Record an explicit fixed-value source for ",n,call.=FALSE)
  }
  f <- estimation_config$free
  if(!f$S0 && !f$I0 && model_config$S0+model_config$I0>1)
    stop("Fixed initial fractions exceed one",call.=FALSE)
  if(f$S0 && !f$I0 && model_config$I0>=1) stop("No mass remains for free S0",call.=FALSE)
  if(f$I0 && !f$S0 && model_config$S0>=1) stop("No mass remains for free I0",call.=FALSE)
  map <- make_parameter_map(model_config,estimation_config,basis)
  context <- list(cases=data$cases,weeks=data$seasonal_week,npieff=data$npieff,
    mask=data$likelihood_include,markers=data$time_markers,fixed=data$fixed_parameters,
    model=model_config,free=estimation_config$free,K=basis$K,lambda=model_config$lambda)
  list(data=data,model=model_config,estimation=estimation_config,basis=basis,map=map,
       comparison_context=context,n_observations=sum(data$likelihood_include),
       ndeps=expand_control_vector(estimation_config$ndeps,length(map$order),"ndeps"),
       parscale=expand_control_vector(estimation_config$parscale,length(map$order),"parscale"))
}
