# Parameter charts only. No Jacobian term belongs in penalized maximum likelihood.
make_parameter_map <- function(model,estimation,basis) {
  constant <- identical(model$transmission,"constant")
  fourier <- identical(model$transmission,"fourier")
  nz <- if(constant) 0L else basis$K-1L
  shape_order <- if(fourier) basis$coefficient_names else if(nz) paste0("z",seq_len(nz)) else character()
  order <- c("log_beta0",shape_order,"log_phi")
  kind <- c("log_positive",rep("spline",nz),"log_positive")
  f <- estimation$free
  if(f$S0) {order <- c(order,if(f$I0) "log_S_over_R" else "logit_S_remaining");kind <- c(kind,"logit")}
  if(f$I0) {order <- c(order,if(f$S0) "log_I_over_R" else "logit_I_remaining");kind <- c(kind,"logit")}
  if(f$rho) {order <- c(order,"logit_rho");kind <- c(kind,"logit")}
  joint <- estimation$joint
  if(!is.null(joint)) for(n in names(joint)) if(joint[[n]]) {
    order <- c(order,paste0("log_",n));kind <- c(kind,"log_positive")
  }
  lower <- vapply(kind,function(k) estimation$bounds[[k]][1],numeric(1))
  upper <- vapply(kind,function(k) estimation$bounds[[k]][2],numeric(1))
  names(lower) <- names(upper) <- order
  if(!is.null(estimation$coordinate_bounds)) for(n in names(estimation$coordinate_bounds)) {
    if(!n %in% order) next
    b <- estimation$coordinate_bounds[[n]];finite_numeric(b,paste("coordinate bounds",n),2L)
    if(b[1]>=b[2]) stop("Coordinate bounds must increase")
    lower[n] <- b[1];upper[n] <- b[2]
  }
  map <- list(order=order,kind=kind,dimension=length(order),lower=lower,upper=upper,
       free=f,fixed=list(S0=model$S0,I0=model$I0,rho=model$rho),K=basis$K,constant=constant,
       fourier=fourier,shape_order=shape_order,
       chart="Joint S/I/R log-ratios; conditional remaining-mass logit for one free fraction")
  if(!is.null(joint)) {
    map$joint <- joint;map$joint_mode <- TRUE
    for(n in names(joint)) map$fixed[[n]] <- model[[n]]
  }
  map
}
open_logit <- function(x,label) {
  scalar_number(x,label,0,1)
  if(x<=0 || x>=1) stop(label," must be strictly interior when free",call.=FALSE)
  stats::qlogis(x)
}
physical_value <- function(parameters,map,name) {
  if(map$free[[name]]) return(parameters[[name]])
  if(!is.null(parameters[[name]]) && !isTRUE(all.equal(unname(parameters[[name]]),unname(map$fixed[[name]]),tolerance=0)))
    stop("Initial parameter conflicts with fixed ",name,call.=FALSE)
  map$fixed[[name]]
}
encode_parameters <- function(parameters,map) {
  scalar_number(parameters$beta0,"beta0",0,strict_lower=TRUE)
  scalar_number(parameters$phi,"phi",0,strict_lower=TRUE)
  if(isTRUE(map$constant)) {
    if(length(parameters$free_coordinates) && any(parameters$free_coordinates!=0))
      stop("M0 has no seasonal coordinates",call.=FALSE)
    z <- numeric()
  } else {
    finite_numeric(parameters$free_coordinates,"free_coordinates",map$K-1L)
    if(isTRUE(map$fourier) && !is.null(names(parameters$free_coordinates)) &&
       !identical(names(parameters$free_coordinates),map$shape_order)) stop("Fourier coefficient ordering mismatch",call.=FALSE)
    z <- as.numeric(parameters$free_coordinates)
  }
  q <- c(log(parameters$beta0),z,log(parameters$phi))
  S <- physical_value(parameters,map,"S0"); I <- physical_value(parameters,map,"I0")
  scalar_number(S,"S0",0,1); scalar_number(I,"I0",0,1)
  R <- 1-S-I
  if(R<0) stop("Initial fractions are incompatible",call.=FALSE)
  if(map$free$S0 && map$free$I0) {
    if(min(S,I,R)<=0) stop("Joint free initial fractions must be interior",call.=FALSE)
    q <- c(q,log(S)-log(R),log(I)-log(R))
  } else {
    if(map$free$S0) q <- c(q,open_logit(S/(1-I),"S0 / remaining mass"))
    if(map$free$I0) q <- c(q,open_logit(I/(1-S),"I0 / remaining mass"))
  }
  rho <- physical_value(parameters,map,"rho"); scalar_number(rho,"rho",0,1)
  if(map$free$rho) q <- c(q,open_logit(rho,"rho"))
  if(isTRUE(map$joint_mode)) for(n in names(map$joint)) {
    if(map$joint[[n]]) {
      scalar_number(parameters[[n]],n,0,strict_lower=TRUE);q <- c(q,log(parameters[[n]]))
    } else if(!is.null(parameters[[n]]) && !isTRUE(all.equal(parameters[[n]],map$fixed[[n]],tolerance=0)))
      stop("Initial parameter conflicts with fixed ",n,call.=FALSE)
  }
  names(q) <- map$order
  check_optimizer_coordinates(q,map)
  q
}
check_optimizer_coordinates <- function(q,map) {
  finite_numeric(q,"optimizer coordinates",map$dimension)
  if(!is.null(names(q)) && !identical(names(q),map$order)) stop("Optimizer coordinate order mismatch",call.=FALSE)
  if(any(q<map$lower | q>map$upper)) stop("Optimizer coordinates outside explicit guard bounds",call.=FALSE)
  invisible(q)
}
decode_parameters <- function(q,map) {
  check_optimizer_coordinates(q,map); names(q) <- map$order
  S <- map$fixed$S0; I <- map$fixed$I0
  if(map$free$S0 && map$free$I0) {
    logs <- c(q[["log_S_over_R"]],q[["log_I_over_R"]],0)
    weights <- exp(logs-max(logs)); weights <- weights/sum(weights)
    S <- weights[1]; I <- weights[2]
  } else {
    if(map$free$S0) S <- (1-I)*stats::plogis(q[["logit_S_remaining"]])
    if(map$free$I0) I <- (1-S)*stats::plogis(q[["logit_I_remaining"]])
  }
  rho <- if(map$free$rho) stats::plogis(q[["logit_rho"]]) else map$fixed$rho
  shape_order <- if(is.null(map$shape_order)) paste0("z",seq_len(map$K-1L)) else map$shape_order
  p <- list(beta0=exp(q[["log_beta0"]]),free_coordinates=if(isTRUE(map$constant)) rep(0,map$K-1L) else unname(q[shape_order]),
            phi=exp(q[["log_phi"]]),S0=unname(S),I0=unname(I),R0=unname(1-S-I),rho=unname(rho))
  scalar_number(p$beta0,"beta0",0,strict_lower=TRUE); scalar_number(p$phi,"phi",0,strict_lower=TRUE)
  if(any(!is.finite(c(p$S0,p$I0,p$R0,p$rho))) || min(p$S0,p$I0,p$R0)<0 || p$rho<0 || p$rho>1)
    stop("Parameter transformation produced invalid fractions; no clipping",call.=FALSE)
  if(isTRUE(map$joint_mode)) {
    if(min(p$S0,p$I0,p$R0,p$rho,1-p$rho)<=0) stop("Joint estimation requires strictly interior fractions")
    for(n in names(map$joint)) p[[n]] <- if(map$joint[[n]]) exp(q[[paste0("log_",n)]]) else map$fixed[[n]]
    for(n in names(map$joint)) scalar_number(p[[n]],n,0,strict_lower=TRUE)
  }
  p
}
parameter_record <- function(problem,p) {
  f <- problem$estimation$free
  joint <- isTRUE(problem$map$joint_mode)
  val <- function(n) if(joint) p[[n]] else if(n %in% c("omega","delta")) problem$model[[n]] else problem$data$fixed_parameters[[n]]
  rows <- data.frame(parameter=c("beta0","phi","S0","I0","R0","rho","omega","delta","pop","mu","gamma"),
    status=c("free","free",if(f$S0) "free" else "fixed",if(f$I0) "free" else "fixed","derived",
      if(f$rho) "free" else "fixed",rep("fixed",5)),
    value=c(p$beta0,p$phi,p$S0,p$I0,p$R0,p$rho,val("omega"),val("delta"),
      problem$data$fixed_parameters$pop,val("mu"),val("gamma")),
    unit=c("week^-1","NB2 size","fraction","fraction","fraction","probability","individuals","week^-1","individuals","week^-1","week^-1"))
  rows$range <- c("positive; log guard","positive; log guard",rep("joint simplex [0,1]",3),"[0,1]",">=0",">=0",">0",">0",">=0")
  rows$source <- vapply(rows$parameter,function(n) {
    if(joint && n %in% names(problem$map$joint)) return(if(problem$map$joint[[n]]) "Independent joint estimation" else problem$model$fixed_sources[[n]])
    if(n %in% names(f) && !f[[n]]) return(problem$model$fixed_sources[[n]])
    if(n %in% c("omega","delta")) return(problem$model$fixed_sources[[n]])
    if(n %in% c("pop","mu","gamma")) return("Phase 1 adapted fixed_parameters")
    if(n=="R0") return("Remaining initial mass, not a separate parameter")
    "Explicit start / fitted value"
  },character(1))
  if(isTRUE(problem$map$joint_mode)) for(n in names(problem$map$joint)) {
    ids <- rows$parameter==n;rows$value[ids] <- p[[n]]
    rows$status[ids] <- if(problem$map$joint[[n]]) "free" else "fixed"
    rows$source[ids] <- if(problem$map$joint[[n]]) "Independent joint estimation" else problem$model$fixed_sources[[n]]
  }
  if(joint) {
    rows$range[rows$parameter %in% c("S0","I0","R0")] <- "(0,1); joint simplex"
    rows$range[rows$parameter=="rho"] <- "(0,1)"
    rows$range[rows$parameter %in% names(problem$map$joint)] <- ">0; numerical guard if free"
  }
  if(isTRUE(problem$map$constant)) return(rows)
  fourier <- isTRUE(problem$map$fourier)
  z <- data.frame(parameter=if(fourier) problem$basis$coefficient_names else paste0("z",seq_len(problem$basis$K-1L)),status="free",
    value=p$free_coordinates,unit=if(fourier) "Fourier log-shape coefficient" else "log seasonal shape",
    range=if(fourier) "explicit Fourier coefficient guard" else "explicit spline coordinate guard",source="Explicit start / fitted value")
  rbind(rows,z)
}
