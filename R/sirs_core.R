# Discrete weekly RSV equations from Park et al.; see references/parksw3/LICENSE.
make_initial_state <- function(S0,I0,pop) {
  scalar_number(S0,"S0",0,1); scalar_number(I0,"I0",0,1)
  scalar_number(pop,"pop",0,strict_lower=TRUE)
  if(S0+I0>1) stop("S0 + I0 must not exceed one",call.=FALSE)
  state <- c(S=unname(S0*pop),I=unname(I0*pop),R=unname((1-S0-I0)*pop))
  if(any(!is.finite(state)) || any(state<0))
    stop("Initial-state arithmetic produced invalid counts; no clipping applied",call.=FALSE)
  state
}
simulate_weekly_sirs <- function(beta_by_week,initial_state,fixed_parameters,time_index,npieff,initial_mass=NULL,
                                 beta_by_transition=NULL) {
  finite_numeric(beta_by_week,"beta_by_week",52L)
  if(any(beta_by_week<0)) stop("beta_by_week must be nonnegative",call.=FALSE)
  w <- seasonal_weeks(time_index); n <- length(w)
  if(!is.null(beta_by_transition)) {
    finite_numeric(beta_by_transition,"beta_by_transition",n)
    if(any(beta_by_transition<0)) stop("beta_by_transition must be nonnegative",call.=FALSE)
  }
  finite_numeric(npieff,"npieff",n)
  if(any(npieff<0)) stop("npieff must be nonnegative",call.=FALSE)
  required <- c("pop","mu","gamma","omega","delta","rho")
  if(!is.list(fixed_parameters) || anyDuplicated(names(fixed_parameters)) ||
     !all(required %in% names(fixed_parameters))) stop("Set pop, mu, gamma, omega, delta and rho explicitly",call.=FALSE)
  p <- fixed_parameters
  scalar_number(p$pop,"pop",0,strict_lower=TRUE)
  scalar_number(p$mu,"mu",0,strict_lower=TRUE)
  for(k in c("gamma","omega","delta")) scalar_number(p[[k]],k,0)
  scalar_number(p$rho,"rho",0,1)
  finite_numeric(initial_state,"initial_state",3L)
  if(!identical(names(initial_state),c("S","I","R")) || any(initial_state<0))
    stop("initial_state must be nonnegative counts named S, I, R in that order",call.=FALSE)
  # Default Phase 1-3 contract is unchanged. Explicit continuation may preserve a
  # stored state's mass; the population denominator and every transition stay fixed.
  target_mass <- if(is.null(initial_mass)) p$pop else initial_mass
  scalar_number(target_mass,"initial_mass",0,strict_lower=TRUE)
  if(abs(sum(initial_state)-target_mass)>1e-10*max(1,target_mass))
    stop("Initial compartment counts must sum to pop or explicit initial_mass; this is not a later-time invariant",call.=FALSE)
  S <- I <- R <- incidence <- expected <- deaths_S <- deaths_I <- deaths_R <- numeric(n)
  S[1] <- initial_state["S"]; I[1] <- initial_state["I"]; R[1] <- initial_state["R"]
  # Initial incidence/expected cases are zero, matching Stan C[1], not the R helper's NA.
  if(n>=2L) for(i in 2:n) {
    beta_i <- if(is.null(beta_by_transition)) beta_by_week[w[i]] else beta_by_transition[i]
    foi <- beta_i*(I[i-1]+p$omega)/p$pop*npieff[i]
    if(!is.finite(foi)) stop("Nonfinite force of infection at row ",i,call.=FALSE)
    Sout <- (1-exp(-(foi+p$mu)))*S[i-1]
    StoI <- foi/(foi+p$mu)*Sout
    Iout <- (1-exp(-(p$gamma+p$mu)))*I[i-1]
    ItoR <- p$gamma/(p$gamma+p$mu)*Iout
    Rout <- (1-exp(-(p$delta+p$mu)))*R[i-1]
    RtoS <- p$delta/(p$delta+p$mu)*Rout
    S[i] <- S[i-1]+RtoS-Sout+p$mu*p$pop
    I[i] <- I[i-1]+StoI-Iout
    R[i] <- R[i-1]+ItoR-Rout
    incidence[i] <- StoI; expected[i] <- StoI*p$rho
    deaths_S[i] <- Sout-StoI; deaths_I[i] <- Iout-ItoR; deaths_R[i] <- Rout-RtoS
    if(any(!is.finite(c(S[i],I[i],R[i],incidence[i],expected[i]))) ||
       any(c(S[i],I[i],R[i],incidence[i],expected[i])<0))
      stop("Invalid state/incidence at row ",i,"; no clipping applied",call.=FALSE)
  }
  data.frame(row=seq_len(n),week=w,S=S,I=I,R=R,incidence=incidence,expected_cases=expected,
             deaths_S=deaths_S,deaths_I=deaths_I,deaths_R=deaths_R)
}
