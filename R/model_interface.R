# Both transmission representations use exactly the same state and observation code.
evaluate_beta_sirs <- function(data,beta_by_week,initial_state,dynamic_parameters,phi) {
  if(!inherits(data,"spline_sirs_data")) stop("Use adapt_preprocessed_island() first",call.=FALSE)
  if(!is.list(dynamic_parameters) || !setequal(names(dynamic_parameters),c("omega","delta","rho")) ||
     anyDuplicated(names(dynamic_parameters)))
    stop("Set omega, delta (=1/tau) and rho explicitly; no estimates are guessed",call.=FALSE)
  n <- length(data$cases)
  nonpi <- data$definition$Nnonpi
  tail_npi <- if(n>nonpi) data$npieff[seq.int(nonpi+1L,n)] else NULL
  expected_npi <- original_npi_effect(n,nonpi,tail_npi)
  if(!identical(as.numeric(data$npieff),as.numeric(expected_npi)))
    stop("NPI prefix must remain one under the original region definition",call.=FALSE)
  if(!identical(data$likelihood_include,c(FALSE,rep(TRUE,n-1L))))
    stop("The original likelihood must exclude only the first observation",call.=FALSE)
  states <- simulate_weekly_sirs(beta_by_week,initial_state,
    c(data$fixed_parameters,dynamic_parameters),data$seasonal_week,data$npieff)
  likelihood <- nb2_loglik(data$cases,states$expected_cases,phi,data$likelihood_include)
  list(states=states,likelihood=likelihood,beta_by_week=beta_by_week,
       observation_order=data$observation_order,time_markers=data$time_markers)
}
evaluate_spline_sirs <- function(data,free_coordinates,beta0,basis,initial_state,dynamic_parameters,phi) {
  profile <- spline_profile(free_coordinates,beta0,basis)
  result <- evaluate_beta_sirs(data,profile$beta_by_week,initial_state,dynamic_parameters,phi)
  result$spline <- profile
  result
}
