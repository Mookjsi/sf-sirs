# Periodize uniform cubic B-splines from the established splines::splineDesign evaluator.
new_cyclic_basis <- function(K=8L,period=52) {
  integer_scalar(K,"K",4L); scalar_number(period,"period",0,strict_lower=TRUE)
  if(period!=52) stop("This weekly model fixes the period at 52",call.=FALSE)
  spacing <- period/K
  # K+7 extended knots yield K+3 local bases. Fold columns modulo K, retaining unity.
  structure(list(K=as.integer(K),period=period,degree=3L,order=4L,spacing=spacing,
    knots=(0:K)*spacing,extended_knots=(-3:(K+3))*spacing,
    local_to_cyclic=((seq_len(K+3)-1L) %% K)+1L,
    weekly_coordinates=0:51,week_coordinate="x_w = w - 1",
    constraint="c = c(z, -sum(z)); z has K-1 elements",
    normalization="log_sum_exp(g(0:51)) - log(52); discrete arithmetic mean",
    evaluator="splines::splineDesign; uniform extension and modulo-K column folding"),
    class="cyclic_cubic_basis")
}
require_cyclic_basis <- function(basis) {
  if(!inherits(basis,"cyclic_cubic_basis")) stop("Construct a basis with new_cyclic_basis()",call.=FALSE)
  invisible(basis)
}
cyclic_basis_matrix <- function(x,basis,derivative=0L) {
  require_cyclic_basis(basis); finite_numeric(x,"x")
  integer_scalar(derivative,"derivative",0L)
  if(derivative>2L) stop("Only value, first and second derivatives are exposed",call.=FALSE)
  wrapped <- x %% basis$period
  local <- splines::splineDesign(basis$extended_knots,wrapped,ord=4L,
                                 derivs=rep(as.integer(derivative),length(x)),outer.ok=FALSE)
  if(ncol(local)!=basis$K+3L) stop("Unexpected splineDesign dimension",call.=FALSE)
  B <- matrix(0,nrow=length(x),ncol=basis$K)
  for(j in seq_len(ncol(local))) {
    k <- basis$local_to_cyclic[j]; B[,k] <- B[,k]+local[,j]
  }
  if(any(!is.finite(B))) stop("Nonfinite basis evaluation",call.=FALSE)
  if(derivative==0L && any(abs(rowSums(B)-1)>1e-10))
    stop("Partition of unity check failed; do not use this basis",call.=FALSE)
  colnames(B) <- paste0("B",seq_len(basis$K))
  B
}
constrained_coefficients <- function(free_coordinates,basis) {
  require_cyclic_basis(basis)
  finite_numeric(free_coordinates,"free_coordinates",basis$K-1L)
  c <- c(as.numeric(free_coordinates),-sum(free_coordinates))
  if(any(!is.finite(c))) stop("Constraint transformation overflow",call.=FALSE)
  c
}
weekly_log_normalizer <- function(coefficients,basis) {
  require_cyclic_basis(basis)
  finite_numeric(coefficients,"coefficients",basis$K)
  g <- drop(cyclic_basis_matrix(basis$weekly_coordinates,basis) %*% coefficients)
  if(any(!is.finite(g))) stop("Nonfinite weekly log transmission",call.=FALSE)
  a <- max(g)
  a + log(sum(exp(g-a))) - log(52)
}
evaluate_spline_beta <- function(x,coefficients,beta0,basis,derivative=0L) {
  require_cyclic_basis(basis); finite_numeric(coefficients,"coefficients",basis$K)
  scalar_number(beta0,"beta0",0,strict_lower=TRUE)
  integer_scalar(derivative,"derivative",0L)
  if(derivative>2L) stop("Derivative must be 0, 1 or 2",call.=FALSE)
  log_normalizer <- weekly_log_normalizer(coefficients,basis)
  g <- drop(cyclic_basis_matrix(x,basis) %*% coefficients)
  beta <- exp(log(beta0)+g-log_normalizer)
  if(any(!is.finite(beta)) || any(beta<=0))
    stop("Transmission overflow/underflow; coefficients are not clipped",call.=FALSE)
  result <- beta
  if(derivative>=1L) {
    d1 <- drop(cyclic_basis_matrix(x,basis,1L) %*% coefficients)
    result <- beta*d1
    if(derivative==2L) {
      d2 <- drop(cyclic_basis_matrix(x,basis,2L) %*% coefficients)
      result <- beta*(d2+d1^2)
    }
  }
  if(any(!is.finite(result))) stop("Nonfinite transmission derivative",call.=FALSE)
  result
}
spline_beta_by_week <- function(free_coordinates,beta0,basis) {
  coefficients <- constrained_coefficients(free_coordinates,basis)
  # One shared transformation is used for continuous values and weekly output.
  evaluate_spline_beta(basis$weekly_coordinates,coefficients,beta0,basis)
}
spline_profile <- function(free_coordinates,beta0,basis) {
  coefficients <- constrained_coefficients(free_coordinates,basis)
  list(beta_by_week=spline_beta_by_week(free_coordinates,beta0,basis),
       coefficients=coefficients,free_coordinates=free_coordinates,beta0=beta0,
       log_normalizer=weekly_log_normalizer(coefficients,basis),basis=basis)
}
cyclic_second_difference <- function(coefficients) {
  finite_numeric(coefficients,"coefficients")
  K <- length(coefficients)
  if(K<4L) stop("Cubic coefficient coordinates require K >= 4",call.=FALSE)
  coefficients[c(2:K,1L)]-2*coefficients+coefficients[c(K,seq_len(K-1L))]
}
cyclic_penalty <- function(coefficients,lambda) {
  scalar_number(lambda,"lambda",0)
  finite_numeric(coefficients,"coefficients")
  if(length(coefficients)<4L) stop("K must be at least 4",call.=FALSE)
  if(lambda==0) return(0)
  d2 <- cyclic_second_difference(coefficients)
  value <- (lambda/2)*sum(d2^2)
  if(!is.finite(value)) stop("Penalty overflow; coefficients are not clipped",call.=FALSE)
  value
}
