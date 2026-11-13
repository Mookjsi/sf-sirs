# Centered log-Fourier forcing on the same x=w-1 calendar as the cyclic spline.
fourier_names <- function(harmonics) {
  integer_scalar(harmonics,"harmonics",1L)
  if(!harmonics %in% 1:2) stop("Supported Fourier orders are 1 (A) and 2 (AS)")
  as.vector(rbind(paste0("a",seq_len(harmonics)),paste0("b",seq_len(harmonics))))
}
new_fourier_basis <- function(harmonics) {
  names <- fourier_names(harmonics)
  structure(list(K=length(names)+1L,harmonics=as.integer(harmonics),period=52,
    coefficient_names=names,weekly_coordinates=0:51,week_coordinate="x_w = w - 1",
    normalization="log_sum_exp(g(0:51)) - log(52); discrete arithmetic mean"),class="fourier_basis")
}
fourier_log_shape <- function(x,coefficients,basis,derivative=0L) {
  if(!inherits(basis,"fourier_basis")) stop("Construct a Fourier basis first")
  finite_numeric(x,"Fourier coordinates");finite_numeric(coefficients,"Fourier coefficients",basis$K-1L)
  integer_scalar(derivative,"Fourier derivative",0L)
  if(derivative>2L) stop("Only Fourier derivatives 0..2 supported")
  if(!is.null(names(coefficients)) && !identical(names(coefficients),basis$coefficient_names))
    stop("Fourier coefficient ordering must be ",paste(basis$coefficient_names,collapse=","))
  x <- x %% 52;g <- numeric(length(x))
  for(h in seq_len(basis$harmonics)) {
    frequency <- 2*pi*h/52;angle <- frequency*x
    a <- coefficients[2L*h-1L];b <- coefficients[2L*h]
    term <- if(derivative==0L) a*cos(angle)+b*sin(angle) else
      if(derivative==1L) frequency*(-a*sin(angle)+b*cos(angle)) else
        -frequency^2*(a*cos(angle)+b*sin(angle))
    g <- g+term
  }
  if(any(!is.finite(g))) stop("Nonfinite Fourier log shape")
  g
}
fourier_profile <- function(coefficients,beta0,basis,x=basis$weekly_coordinates,derivative=0L) {
  scalar_number(beta0,"beta0",0,strict_lower=TRUE)
  weekly <- fourier_log_shape(basis$weekly_coordinates,coefficients,basis)
  maximum <- max(weekly)
  normalizer <- maximum+log(sum(exp(weekly-maximum)))-log(52)
  beta <- exp(log(beta0)+fourier_log_shape(x,coefficients,basis)-normalizer)
  if(any(!is.finite(beta)) || any(beta<=0)) stop("Fourier transmission overflow/underflow; no clipping")
  integer_scalar(derivative,"beta derivative",0L)
  if(derivative>2L) stop("Only beta derivatives 0..2 supported")
  values <- beta
  if(derivative>=1L) {
    first <- fourier_log_shape(x,coefficients,basis,1L)
    values <- if(derivative==1L) beta*first else beta*(first^2+fourier_log_shape(x,coefficients,basis,2L))
  }
  if(any(!is.finite(values))) stop("Nonfinite Fourier derivative")
  list(beta_by_week=values,coefficients=coefficients,log_normalizer=normalizer,
    a0=log(beta0)-normalizer,beta0=beta0,basis=basis)
}
