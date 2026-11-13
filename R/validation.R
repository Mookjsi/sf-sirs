# Shared input checks. Definitions only; no code is evaluated on source().
finite_numeric <- function(x, name, length_required=NULL) {
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || any(!is.finite(x)) ||
      (!is.null(length_required) && length(x)!=length_required))
    stop(name, " must be finite numeric with the required length", call.=FALSE)
  invisible(x)
}
scalar_number <- function(x,name,lower=-Inf,upper=Inf,strict_lower=FALSE) {
  finite_numeric(x,name,1L)
  below <- if(strict_lower) x<=lower else x<lower
  if(x>upper || below)
    stop(name," is outside its admissible range",call.=FALSE)
  invisible(x)
}
integer_scalar <- function(x,name,lower=1L) {
  scalar_number(x,name,lower)
  if(x!=floor(x) || x>.Machine$integer.max) stop(name," must be an integer",call.=FALSE)
  invisible(x)
}
seasonal_weeks <- function(week) {
  finite_numeric(week,"seasonal week")
  if(any(week!=floor(week) | week<1 | week>52))
    stop("Seasonal weeks must be integers 1..52. Run the original week-53 mean/round preprocessing first.",call.=FALSE)
  as.integer(week)
}
