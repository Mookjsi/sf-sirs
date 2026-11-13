# NB2 with mean m and size phi, including the exact zero-mean boundaries.
nb2_loglik <- function(y,mean,phi,include=NULL) {
  finite_numeric(y,"y"); finite_numeric(mean,"mean",length(y))
  scalar_number(phi,"phi",0,strict_lower=TRUE)
  if(any(y<0 | y!=floor(y)) || any(mean<0)) stop("NB2 requires nonnegative integer counts and nonnegative means",call.=FALSE)
  if(!is.logical(include) || length(include)!=length(y) || anyNA(include))
    stop("Supply the explicit likelihood inclusion mask (original model excludes row 1)",call.=FALSE)
  logp <- numeric(length(y))
  zero <- mean==0
  logp[zero & y>0] <- -Inf
  k <- which(!zero)
  if(length(k)) {
    lm <- log(mean[k]); lp <- log(phi)
    ld <- pmax(lm,lp)+log1p(exp(-abs(lm-lp)))
    logp[k] <- lgamma(y[k]+phi)-lgamma(phi)-lgamma(y[k]+1)+
      phi*(lp-ld)+y[k]*(lm-ld)
  }
  if(anyNA(logp) || any(logp==Inf)) stop("NB2 arithmetic produced a nonfinite positive/undefined log probability; inputs are not clipped",call.=FALSE)
  contribution <- logp; contribution[!include] <- 0
  list(per_week=logp,included_per_week=contribution,total=sum(contribution),include=include,
       parameterization="NB2: mean=m, variance=m+m^2/phi; phi is not rho")
}
