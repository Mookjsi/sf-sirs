simulate_sirs <- function(beta, 
                          week,
                          mu,
                          gamma,
                          omega,
                          delta,
                          rho,
                          npieff,
                          pop,
                          S0,
                          I0,
                          tmax=572) {
  Svec <- Ivec <- Rvec <- Cvec <- rep(0, tmax)
  
  Svec[1] <- S0 * pop
  Ivec[1] <- I0 * pop
  Rvec[1] <- (1 - S0 - I0) * pop
  Cvec[1] <- NA
  
  for (i in 2:tmax) {
    foi <- beta[[week[i]]] * (Ivec[i-1]+omega)/pop * npieff[i]
    
    Sout <- (1-exp(-(foi+mu))) * Svec[i-1]
    StoI <- foi/(foi+mu) * Sout
    Iout <- (1-exp(-(gamma+mu))) * Ivec[i-1]
    ItoR <- gamma/(gamma+mu) * Iout
    Rout <- (1-exp(-(delta+mu))) * Rvec[i-1]
    RtoS <- delta/(delta+mu) * Rout
    
    Svec[i] <- Svec[i-1] + RtoS - Sout + mu * pop
    Ivec[i] <- Ivec[i-1] + StoI - Iout
    Rvec[i] <- Rvec[i-1] + ItoR - Rout
    Cvec[i] <- StoI * rho
  }
  
  data.frame(
    S=Svec,
    I=Ivec,
    R=Rvec,
    cases=Cvec
  )
}

simulate_sirs_stoch <- function(beta, 
                                week,
                                mu,
                                gamma,
                                omega,
                                delta,
                                rho,
                                npieff,
                                pop,
                                S0,
                                I0,
                                phi,
                                tmax=572) {
  Svec <- Ivec <- Rvec <- Cvec_true <- Cvec_binom <- Cvec_nbinom <- rep(0, tmax)
  
  Svec[1] <- round(S0 * pop)
  Ivec[1] <- round(I0 * pop)
  Rvec[1] <- round((1 - S0 - I0) * pop)
  
  for (i in 2:tmax) {
    foi <- beta[[week[i]]] * (Ivec[i-1]+omega)/pop * npieff[i]
    
    Sout <- rbinom(1, Svec[i-1], (1-exp(-(foi+mu))))
    StoI <- rbinom(1, Sout, foi/(foi+mu))
    Iout <- rbinom(1, Ivec[i-1], (1-exp(-(gamma+mu))))
    ItoR <- rbinom(1, Iout, gamma/(gamma+mu))
    Rout <- rbinom(1, Rvec[i-1], (1-exp(-(delta+mu))))
    RtoS <- rbinom(1, Rout, delta/(delta + mu))
    birth <- rpois(1, mu * pop)
    
    Svec[i] <- Svec[i-1] + RtoS - Sout + birth
    Ivec[i] <- Ivec[i-1] + StoI - Iout
    Rvec[i] <- Rvec[i-1] + ItoR - Rout
    Cvec_true[i] <- StoI * rho
    Cvec_binom[i] <- rbinom(1, StoI, rho)
    Cvec_nbinom[i] <- rnbinom(1, mu=StoI*rho, size=phi)
  }
  
  Cvec_true[1] <- NA
  Cvec_binom[1] <- NA
  Cvec_nbinom[1] <- NA
  
  data.frame(
    S=Svec,
    I=Ivec,
    R=Rvec,
    cases_true=Cvec_true,
    cases_binom=Cvec_binom,
    cases_nbinom=Cvec_nbinom
  )
}

simulate_sirs_stoch_extra <- function(beta, 
                                week,
                                mu,
                                gamma,
                                omega,
                                delta,
                                rho,
                                npieff,
                                pop,
                                S0,
                                I0,
                                phi,
                                sigma2,
                                corr,
                                tmax=572) {
  Svec <- Ivec <- Rvec <- Cvec <- rep(0, tmax)
  
  Svec[1] <- round(S0 * pop)
  Ivec[1] <- round(I0 * pop)
  Rvec[1] <- round((1 - S0 - I0) * pop)
  
  sigmamat <- matrix(0, tmax, tmax) 
  diag(sigmamat) <- 1
  
  for (x in 1:(tmax-1)) {
    sigmamat[x,x+1] <- corr
    sigmamat[x+1,x] <- corr
  }
  
  sigmamat <- sigmamat * sigma2
  
  epsilon <- c(mvtnorm::rmvnorm(1, sigma=sigmamat))
  
  for (i in 2:tmax) {
    foi <- beta[[week[i]]] * (Ivec[i-1]+omega)/pop * npieff[i] * exp(epsilon[i])
    
    Sout <- rbinom(1, Svec[i-1], (1-exp(-(foi+mu))))
    StoI <- rbinom(1, Sout, foi/(foi+mu))
    Iout <- rbinom(1, Ivec[i-1], (1-exp(-(gamma+mu))))
    ItoR <- rbinom(1, Iout, gamma/(gamma+mu))
    Rout <- rbinom(1, Rvec[i-1], (1-exp(-(delta+mu))))
    RtoS <- rbinom(1, Rout, delta/(delta + mu))
    birth <- rpois(1, mu * pop)
    
    Svec[i] <- Svec[i-1] + RtoS - Sout + birth
    Ivec[i] <- Ivec[i-1] + StoI - Iout
    Rvec[i] <- Rvec[i-1] + ItoR - Rout
    Cvec[i] <- rbinom(1, StoI, rho)
  }
  
  Cvec[1] <- NA
  
  data.frame(
    S=Svec,
    I=Ivec,
    R=Rvec,
    cases=Cvec
  )
}
