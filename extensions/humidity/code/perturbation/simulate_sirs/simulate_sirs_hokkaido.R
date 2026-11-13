# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("simulate_sirs")

library(rstan)
source("../script/script_data.R")
source("../R/simulate_sirs.R")

load("../stanfit_sirs/stanfit_sirs_hokkaido_npi.rda")

ss <- summary(stanfit_sirs_hokkaido_npi)

pop <- filter(japan_pop_island, island=="Hokkaido")$pop

week <- c(25:52, 1:25)
mu <- 1/80/52
gamma <- 1
beta <- ss$summary[grepl("beta\\[",rownames(ss$summary)),6]
omega <- ss$summary[grepl("omega",rownames(ss$summary)),6]
delta <- ss$summary[grepl("delta",rownames(ss$summary)),6]
rho <- ss$summary[grepl("rho",rownames(ss$summary)),6]
npieff <- ss$summary[grepl("npieff\\[",rownames(ss$summary)),6]
S0 <- ss$summary[grepl("S0",rownames(ss$summary)),6]
I0 <- ss$summary[grepl("I0",rownames(ss$summary)),6]

S0vec <- seq(0.078, 0.09, length.out=21)
I0vec <- exp(seq(log(3e-4), log(3e-3), length.out=21))

paramdata <- expand.grid(S0vec, I0vec)

simlist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
  pp <- paramdata[i,]
  
  sim <- simulate_sirs(beta=beta, 
                week=week,
                mu=mu,
                gamma=gamma,
                omega=omega,
                delta=delta,
                rho=rho,
                npieff=npieff,
                pop=pop,
                S0=pp[[1]],
                I0=pp[[2]],
                tmax=53)
  
  cases <- sim$cases[-1]
  p <- cases/sum(cases)
  
  cog <- sum(cases*1:52)/sum(cases)+26
  intensity <- -1/sum(p * log(p))
  
  simlist[[i]] <- data.frame(
    week=c(26:52, 1:25),
    cases=cases,
    S0=pp[[1]],
    I0=pp[[2]],
    island="Hokkaido"
  )
}

simulate_sirs_hokkaido <- simlist %>%
  bind_rows()

save("simulate_sirs_hokkaido", file="simulate_sirs_hokkaido.rda")
