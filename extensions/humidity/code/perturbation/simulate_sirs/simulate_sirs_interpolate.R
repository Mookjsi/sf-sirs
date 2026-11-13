# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("simulate_sirs")

library(rstan)
library(mgcv)
source("../script/script_data.R")
source("../R/simulate_sirs.R")

load("../stanfit_sirs/stanfit_sirs_honshu_npi.rda")
load("../stanfit_sirs/stanfit_sirs_kyushu_npi.rda")

ss_honshu <- summary(stanfit_sirs_honshu_npi)
ss_kyushu <- summary(stanfit_sirs_kyushu_npi)

beta_honshu <- ss_honshu$summary[grepl("beta\\[",rownames(ss_honshu$summary)),6]
beta_kyushu <- ss_kyushu$summary[grepl("beta\\[",rownames(ss_kyushu$summary)),6]

week_gam <- 1:52

gfit_honshu <- gam(beta_honshu~s(week_gam, bs="cc"))
gfit_kyushu <- gam(beta_kyushu~s(week_gam, bs="cc"))

beta_smooth_honshu <- predict(gfit_honshu)
beta_smooth_kyushu <- predict(gfit_kyushu)

seas_smooth_honshu <- beta_smooth_honshu/mean(beta_smooth_honshu)-1
seas_smooth_kyushu <- beta_smooth_kyushu/mean(beta_smooth_kyushu)-1

amp_honshu <- (max(seas_smooth_honshu)-min(seas_smooth_honshu))/2
amp_kyushu <- (max(seas_smooth_kyushu)-min(seas_smooth_kyushu))/2

pop <- filter(japan_pop_island, island=="Honshu")$pop

week <- c(25:52, 1:25)
mu <- 1/80/52
gamma <- 1

omega <- ss_honshu$summary[grepl("omega",rownames(ss_honshu$summary)),6]
delta <- ss_honshu$summary[grepl("delta",rownames(ss_honshu$summary)),6]
rho <- ss_honshu$summary[grepl("rho",rownames(ss_honshu$summary)),6]
npieff <- ss_honshu$summary[grepl("npieff\\[",rownames(ss_honshu$summary)),6]
# S0 <- ss$summary[grepl("S0",rownames(ss$summary)),6]
# I0 <- ss$summary[grepl("I0",rownames(ss$summary)),6]

interp_vec <- seq(0, 1, length.out=41)
amp_vec <- seq(amp_honshu, amp_kyushu, length.out=41)
I0 <- exp(seq(log(1e-4), log(3e-3), length.out=21))[9]

paramdata <- expand.grid(interp_vec, amp_vec)

simlist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
  pp <- paramdata[i,]
  
  interp <- pp[[1]]
  amp <- pp[[2]]
  
  seas <- (seas_smooth_honshu/amp_honshu * (1-interp) + seas_smooth_kyushu/amp_kyushu * interp)
  
  amp_new <- (max(seas)-min(seas))/2
  
  beta <- (seas/amp_new * amp + 1) * mean(beta_smooth_honshu)
  
  sim1 <- simulate_sirs(beta=beta, 
                       week=week,
                       mu=mu,
                       gamma=gamma,
                       omega=omega,
                       delta=delta,
                       rho=rho,
                       npieff=npieff,
                       pop=pop,
                       S0=0.078,
                       I0=I0,
                       tmax=53)
  
  sim2 <- simulate_sirs(beta=beta, 
                        week=week,
                        mu=mu,
                        gamma=gamma,
                        omega=omega,
                        delta=delta,
                        rho=rho,
                        npieff=npieff,
                        pop=pop,
                        S0=0.105,
                        I0=I0,
                        tmax=53)
  
  cases1 <- sim1$cases[-1]
  cases2 <- sim2$cases[-1]
  
  simlist[[i]] <- data.frame(
    week=c(26:52, 1:25),
    shift_trans=which.max(c(seas, seas)[40:60])[[1]]+40-32,
    shift_cases=which.max(cases1)-which.max(cases2),
    interp=pp[[1]],
    amp=pp[[2]]
  )
}

simulate_sirs_interpolate <- simlist %>%
  bind_rows()

save("simulate_sirs_interpolate", file="simulate_sirs_interpolate.rda")
