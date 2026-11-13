# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("stanfit_sirs")

library(dplyr)
library(rstan)
rstan_options(auto_write=TRUE)
source("../script/script_data.R")

ii <- "Ryukyu"

cases <- filter(japan_rsv_island, island==ii)$cases
pop <- filter(japan_pop_island, island==ii)$pop

model <- stan_model("../stanmodel/sirs_npi.stan")

standata <- list(
  N=572,
  Nnonpi=312,
  week=filter(japan_rsv_island, island==ii)$week,
  cases=cases,
  mu=1/80/52,
  pop=pop,
  gamma=1
)

perturbation_validate(standata)

stanfit_sirs_ryukyu_npi <- sampling(model,
                                    data = standata,
                                    seed=102,
                                    iter=4000,
                                    chains=4,
                                    cores=4,
                                    control=list(max_treedepth=15,
                                                 adapt_delta=0.99))

check_hmc_diagnostics(stanfit_sirs_ryukyu_npi)
get_num_divergent(stanfit_sirs_ryukyu_npi)
get_num_max_treedepth(stanfit_sirs_ryukyu_npi)
get_low_bfmi_chains(stanfit_sirs_ryukyu_npi)

save("stanfit_sirs_ryukyu_npi", file="stanfit_sirs_ryukyu_npi.rda")

ss <- summary(stanfit_sirs_ryukyu_npi)

max(ss$summary[which(!is.na(ss$summary[,10])),10]) ## 1.004
min(ss$summary[which(!is.na(ss$summary[,10])),9]) ## 604

plot(cases)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),6])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("npieff\\[", rownames(ss$summary)),6], type="l")
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("beta\\[", rownames(ss$summary)),6], type="l")

plot(ss$summary[grepl("S\\[", rownames(ss$summary)),6]/pop)
