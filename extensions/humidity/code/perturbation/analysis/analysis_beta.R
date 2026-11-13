# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("analysis")

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="serif"))
library(egg)
library(rstan)

load("../stanfit_sirs2/stanfit_sirs_honshu_npi2.rda")
load("../stanfit_sirs2/stanfit_sirs_kyushu_npi2.rda")
load("../stanfit_sirs2/stanfit_sirs_shikoku_npi2.rda")

ee_honshu <- rstan::extract(stanfit_sirs_honshu_npi2)
ee_kyushu <- rstan::extract(stanfit_sirs_kyushu_npi2)
ee_shikoku <- rstan::extract(stanfit_sirs_shikoku_npi2)

quantile(sapply(1:length(ee_honshu$S0), function(i) {
  mean(ee_honshu$beta2[i,])/mean(ee_honshu$beta[i,])
}), c(0.025, 0.5, 0.975))

quantile(sapply(1:length(ee_kyushu$S0), function(i) {
  mean(ee_kyushu$beta2[i,])/mean(ee_kyushu$beta[i,])
}), c(0.025, 0.5, 0.975))

quantile(sapply(1:length(ee_shikoku$S0), function(i) {
  mean(ee_shikoku$beta2[i,])/mean(ee_shikoku$beta[i,])
}), c(0.025, 0.5, 0.975))
