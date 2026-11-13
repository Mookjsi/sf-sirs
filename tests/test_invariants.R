script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."))
setwd(root)
sys.source("R/load_phase6.R", globalenv())
load_spline_sirs_phase6(root, globalenv())
for (f in c("R/phase5_reporting.R", "R/phase10_peak_identity.R",
            "R/result_storage.R", "R/accepted_extension.R", "R/phase21.R",
            "R/phase22.R", "R/phase22_analysis.R", "R/phase24.R",
            "R/phase25.R", "R/phase26.R", "R/phase16_targeted.R",
            "R/phase27.R")) sys.source(f, globalenv())

saved <- c(S=450, I=50, R=500)
s <- p27_restart_state(saved, u=log(1.2))
i <- p27_restart_state(saved, z=log(1.2))
stopifnot(identical(unname(s["I"]), unname(saved["I"])),
          identical(unname(i["S"]), unname(saved["S"])),
          abs(sum(s)-sum(saved)) < 1e-12,
          abs(sum(i)-sum(saved)) < 1e-12,
          abs(s["S"]/saved["S"]-1.2) < 1e-12,
          abs(i["I"]/saved["I"]-1.2) < 1e-12)

fixed <- list(pop=1000, mu=0.001, gamma=0.5,
              omega=1, delta=0.01, rho=0.2)
beta <- rep(0.5, 52)
time <- 1:6
baseline <- simulate_weekly_sirs(beta, saved, fixed, time,
                                  rep(1, length(time)))
perturbed <- simulate_weekly_sirs(beta, saved, fixed, time,
                                   rep(1, length(time)),
                                   beta_by_transition=c(beta[1], rep(1.2*beta[1], 5)))
stopifnot(all(baseline[c("S", "I", "R")] >= 0),
          all(perturbed[c("S", "I", "R")] >= 0),
          all(abs(baseline$expected_cases - baseline$incidence*fixed$rho) < 1e-12),
          perturbed$incidence[2] > baseline$incidence[2],
          identical(unlist(perturbed[1, c("S", "I", "R")], use.names=TRUE), saved))
cat("Scientific invariant checks passed\n")
