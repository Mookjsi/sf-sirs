NPIEFF PAIRED COMPARISON: 2013-2017 training / 2018-2019 conditional hindcast

Run:
1. C:/Users/c/Desktop/analysis/fit_baseline_and_humidity_npieff_2013_2017.R
   Fits Baseline and humidity_seasonal independently for Hokkaido, Shikoku, Ryukyu.
   A matching saved fit is reused. force_refit=TRUE requests a new MCMC fit.
   No automatic retries or tuning. Failed diagnostics stop that job's prediction
   but the remaining jobs continue. Existing caches/results are backed up on write.
2. C:/Users/c/Desktop/analysis/Baseline vs Humidity anomaly.R
   Reads only this new workflow's six successful outputs. Displays g_prediction.
   A failed or interrupted latest attempt prevents reuse of an older result.

Individual entrypoints:
- fit_baseline_three_islands_2013_2017.R: all three baselines.
- 2013-2017 training/humidity_seasonal_hokkaido.R: Hokkaido humidity fit.
- 2013-2017 training/humidity_seasonal_sikoku.R: Shikoku humidity fit.
- 2013-2017 training/humidity_seasonal_ryukyu.R: Ryukyu humidity fit.
- 2013-2017 training/plot_hokkaido_shikoku_ryukyu_2018_2019.R:
  the three new humidity fits only.

Input validation without sampling/compilation/output writes:
options(npieff.validation_only = TRUE)
source("C:/Users/c/Desktop/analysis/fit_baseline_and_humidity_npieff_2013_2017.R",
       encoding = "UTF-8")
options(npieff.validation_only = FALSE)

Model:
Baseline effective transmission: beta[week] * npieff[t].
Humidity effective transmission: beta[week] * npieff[t] *
 exp((theta0 + theta_s*sin(2*pi*week/52) + theta_c*cos(2*pi*week/52))*z[t]).
Both use the original SIRS transition, cyclic beta prior, common parameter priors,
and negative-binomial observation model. Humidity means and anomaly scale use
2013-2017 only; actual 2018-2019 humidity is a known conditional forecast input.

Original island-specific Nnonpi settings are retained:
Hokkaido: Nnonpi=208, npieff varies from 2017; 52 active training weeks and
104 active forecast weeks.
Ryukyu: Nnonpi=312, npieff varies from 2019; no active training weeks and
52 active forecast weeks. The 2019 path/scale is prior driven, not estimated
from 2019 reported cases.
Shikoku: Nnonpi=364, npieff varies from 2020. It stays exactly 1 throughout
this analysis; the original Shikoku no-anchor prior exception is retained
in the code but no active npi coefficients are required for 2013-2019.

Future latent npieff values are sampled jointly with training parameters under
the original finite-path prior: positive support, per-week Normal(1,0.2)
anchor for Hokkaido/Ryukyu, adjacent Normal(npi_previous,sigma_npi) factors,
and sigma_npi positive Normal(0,0.1). They are NOT fitted to test cases.
Noncentering the increments preserves that joint density, including positivity.
It does not add sequential truncated-normal normalizing constants.

The path prior horizon is explicitly 2013-2019 (364 weeks), matching this task.
The original study used 2013-2023 (572 weeks). The functional form is retained;
this is not the marginal of the longer 2023-horizon prior and is not a numerical
reproduction of the original full-data posterior. The finite-horizon anchors
also regularize the latent path; future uncertainty is not an unconstrained
random walk appended after training.

Stan receives ONLY 260 training cases, weeks, and humidity design rows.
The held-out SIRS paths are propagated separately from each fit's own 2017-end
S/I/R states, each draw's parameters and future npieff, continuously for 104 weeks.
Prediction intervals include negative-binomial observation noise.
The diagnostics gate checks all nonconstant fitted quantities, including latent
future npieff; excludes ONLY C[1] and the npieff values fixed mathematically to 1.
Required: Rhat <= 1.01, n_eff >=400, no divergences or maximum-depth hits,
and each chain E-BFMI >=0.3. These thresholds are checks, not a convergence guarantee.

New outputs:
PROJECT/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/
  {hokkaido,shikoku,ryukyu}/{baseline,humidity_seasonal}/
    fit_cache.rds, result.rds, latest_attempt.rds
Old no-npieff result folders are never used or overwritten.
Existing A/B peak-shift scripts have not been migrated to these new caches.

Setup was verified without MCMC or C++ compilation.

