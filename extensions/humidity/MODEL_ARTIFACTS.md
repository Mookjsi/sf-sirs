# Model artifacts and frozen-state limitations

All paths below are relative to the original `perturbation` source root unless explicitly labelled `analysis`. Use file_inventory.csv and the artifact/archive manifests to resolve the delivered path and original archival path. Nothing was refit, deserialized through R, rounded, or numerically extracted during cleanup.

## What primary delivery supports

The primary package retains code, input tables, saved numerical summaries/figures, the compact result objects listed below, and six intact no-npieff humidity objects: a fitted cache plus result for Hokkaido, Shikoku and Ryukyu. All retained objects preserve their original bytes. The humidity caches retain full posterior draws and reference states according to their producing code; the owner has not certified them as final accepted fits. Trajectories and posterior curve summaries remain distinct from fitted parameters and exact initial states.

No standalone approved compact parameter-vector export was identified. Instead, Package A preserves three original no-npieff humidity fitted caches and their three saved results, supporting the recorded humidity replay and peak-analysis paths without archival fit files, subject to owner runtime validation, source/cache compatibility and the installed environment. Full Baseline + npieff versus humidity comparison and original manuscript fitted-model simulations still require specific archived fitted objects. No compact state extraction was performed.

## Current three-island comparison: primary and archival requirements

The available evidence identifies a 2013–2017 training / 2018–2019 conditional hindcast comparison of **Baseline + npieff** against **Humidity seasonal without npieff**. No owner acceptance declaration was found. Keep this variant distinct from older paired-npieff or other training windows.

Let `B` be `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019`. Baseline collection needs the following exact sets:

- `B/hokkaido/baseline/fit_cache.rds`, `B/hokkaido/baseline/result.rds`, `B/hokkaido/baseline/latest_attempt.rds`.
- `B/shikoku/baseline/fit_cache.rds`, `B/shikoku/baseline/result.rds`, `B/shikoku/baseline/latest_attempt.rds`.
- `B/ryukyu/baseline/fit_cache.rds`, `B/ryukyu/baseline/result.rds`, `B/ryukyu/baseline/latest_attempt.rds`.

No-npieff humidity collection uses these six intact files retained in Package A:

- `analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds` and sibling `humidity_seasonal_result.rds`.
- `analysis/humidity_seasonal_shikoku/humidity_seasonal_fit_cache.rds` and sibling `humidity_seasonal_result.rds`.
- `analysis/humidity_seasonal_ryukyu/humidity_seasonal_fit_cache.rds` and sibling `humidity_seasonal_result.rds`.

The optional sibling `latest_attempt_no_npieff.rds` is read if present; those filenames were absent in the initial inventory. All six model caches together occupy 177092683 bytes. Package A retains the three no-npieff humidity caches and their corresponding results (58022021 bytes combined), each below 20 MB. The three baseline cache/result sets remain archival dependencies; the Hokkaido baseline cache alone is 103208970 bytes and exceeds the requested 100 MB ordinary-tracking threshold.

`npieff_support/npieff_workflow.R` records schema, island, model, model_md5, standata, sampling_settings, draw_array, draws and diagnostics in baseline caches. `np_predict` requires the joint draw matrix containing 52 beta coefficients, the full npieff path through N_total, terminal S/I/R at N_train, rho/tau/omega/phi, plus standata, forecast weeks and X_test. Humidity prediction similarly uses beta/theta, terminal S/I/R, rho/tau/omega/phi, standata and X_test. Preserve full-precision joint arrays; marginal summaries cannot substitute for them.

The collectors check stored model and fit/result file MD5 links. Keep object bytes intact and do not substitute smaller backup versions. Full comparison plotting invokes replay of the saved humidity draws and also reads the archived baseline fit/result sets; it therefore requires Package B dependencies. Humidity-only replay and peak analysis can use the fitted objects included in Package A, subject to the documented static compatibility limitations. The separate saved-CSV metric plot workflow is described in REPRODUCIBILITY.md.

## Original manuscript fitted simulations

Original simulation scripts load `stanfit_sirs/stanfit_sirs_<island>_npi.rda`, extract posterior summaries and simulate from defined state grids. Three such fits are supplied (Honshu, Kyushu, Ryukyu); the required Hokkaido and Shikoku files are absent. The original alternative npi2 fit dependencies are also incomplete as recorded in PROVENANCE.md. No arbitrary parameters or refits replace missing files.

The Honshu trajectories below preserve existing outcomes but do not contain a documented portable frozen-state export sufficient to reproduce all original simulation workflows. `simulate_sirs_honshu_red.R` additionally has an unresolved `honshu_data$beta_red` dependency identified by static inspection.

## Frozen humidity objects retained in primary delivery

| Original relative path | Bytes |
|---|---:|
| `analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds` | 15560213 |
| `analysis/humidity_seasonal_hokkaido/humidity_seasonal_result.rds` | 13235724 |
| `analysis/humidity_seasonal_shikoku/humidity_seasonal_fit_cache.rds` | 7827491 |
| `analysis/humidity_seasonal_shikoku/humidity_seasonal_result.rds` | 6809064 |
| `analysis/humidity_seasonal_ryukyu/humidity_seasonal_fit_cache.rds` | 7828809 |
| `analysis/humidity_seasonal_ryukyu/humidity_seasonal_result.rds` | 6760720 |

## Compact result artifacts retained

This list identifies existing result files no larger than 1000000 bytes and the two existing Honshu trajectory objects. Their scientific status is preserved as supplied; see PROVENANCE.md for producer/consumer and ambiguity mapping.

| Original relative path | Bytes |
|---|---:|
| `analysis/baseline_humidity_seasonal_three_islands/hokkaido/baseline_result.rds` | 41339 |
| `analysis/baseline_humidity_seasonal_three_islands/shikoku/baseline_result.rds` | 41183 |
| `analysis/five_island_tuning_run1/comparison_results.rds` | 40849 |
| `analysis/five_island_tuning_run1/results/holdout__Hokkaido__existing_seasonal__v001.rds` | 35434 |
| `analysis/five_island_tuning_run1/results/holdout__Honshu__existing_seasonal__v001.rds` | 36352 |
| `analysis/five_island_tuning_run1/results/holdout__Kyushu__existing_seasonal__v001.rds` | 35704 |
| `analysis/five_island_tuning_run1/results/holdout__Ryukyu__existing_seasonal__v001.rds` | 35630 |
| `analysis/five_island_tuning_run1/results/holdout__Shikoku__existing_seasonal__v001.rds` | 35709 |
| `analysis/five_island_tuning_run1/results/reconstruction__Hokkaido__baseline__v001.rds` | 160143 |
| `analysis/five_island_tuning_run1/results/reconstruction__Honshu__baseline__v001.rds` | 163191 |
| `analysis/five_island_tuning_run1/results/reconstruction__Shikoku__baseline__v001.rds` | 162741 |
| `analysis/humidity_tuning_run1/comparison_results.rds` | 552472 |
| `analysis/humidity_tuning_run1/results/humidity_asymmetric.rds` | 219949 |
| `analysis/humidity_tuning_run1/results/humidity_quadratic.rds` | 222584 |
| `analysis/humidity_tuning_run1/results/humidity_seasonal.rds` | 230957 |
| `analysis/humidity_tuning_run1/results/original_baseline.rds` | 220423 |
| `analysis/humidity_tuning_run1/results/original_humidity.rds` | 218775 |
| `simulate_sirs/simulate_sirs_honshu_red.rda` | 178200 |
| `simulate_sirs/simulate_sirs_honshu.rda` | 178383 |

## Fitted and checkpoint files preserved as archival dependencies

Except for the six primary humidity objects listed above, fitted/checkpoint objects and large saved results remain in the archival delivery inventory. Historical versions remain distinct. A fit-cache suffix does not make a saved posterior disposable. The table below lists fitted/checkpoint objects; the archive manifest additionally includes large result files and exact checksums. Original workspaces remain untouched.

| Original relative path | Bytes |
|---|---:|
| `analysis/baseline_humidity_seasonal_three_islands/hokkaido/baseline_fit_cache.rds` | 14810975 |
| `analysis/baseline_humidity_seasonal_three_islands/ryukyu/baseline_fit_cache.rds` | 14890096 |
| `analysis/baseline_humidity_seasonal_three_islands/shikoku/baseline_fit_cache.rds` | 14896642 |
| `analysis/five_island_tuning_run1/fits/Hokkaido_baseline_2019__55dd4adf92b4c57b45ca405d8faf9a97/attempt_001.rds` | 39692571 |
| `analysis/five_island_tuning_run1/fits/Hokkaido_baseline_2019__55dd4adf92b4c57b45ca405d8faf9a97/attempt_002.rds` | 39706752 |
| `analysis/five_island_tuning_run1/fits/Honshu_baseline_2019__8c636ab2d0804305066d4fe67577a74c/attempt_001.rds` | 19903504 |
| `analysis/five_island_tuning_run1/fits/Honshu_baseline_2019__8c636ab2d0804305066d4fe67577a74c/attempt_002.rds` | 39693233 |
| `analysis/five_island_tuning_run1/fits/Kyushu_baseline_2019__39c5fb002ef4c4e44f3af700a7f5edc6/attempt_001.rds` | 19901927 |
| `analysis/five_island_tuning_run1/fits/Shikoku_baseline_2019__360d4b172783cde93ad364f04fd0cb53/attempt_001.rds` | 19952066 |
| `analysis/honshu_humidity_checkpoint.rda` | 26193444 |
| `analysis/humidity_seasonal_hokkaido_train2013_2015/humidity_seasonal_fit_cache.rds` | 15581862 |
| `analysis/humidity_seasonal_hokkaido_train2013_2016/humidity_seasonal_fit_cache.rds` | 38784595 |
| `analysis/humidity_seasonal_hokkaido_train2013_2016/humidity_seasonal_fit_cache.rds.before_refit_20261112_122720` | 15528307 |
| `analysis/humidity_seasonal_hokkaido_train2013_2016/humidity_seasonal_fit_cache.rds.before_refit_20261112_125936` | 38769506 |

| `analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds.before_refit_20261110_231900` | 7808608 |
| `analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds.before_refit_20261110_232958` | 15568498 |
| `analysis/humidity_seasonal_honshu_train2013_2015/humidity_seasonal_fit_cache.rds` | 7811634 |
| `analysis/humidity_seasonal_honshu_train2013_2015/humidity_seasonal_fit_cache.rds.before_refit_20261113_011305` | 7811519 |
| `analysis/humidity_seasonal_honshu_train2013_2016/humidity_seasonal_fit_cache.rds` | 7793766 |
| `analysis/humidity_seasonal_kyushu_train2013_2015/humidity_seasonal_fit_cache.rds` | 7810528 |
| `analysis/humidity_seasonal_kyushu_train2013_2016/humidity_seasonal_fit_cache.rds` | 7812088 |
| `analysis/humidity_seasonal_kyushu/humidity_seasonal_fit_cache.rds` | 23391105 |
| `analysis/humidity_seasonal_kyushu/humidity_seasonal_fit_cache.rds.before_refit_20261112_222847` | 7818011 |
| `analysis/humidity_seasonal_kyushu/humidity_seasonal_fit_cache.rds.before_refit_20261112_235800` | 7821264 |
| `analysis/humidity_seasonal_kyushu/humidity_seasonal_fit_cache.rds.before_refit_20261113_001810` | 7845282 |
| `analysis/humidity_seasonal_ryukyu_train2013_2015/humidity_seasonal_fit_cache.rds` | 7840187 |
| `analysis/humidity_seasonal_ryukyu_train2013_2016/humidity_seasonal_fit_cache.rds` | 7831237 |

| `analysis/humidity_seasonal_shikoku_train2013_2015/humidity_seasonal_fit_cache.rds` | 7849838 |

| `analysis/humidity_seasonal_standalone/humidity_seasonal_fit_cache.rds` | 7787586 |
| `analysis/humidity_seasonal_standalone/humidity_seasonal_fit_cache.rds.before_refit_20261112_222617` | 7791712 |
| `analysis/humidity_seasonal_standalone/humidity_seasonal_fit_cache.rds.before_refit_20261112_235240` | 7789780 |
| `analysis/humidity_tuning_run1/fit_cache/humidity_asymmetric.rds` | 7674238 |
| `analysis/humidity_tuning_run1/fit_cache/humidity_quadratic.rds` | 7679085 |
| `analysis/humidity_tuning_run1/fit_cache/humidity_seasonal.rds` | 7791576 |
| `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/hokkaido/baseline/fit_cache.rds` | 103208970 |
| `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/hokkaido/baseline/fit_cache.rds.before_write_20261113_144810` | 51655213 |
| `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/hokkaido/baseline/result.rds.before_write_20261113_144824` | 47700 |
| `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/ryukyu/baseline/fit_cache.rds` | 27392125 |
| `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/shikoku/baseline/fit_cache.rds` | 15275075 |
| `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/shikoku/humidity_seasonal/fit_cache.rds` | 15996825 |
| `analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/shikoku/humidity_seasonal/result.rds.before_write_20261113_135650` | 12882829 |
| `stanfit_sirs/honshu_baseline_train2013_2017.rda` | 129627703 |
| `stanfit_sirs/honshu_humidity_train2013_2017.rda` | 129977239 |
| `stanfit_sirs/stanfit_sirs_honshu_npi.rda` | 323444516 |
| `stanfit_sirs/stanfit_sirs_kyushu_npi.rda` | 323576439 |
| `stanfit_sirs/stanfit_sirs_ryukyu_npi.rda` | 336241034 |

Unknown workspaces `perturbation/analysis/.RData` (940970573 bytes), `perturbation/stanfit_sirs/.RData` (663802351 bytes), and `analysis/.RData` (709299 bytes) are archival artifacts with contents not interpreted. They are not assumed redundant.

The `stanmodel/*.rds` objects are inferred from their location and source usage to be compiled Stan model caches with corresponding `.stan` files. They were not deserialized, so their full contents were not verified. Conservatively retain them as Tier B archival artifacts rather than disposable development files. Source Stan models and required R implementations remain in Package A.

## Owner decisions

Identify which training window and model comparison constitute the accepted contribution; verify frozen cache/source links and recorded diagnostics; approve any future compact state export; resolve unavailable original fitted files; and decide whether historical tuning alternatives belong in eventual public release. No runtime or numerical-equivalence claim is made by this handoff.


