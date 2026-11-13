# Reproduction map and limits

## Source and environment

The source is the 48-page `final_2.pdf`, SHA-256 `50d65804d4e7fdfa2b10783e5536525a2dc431d1904ac63fca3de326af957a7e`. No TeX source was found in the source working tree. [Source checksums](source-checksums.csv) identify exact local inputs and completed results; those files are excluded from this public repository. Code was copied as individual files, without source Git objects or configuration.

The source used R 4.5.3, with package versions in `requirements-R.txt`. Phase 22 seed: `20261111`; Phase 24/26 seed: `20261112`. Park used seed `102`, four chains, 4,000 iterations, `adapt_delta=0.9`, and tree depth 15. Use fresh output directories under ignored `artifacts/`.

| Report | Completed source evidence | Producer and needed inputs |
|---|---|---|
| Ch. 2, App. A | `phase26-m2-fit-001/result.rds` | `scripts/phase26.R`, `config/phase26.R`, common `R/` modules; earlier Phase 22/24/25 checkpoints |
| Table 3.1, Figs. B.1–B.3 | `phase29-evaluate-003/{fitting_metrics,timing_shift,peak_timing_by_season}.csv` and `figures/{case_reconstruction,epidemic_timing,timing_shift}.pdf` | `scripts/phase29.R evaluate`, `R/phase29.R`, `config/phase29.R`; its 12 listed Phase 22/26/27 inputs |
| Tables 3.2–3.3, Fig. 3.1 | `m2-full-h1-001`, `m2-full-h2-001`, corrected H3 `m2-full-h3-002`, `m2-full-report-001` | `scripts/m2_hypotheses.R`; Phase 27 replay/mechanism/sensitivity and M2 checkpoints |
| Tables 4.1–4.3, Figs. 4.1–4.2 | `phase28-compare-003/{fit_metrics,peak_metrics,movement_metrics,runtime_metrics,parameter_inventory}.csv`, `observed_fitted.pdf`, `accuracy_runtime_coordinates.pdf` | `scripts/phase28.R prepare → verify → fit M2/PARK → compare`; panel, previous fits, Park Stan and posterior results |
| Table 5.1, Figs. 5.1–5.3 | `forcast-2024-verified`, `forcast-2026-extended/{metrics,scored_weekly,forecast_weekly}.csv`, `forecast_2024.png` through `forecast_2026.png` | `scripts/forcast/{prepare.py,forecast.R,score.R}`; 2023-W52 M2 fit and exact JIHS releases |

All paths in the table are under `artifacts/` in the source run and are intentionally absent here. Phase 29 `report` copies an evaluation's figures; `evaluate` creates them. The common panel has 573 rows and 572 likelihood rows. Seasons span W27 to the next W26; earliest ties are chosen and tie/boundary flags retained. The source M2 MAE/RMSE is 264.74032641234/533.792847544409. Against the observed 13-week advance, M0/M1/M2/A_M2/AS_M2 reproduce 5/9/11/9/8 weeks. Fourier inputs are the raw A and AS fit paths in `config/phase29.R`.

The H1 final run uses β log multipliers `[-log(1.2), log(1.2)]`, wider than the historical `config/phase27.R` default. Its raw parent grid has 72 resolved displacements among 80 timing events; blend has 69 among 77. Mixed and unresolved intervals remain separate. H2 uses the median eligible one of five steps: 7/8 seasons have `|Eβ|>|ES|`; 2021/22 has ES 11.101892 and Eβ 7.020422. H3 holds initial S and β fixed, balances R, and compares I/1.2 with 1.2I. Use corrected H3 `-002`, not the earlier reporting pass. The original experiments were not rerun here.

Park MAE/WAPE use posterior median expected cases; RMSE and the primary Table 4.2 trajectory use posterior mean. Original Park MAE/RMSE is 419.719249291797/744.741299455345. Both M2 and Park advance the peak 11 weeks. Free coordinates are 64/269. Recorded single-run wall times are 608.23056435585/2153.27476191521 seconds, not timing targets for a new fit. The complete Park run is `phase28-park-fit-003`; `-001` is incomplete. Its approximately 452 MB posterior draw file is excluded, while the fit code retains posterior trajectory summarization.

## Commands with supplied inputs

Place the 12 files listed in `config/phase29.R` at those paths, checking their hashes against `source-checksums.csv`. Then run `Rscript scripts/phase29.R evaluate --out artifacts/phase29-evaluate-new`. Earlier fit initialization needs saved joint-fit coordinates; the bare clone cannot refit. The Phase 28/H1–H3 code likewise expects earlier checkpoints named in `config/phase22.R`, `phase24.R`, `phase25.R`, `phase26.R`, and `phase27.R`. These are explicit external inputs; private repository access is not a runtime assumption.

```sh
Rscript scripts/phase28.R prepare --out artifacts/phase28-prepare-new
Rscript scripts/phase28.R verify --out artifacts/phase28-verify-new --prepare artifacts/phase28-prepare-new/result.rds
Rscript scripts/phase28.R fit --model M2 --out artifacts/phase28-m2-fit-new --prepare artifacts/phase28-prepare-new/result.rds --verify artifacts/phase28-verify-new/result.rds
Rscript scripts/phase28.R fit --model PARK --out artifacts/phase28-park-fit-new --prepare artifacts/phase28-prepare-new/result.rds --verify artifacts/phase28-verify-new/result.rds
Rscript scripts/phase28.R compare --out artifacts/phase28-compare-new --prepare artifacts/phase28-prepare-new/result.rds --verify artifacts/phase28-verify-new/result.rds --m2 artifacts/phase28-m2-fit-new/result.rds --park artifacts/phase28-park-fit-new/result.rds
```

For H1–H3 run `plan`, `h1-full`, `summarize-h2`, `state-memory`, then `report` with `scripts/m2_hypotheses.R`. Each stage takes `--out NEW_DIR --replay DIR --mechanism DIR --sensitivity DIR`; `report` takes `--inventory`, `--h1`, `--h2`, `--state`. The source runs recorded `STORAGE_COMPLETED` but their checkpoints remain external.

For Chapter 5, obtain the finalized [2024 NESID annual Table 8-1](https://id-info.jihs.go.jp/en/surveillance/idwr/rapid/index.html) as `Syu_08_1.xlsx`, plus JIHS IDWR 2025 W52 and 2026 W44 year-to-date CSVs. Check the exact filenames and SHA-256 values in `source-checksums.csv`. Supply an authorized 47-prefecture population mapping, 2013–2023 case history, Phase 26 M2 fit, original 66-week forecast, and its fit hash.

```sh
python3 scripts/forcast/prepare.py --annual INPUT/Syu_08_1.xlsx --weekly INPUT/2025-52-teiten-tougai.csv --weekly INPUT/2026-44-teiten-tougai.csv --population INPUT/japan_pop_prefecture.csv --cases INPUT/japan_rsv_prefecture.csv --out artifacts/forcast-new
Rscript scripts/forcast/forecast.R INPUT/phase26-m2-result.rds artifacts/forcast-new --reference INPUT/forecast-66-weeks.csv --hash INPUT/forecast_fit_hash.txt --population INPUT/japan_pop_prefecture.csv
Rscript scripts/forcast/score.R artifacts/forcast-new
```

The forecast fixes the 2023-W52 nonseasonal multiplier, repeats seasonal β, propagates S/I/R without year resets, and uses origin S+I+R as the denominator. 2026 W53 uses W52 β. The first 66 original predictions are checked within 1e-6. Scoring requires all 34 Honshu prefectures per week and splits at the 2025 W15 surveillance change. The original four M2/naïve MAE pairs are 1157.5/1132.0, 712.7/1470.6, 526.6/2629.0, and 643.3/2187.3. There are 148 observed weeks, ending at 2026 W44. No pediatric sentinel denominators were supplied for per-sentinel analysis.

## Verification and release boundary

The public checkout alone cannot recompute report values because fitted states and observations are absent. In an isolated copy of this package with separately supplied source inputs, Phase 29 `evaluate` regenerated all seven CSV files with every row exactly equal to the completed source run; it also generated the three Appendix B figure PDFs. The Chapter 5 prepare/forecast/score sequence regenerated the 157-week forecast, 148-week scored panel, four metric rows, denominator sensitivity and other CSV outputs with every row exactly equal to the completed source run; it also generated all three yearly figures. The comparison ignored source manifest absolute paths. Refitting, H1/H3 simulation, and Park MCMC were not rerun. Do not substitute Phase 30 forecasts or a newer JIHS release.

The [upstream Park repository](https://github.com/parksw3/perturbation) publicly distributes code and data, but reuse rights for the competition training file, JIHS snapshots, and prepared/result files were not established for this release. None are committed. Generated RDS, CSV, logs, and figures can contain observed data or absolute paths; output directories are ignored.
