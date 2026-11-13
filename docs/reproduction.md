# Reproducing the SF-SIRS report

Run commands from the repository root. The reference is [final_2.pdf](../final_2.pdf), SHA-256 `50d65804d4e7fdfa2b10783e5536525a2dc431d1904ac63fca3de326af957a7e`. Use fresh output directories under the ignored `artifacts/` path.

## Environment and inputs

Use R 4.5.3 with the packages in [requirements-R.txt](../requirements-R.txt), plus Python 3. The forecast preparation and numerical verification scripts use the Python standard library. [Source checksums](source-checksums.csv) record the included inputs. The model comparison reads the twelve frozen inputs in `config/phase29.R`. H1–H3 use `artifacts/phase26-m2-evaluate-001/result.rds`; H1 also uses the grid and parent events in `artifacts/m2-full-h1-001/`. Park comparison reads `artifacts/phase28-park-fit-003/weekly_fit.csv`. Chapter 5 reads the NESID and JIHS releases in `data/forcast/`, the Honshu case and population files in `data/processed/original-island-v1-traceable/`, the frozen M2 fit, and the 66-week forecast reference.

## Model comparison and H1–H3

```sh
Rscript scripts/phase29.R evaluate --out artifacts/phase29-evaluate-new
Rscript scripts/reproduce_frozen.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --mode mechanisms --out artifacts/mechanisms-new
Rscript scripts/reproduce_h1.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --grid artifacts/m2-full-h1-001/peak_path_grid.csv --parents artifacts/m2-full-h1-001/h1_event_summary.csv --out artifacts/h1-new
Rscript scripts/render_figure3.R --h1 artifacts/h1-new --mechanisms artifacts/mechanisms-new --out artifacts/figure3-new
```

Phase 29 writes seven CSVs and three figures under its output directory. Its 572 likelihood rows give M2 MAE `264.74032641234` and RMSE `533.792847544409`; M0/M1/M2/A_M2/AS_M2 reproduce peak advances of 5/9/11/9/8 weeks against the observed 13-week advance. H1 writes parent classifications and edges: 72 of 80 raw and 69 of 77 blended events are pure displacements. H2 and H3 write seasonal sensitivities under `artifacts/mechanisms-new/`. In H2, transmission has the larger absolute sensitivity in seven of eight eligible seasons. In H3, eight eligible local slopes are negative and seven seasons have eligible finite contrasts. Figure 3.1 is written under `artifacts/figure3-new/`.

H1 uses transmission log multipliers `[-log(1.2), log(1.2)]`; mixed and unresolved intervals remain distinct. H2 compares the median eligible one of five steps. H3 varies initial infectious prevalence between `I/1.2` and `1.2I`, holds initial susceptibility and transmission fixed, and balances the recovered state. The report's H3 reference is `artifacts/m2-full-h3-002/`.

## Park comparison

```sh
Rscript scripts/reproduce_comparison.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --park-weekly artifacts/phase28-park-fit-003/weekly_fit.csv --out artifacts/comparison-new --history artifacts/phase28-compare-003
```

The output contains reconstruction, peak, and movement metrics and figures. Park MAE/WAPE use posterior median expected cases; RMSE and the primary timing trajectory use the posterior mean. Park MAE/RMSE are `419.719249291797`/`744.741299455345`; both models advance the peak by 11 weeks. Runtime and coordinate counts are stored observations from distinct inferential workflows, not benchmarks for this replay. Full Park MCMC requires its original sampling inputs and computation; this command uses stored weekly summaries.

## Forecast and scoring

```sh
python3 scripts/forcast/prepare.py --annual data/forcast/Syu_08_1.xlsx --weekly data/forcast/2025-52-teiten-tougai.csv --weekly data/forcast/2026-44-teiten-tougai.csv --population data/processed/original-island-v1-traceable/japan_pop_prefecture.csv --cases data/processed/original-island-v1-traceable/japan_rsv_prefecture.csv --out artifacts/forecast-new
Rscript scripts/forcast/forecast.R artifacts/phase26-m2-fit-001/result.rds artifacts/forecast-new --reference artifacts/forcast-2024-verified/forecast_weekly.csv --hash artifacts/forcast-2024-verified/forecast_fit_hash.txt --population data/processed/original-island-v1-traceable/japan_pop_prefecture.csv
Rscript scripts/forcast/score.R artifacts/forecast-new
```

The output has 157 forecast weeks, 148 scored weeks through 2026 W44, four metric rows, and yearly figures. The four M2/naïve MAE pairs are `1157.5/1132.0`, `712.7/1470.6`, `526.6/2629.0`, and `643.3/2187.3`. The forecast fixes the 2023 W52 origin and compares its first 66 predictions with the included reference within `1e-6`. Scoring requires all 34 Honshu prefectures and separates observations at the 2025 W15 surveillance change. Per-sentinel analysis needs pediatric sentinel denominators, which are unavailable here.

## Numerical verification and interpretation

```sh
Rscript tests/test_invariants.R
python3 scripts/verify_report.py --mechanisms artifacts/mechanisms-new --h1 artifacts/h1-new --phase29 artifacts/phase29-evaluate-new --phase28 artifacts/comparison-new --phase28-history artifacts/phase28-compare-003 --forecast artifacts/forecast-new --out artifacts/verification-new
```

The verifier writes `verification_results.csv`, `verification_summary.md`, and `manifest.json`; missing inputs are marked `BLOCKED`. It compares calculated metrics with the report at the stated display precision.

These commands reproduce calculations from frozen fitted states, experiment designs, and stored Park posterior summaries. They do not perform full model refitting or Park MCMC. The fitted decomposition and controlled perturbations do not uniquely identify climatic, behavioral, reporting, or other historical causes of the observed shift.

## Exploratory climate analysis

With the four original CSV inputs under `data/climate/source/`, run `python3 scripts/reproduce_exploratory.py --input-dir data/climate/source --out artifacts/exploratory-climate-new`. Select one component using `--analysis NAME`. See the [module guide](../extensions/humidity/exploratory/README.md) and [figure provenance](exploratory-provenance.md). These descriptive plots do not reproduce Chapter 6 Bayesian model fitting; the remaining fitted-object gap is recorded in the provenance note.
