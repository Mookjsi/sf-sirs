# SF-SIRS: Spline-Fourier SIRS Modeling of RSV Transmission in Japan

## Overview

Japan's respiratory syncytial virus (RSV) epidemics shifted from winter toward summer. SF-SIRS is a finite-dimensional susceptible–infectious–recovered–susceptible transmission model with a spline and Fourier representation. Building on the mechanistic questions of Park et al. (2025), this study compares fitted epidemic trajectories and tests how peak timing responds to transmission, initial susceptibility, and initial infectious prevalence. See the [report](final_2.pdf) for methods and interpretation.

## Main findings

- In the selected M2 transmission experiments, peak shifts are predominantly classified as displacement rather than dominance exchange: 72 of 80 raw-grid events and 69 of 77 boundary-blended events are pure displacements.
- Transmission perturbations have larger peak-timing sensitivity than initial susceptibility in seven of eight eligible seasons; 2021/22 is the exception.
- Initial infectious prevalence provides an additional timing-response mechanism when initial susceptibility and transmission are held fixed.

These controlled model results do not identify the historical cause of Japan's RSV timing shift.

## Repository structure

| Path | Contents |
|---|---|
| `R/` | Model, estimation, analysis, and frozen-result reproduction functions |
| `config/` | Model and experiment configurations |
| `scripts/` | Executable R and Python workflows |
| `data/` | Included case, population, and surveillance inputs |
| `artifacts/` | Frozen fits, experiment designs, and reference outputs |
| `tests/` | Scientific invariant tests |
| `references/` | Attributed Park implementation and license |

## Requirements

Run from the repository root with R 4.5.3 and the packages in [requirements-R.txt](requirements-R.txt). Python 3 is required for forecast preparation and numerical verification; both scripts use the standard library. Use a fresh output directory for each command. The [reproduction guide](docs/reproduction.md) gives inputs, expected values, and limitations.

## Reproducing the results

```sh
Rscript scripts/phase29.R evaluate --out artifacts/phase29-evaluate-new
Rscript scripts/reproduce_frozen.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --mode mechanisms --out artifacts/mechanisms-new
Rscript scripts/reproduce_h1.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --grid artifacts/m2-full-h1-001/peak_path_grid.csv --parents artifacts/m2-full-h1-001/h1_event_summary.csv --out artifacts/h1-new
Rscript scripts/reproduce_comparison.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --park-weekly artifacts/phase28-park-fit-003/weekly_fit.csv --out artifacts/comparison-new --history artifacts/phase28-compare-003
python3 scripts/forcast/prepare.py --annual data/forcast/Syu_08_1.xlsx --weekly data/forcast/2025-52-teiten-tougai.csv --weekly data/forcast/2026-44-teiten-tougai.csv --population data/processed/original-island-v1-traceable/japan_pop_prefecture.csv --cases data/processed/original-island-v1-traceable/japan_rsv_prefecture.csv --out artifacts/forecast-new
Rscript scripts/forcast/forecast.R artifacts/phase26-m2-fit-001/result.rds artifacts/forecast-new --reference artifacts/forcast-2024-verified/forecast_weekly.csv --hash artifacts/forcast-2024-verified/forecast_fit_hash.txt --population data/processed/original-island-v1-traceable/japan_pop_prefecture.csv
Rscript scripts/forcast/score.R artifacts/forecast-new
python3 scripts/verify_report.py --mechanisms artifacts/mechanisms-new --h1 artifacts/h1-new --phase29 artifacts/phase29-evaluate-new --phase28 artifacts/comparison-new --phase28-history artifacts/phase28-compare-003 --forecast artifacts/forecast-new --out artifacts/verification-new
```

H2 and H3 are produced together by `reproduce_frozen.R`. The Park comparison uses stored posterior weekly summaries.

## Data and references

Included inputs and their SHA-256 hashes are listed in [source-checksums.csv](docs/source-checksums.csv). The original Park implementation is from [`parksw3/perturbation`](https://github.com/parksw3/perturbation), commit `475c7dcdadc12987033d332b101f025625074301`. Cite the included report and Park et al. (2025) when using these results.

## License and attribution

The Park source retains its [MIT license](references/parksw3/LICENSE) and copyright notice. [Third-party notices](THIRD_PARTY_NOTICES.md) describe other attribution and reuse boundaries; no whole-project license is asserted there.
