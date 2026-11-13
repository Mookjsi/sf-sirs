# SF-SIRS report reproduction code

Selected R and Python code paths plus the frozen inputs and numerical evidence for the 48-page [SF-SIRS report](final_2.pdf). The project owner confirmed redistribution permission for the selected inputs in this release. The large Park posterior objects and exploratory refitting chain are excluded. The [reproduction guide](docs/reproduction.md) maps report tables and figures to their producers, inputs, and verification status.

The source report has SHA-256 `50d65804d4e7fdfa2b10783e5536525a2dc431d1904ac63fca3de326af957a7e`. [Source checksums](docs/source-checksums.csv) identify selected files without publishing their contents.

Use R 4.5.3; versions are in [requirements-R.txt](requirements-R.txt). Run from this repository root into fresh output directories. The retrospective comparison is:

```sh
Rscript scripts/phase29.R evaluate --out artifacts/phase29-evaluate-new
```

`phase29.R report` copies previously generated figures; `evaluate` calculates them. See the guide for H1–H3, Park, and the fixed 2023-W52 forecast commands. The Phase 30 multiple-origin experiment is not Chapter 5's execution path.

Park's original Stan model and caller/preprocessing references come from [`parksw3/perturbation`](https://github.com/parksw3/perturbation) at commit `475c7dcdadc12987033d332b101f025625074301`. They retain the [MIT license](references/parksw3/LICENSE) and copyright notice. See [third-party notices](THIRD_PARTY_NOTICES.md).

## Direct frozen M2 verification

A reviewer can run H1–H3 without the exploratory Phase 7C–27 fit checkpoints. All inputs below are included in a clean clone:

```sh
Rscript scripts/reproduce_frozen.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --mode mechanisms --out reproduction/mechanisms
Rscript scripts/reproduce_h1.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --grid artifacts/m2-full-h1-001/peak_path_grid.csv --parents artifacts/m2-full-h1-001/h1_event_summary.csv --out reproduction/h1
Rscript scripts/render_figure3.R --h1 reproduction/h1 --mechanisms reproduction/mechanisms --out reproduction/figure3
Rscript scripts/reproduce_comparison.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --park-weekly artifacts/phase28-park-fit-003/weekly_fit.csv --out reproduction/comparison --history artifacts/phase28-compare-003
python3 scripts/verify_report.py --mechanisms reproduction/mechanisms --h1 reproduction/h1 --phase29 artifacts/phase29-evaluate-003 --phase28 reproduction/comparison --phase28-history artifacts/phase28-compare-003 --forecast artifacts/forcast-2026-extended --out reproduction/verification
```

The first command checks all 573 stored S/I/R and expected-case values by rerunning the SIRS recurrence, then recomputes raw and boundary-blended H2/H3 timing experiments. The second reruns every H1 grid simulation and reclassifies parent events. The comparison command recalculates Park metrics and timing from the saved weekly posterior mean and median without fresh MCMC. `--phase28-history` reads the original runtime and parameter inventory as Level 1 records. The Python verifier writes `verification_results.csv`, `verification_summary.md`, and `manifest.json`; missing chapter inputs appear as `BLOCKED`. Pass fresh `--phase29` and `--forecast` output directories to verify their tables as well. Each command requires a new output directory and never overwrites source artifacts.

See the [claim registry](reproducibility/claims.yaml), [audit](docs/reproducibility-audit.md), and [report evidence map](docs/report-evidence-map.md). The direct H1–H3 path was verified against the archived results. Table 3.1 has two displayed RMSE discrepancies of 0.01 relative to its saved Phase 29 output; the verifier reports them as failures.
