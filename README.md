# SF-SIRS report reproduction code

Selected R and Python code paths for the 48-page SF-SIRS report. Source observations, fitted RDS files, posterior draws, and the report PDF are excluded because their redistribution status was not established. Supply inputs separately at the documented paths. The [reproduction guide](docs/reproduction.md) maps report tables and figures to their producers, inputs, and verification status.

The source report has SHA-256 `50d65804d4e7fdfa2b10783e5536525a2dc431d1904ac63fca3de326af957a7e`. [Source checksums](docs/source-checksums.csv) identify selected files without publishing their contents.

Use R 4.5.3; versions are in [requirements-R.txt](requirements-R.txt). Run from this repository root into fresh output directories. With the saved Phase 22/26/27 inputs at paths in `config/phase29.R`, the retrospective comparison is:

```sh
Rscript scripts/phase29.R evaluate --out artifacts/phase29-evaluate-new
```

`phase29.R report` copies previously generated figures; `evaluate` calculates them. See the guide for H1–H3, Park, and the fixed 2023-W52 forecast commands. The Phase 30 multiple-origin experiment is not Chapter 5's execution path.

Park's original Stan model and caller/preprocessing references come from [`parksw3/perturbation`](https://github.com/parksw3/perturbation) at commit `475c7dcdadc12987033d332b101f025625074301`. They retain the [MIT license](references/parksw3/LICENSE) and copyright notice. See [third-party notices](THIRD_PARTY_NOTICES.md).

## Direct frozen M2 verification

An authorized reviewer can run H1–H3 without the exploratory Phase 7C–27 fit checkpoints. Place the three exact files and hashes listed in [data access](docs/data-access.md) under an external `INPUT/` directory, then run from a clean clone:

```sh
Rscript scripts/reproduce_frozen.R --evaluation INPUT/m2-evaluation.rds --mode mechanisms --out reproduction/mechanisms
Rscript scripts/reproduce_h1.R --evaluation INPUT/m2-evaluation.rds --grid INPUT/peak_path_grid.csv --parents INPUT/h1_event_summary.csv --out reproduction/h1
python3 scripts/verify_report.py --mechanisms reproduction/mechanisms --h1 reproduction/h1 --out reproduction/verification
```

The first command checks all 573 stored S/I/R and expected-case values by rerunning the SIRS recurrence, then recomputes raw and boundary-blended H2/H3 timing experiments. The second reruns every H1 grid simulation and reclassifies parent events. The Python verifier writes `verification_results.csv`, `verification_summary.md`, and `manifest.json`; missing chapter inputs appear as `BLOCKED`. Pass fresh `--phase29`, `--phase28`, and `--forecast` output directories to verify their tables as well. Each command requires a new output directory and never overwrites source artifacts.

See the [claim registry](reproducibility/claims.yaml), [audit](docs/reproducibility-audit.md), and [report evidence map](docs/report-evidence-map.md). The direct H1–H3 path was verified against the archived results. Table 3.1 has two displayed RMSE discrepancies of 0.01 relative to its saved Phase 29 output; the verifier reports them as failures. The public clone does not contain the fitted state or observations because redistribution permission has not been established.
