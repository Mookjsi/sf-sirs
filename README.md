# SF-SIRS report reproduction code

Selected R and Python code paths for the 48-page SF-SIRS report. Source observations, fitted RDS files, posterior draws, and the report PDF are excluded because their redistribution status was not established. Supply inputs separately at the documented paths. The [reproduction guide](docs/reproduction.md) maps report tables and figures to their producers, inputs, and verification status.

The source report has SHA-256 `50d65804d4e7fdfa2b10783e5536525a2dc431d1904ac63fca3de326af957a7e`. [Source checksums](docs/source-checksums.csv) identify selected files without publishing their contents.

Use R 4.5.3; versions are in [requirements-R.txt](requirements-R.txt). Run from this repository root into fresh output directories. With the saved Phase 22/26/27 inputs at paths in `config/phase29.R`, the retrospective comparison is:

```sh
Rscript scripts/phase29.R evaluate --out artifacts/phase29-evaluate-new
```

`phase29.R report` copies previously generated figures; `evaluate` calculates them. See the guide for H1–H3, Park, and the fixed 2023-W52 forecast commands. The Phase 30 multiple-origin experiment is not Chapter 5's execution path.

Park's original Stan model and caller/preprocessing references come from [`parksw3/perturbation`](https://github.com/parksw3/perturbation) at commit `475c7dcdadc12987033d332b101f025625074301`. They retain the [MIT license](references/parksw3/LICENSE) and copyright notice. See [third-party notices](THIRD_PARTY_NOTICES.md).
