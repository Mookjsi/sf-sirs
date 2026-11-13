# Included reviewer inputs

The project owner confirmed redistribution permission for the selected report inputs and evidence in this release. Exact source files are tracked under `artifacts/` and `data/`; their SHA-256 hashes are in `docs/source-checksums.csv`. The reference report is `final_2.pdf`, SHA-256 `50d65804d4e7fdfa2b10783e5536525a2dc431d1904ac63fca3de326af957a7e`.

| Purpose | Included files |
|---|---|
| Frozen M2 and Phase 29 | The 12 exact inputs in `config/phase29.R`, including selected `phase26-m2-fit-001/result.rds`, its evaluation, M0/M1 checkpoints, and raw A/AS Fourier fits |
| H1 | `artifacts/m2-full-h1-001/peak_path_grid.csv`, `h1_event_summary.csv`, and `classification_audit.csv` |
| H2 and H3 reference | `artifacts/m2-full-h2-001/` seasonal summary and step audit; corrected `artifacts/m2-full-h3-002/` local summary and finite effects |
| Park comparison | `artifacts/phase28-park-fit-003/weekly_fit.csv` with pointwise posterior mean and median; `artifacts/phase28-compare-003/` metrics, parameter inventory, runtime records, and figures |
| Chapter 5 | Exact 2024 NESID and 2025/2026 JIHS releases under `data/forcast/`; population and cases under `data/processed/original-island-v1-traceable/`; frozen forecast reference and fit hash under `artifacts/forcast-2024-verified/` |

The full Park posterior objects (roughly 452 MB) and exploratory refitting chain are not needed for the fast verification path and are not included. Park code copied from `parksw3/perturbation` remains pinned to commit `475c7dcdadc12987033d332b101f025625074301` and retains its MIT license and attribution in `references/parksw3/LICENSE`.

Saved RDS metadata records historical paths and earlier phase hashes for provenance. The direct H1–H3 and Park commands use repository-relative input paths and do not dereference those historical paths. New output directories remain ignored by Git; create fresh names rather than overwriting tracked reference artifacts.
