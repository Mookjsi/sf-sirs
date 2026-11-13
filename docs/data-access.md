# Reviewer inputs and redistribution boundary

The public code can be cloned independently. The exact numerical inputs are supplied separately by a reviewer with permission to use them. `docs/source-checksums.csv` lists SHA-256 hashes; the direct scripts reject an unexpected M2 evaluation or H1 design. No runtime path points to the source author's archive.

| File in reviewer `INPUT/` | Source provenance | Status |
|---|---|---|
| `m2-evaluation.rds` | `artifacts/phase26-m2-evaluate-001/result.rds`; selected M2 573-week trajectory and restart states; SHA-256 `45c5e155d900da17637dbb26a8d562af2315c2d863a78df05679b10728c65e18` | Redistribution permission for fitted state and embedded observations not established; keep outside public Git |
| `peak_path_grid.csv` | `artifacts/m2-full-h1-001/peak_path_grid.csv`; H1 declared grid; SHA-256 `7a8913d47d1a98be409981709d23e5530f7130f83079707280c37be8e566c6dd` | Derived result redistribution unconfirmed; keep outside public Git |
| `h1_event_summary.csv` | `artifacts/m2-full-h1-001/h1_event_summary.csv`; only parent interval bounds are used by the new H1 simulation; SHA-256 `5fc580cce9e481694eb1f9d2fe30dc30081df2dd004b687e824e48ef1f68880e` | Derived result redistribution unconfirmed; keep outside public Git |
| Phase 29 fit/checkpoint files | Exact 12 paths in `config/phase29.R` and hashes in `source-checksums.csv` | Fit states and source observations are not cleared for publication |
| Phase 28 panel, M2 fit, Park result | Exact paths in `config/phase28.R`; Park `phase28-park-fit-003/weekly_fit.csv` contains pointwise mean and median | Posterior summaries and panel are not cleared for publication; full posterior object is approximately 452 MB and unnecessary for the fast compare path |
| Chapter 5 forecast inputs | Final 2024 NESID annual Table 8-1; JIHS 2025-W52 and 2026-W44 CSV releases; M2 fit; original 66-week forecast, population and historical cases | Source release is publicly viewable, but redistribution authorization for these snapshots and derived panels was not established; supply exact hashed files separately |

The selected Park implementation is pinned to `parksw3/perturbation` commit `475c7dcdadc12987033d332b101f025625074301`. Its copied code retains the MIT license in `references/parksw3/LICENSE`; this license does not itself establish rights to unrelated competition data, posterior draws, or surveillance releases. For external input rights, obtain the applicable source terms or authorization before making a public companion package.

Keep `INPUT/` outside the repository or in ignored `data/`. Generated `artifacts/` and `reproduction/` directories are ignored. A private reviewer package can contain the exact files above after rights and access are confirmed; include their SHA-256 values and the package recipient scope.
