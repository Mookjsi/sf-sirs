# Report evidence map

The exact input files are included; see `docs/data-access.md`. Commands write only to new directories. `scripts/verify_report.py` records claim values and pass/fail status; a saved result is Level 1 evidence, while `reproduce_h1.R` and `reproduce_frozen.R` execute the original SIRS recurrence at each perturbation.

| Report result | Producer / reviewer command | Output | Verification level and status |
|---|---|---|---|
| Table 3.1 | `Rscript scripts/phase29.R evaluate --out artifacts/phase29-new` | `fitting_metrics.csv` | Previous independent fit-checkpoint replay; two displayed RMSE discrepancies |
| Supplement B.1–B.3 | Same Phase 29 `evaluate` command | `figures/case_reconstruction.pdf`, `epidemic_timing.pdf`, `timing_shift.pdf`; `peak_timing_by_season.csv` | Previously regenerated from fitted inputs |
| Table 3.2 | `Rscript scripts/reproduce_h1.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --grid artifacts/m2-full-h1-001/peak_path_grid.csv --parents artifacts/m2-full-h1-001/h1_event_summary.csv --out reproduction/h1` | `h1_edges.csv` | Level 2, verified from new simulations |
| Figure 3.1a | H1 command, then `Rscript scripts/render_figure3.R --h1 reproduction/h1 --mechanisms reproduction/mechanisms --out reproduction/figure3` | `figure_3_1.pdf`, `figure_3_1_data.csv` | Level 2 data independently simulated and plotted |
| Table 3.3 and Figure 3.1b–c | `Rscript scripts/reproduce_frozen.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --mode mechanisms --out reproduction/mechanisms`, then Figure 3 rendering command above | `h2_*`, `h3_*`, `model_integrity.csv`, `figure_3_1.pdf` | Level 2 raw and blend values independently simulated and plotted |
| Table 4.1, Figures 4.1–4.2 | `Rscript scripts/reproduce_comparison.R --evaluation artifacts/phase26-m2-evaluate-001/result.rds --park-weekly artifacts/phase28-park-fit-003/weekly_fit.csv --out reproduction/comparison --history artifacts/phase28-compare-003` | `fit_metrics.csv`, `observed_fitted.pdf`, `accuracy_runtime_coordinates.pdf` | Level 2 accuracy and trajectory recomputation; Figure 4.2 includes historical Level 1 runtime records |
| Table 4.2 | Same direct comparison command | `peak_metrics.csv`, `movement_metrics.csv` | Level 2 recomputation verified |
| Table 4.3 | Same direct comparison command with `--history` | `runtime_metrics.csv`, `parameter_inventory.csv` | Historical runtime and coordinate evidence checked |
| Table 5.1, Figures 5.1–5.3 | Chapter 5 prepare → forecast → score, as in `docs/reproduction.md` | `metrics.csv`, `scored_weekly.csv`, yearly PNGs | Prior isolated Level 2 run matched saved results |

The expected report values, numerical tolerances, inputs, and statuses are in `reproducibility/claims.yaml`. See `docs/reproducibility-audit.md` for the complete claim and dependency audit.
