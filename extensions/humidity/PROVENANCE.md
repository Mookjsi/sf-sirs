# Result provenance and acceptance boundaries

No result was regenerated. Every retained scientific artifact has a source identifier, SHA-256 and retention tier. `RESULT_PROVENANCE.csv` provides per-file known producers or explicit unresolved status; `ARTIFACT_MANIFEST.csv` and `file_inventory.csv` give exact original/delivery/restore paths. Dynamic file references are not proven by basename matching: `docs/STATIC_REFERENCES.csv` is an index of text evidence, not an execution trace.

The original README and manuscript sources establish the original pipeline. Current comparison script labels and saved CSVs establish the most evident later candidate analysis, but no explicit acceptance record was found. No newest-file heuristic was used. Existing diagnostic/ready fields were preserved without endorsing them.

## Current analysis pipeline

1. Retained CSV inputs: population, RSV cases, processed humidity; original common preparation additionally reads processed temperature. Raw daily ERA5 inputs are absent.
2. Baseline Stan specification is `code/perturbation/stanmodel/sirs_npieff_holdout.stan`; no-npieff humidity specifications are embedded in the three current humidity scripts. Their bytes/model strings were preserved.
3. Frozen current humidity cache/result pairs are A. Baseline pairs are B; small attempt/hash metadata are A. `MODEL_ARTIFACTS.md` lists the full state interfaces and exact paths.
4. Collectors use saved inputs, terminal states, parameters and joint draws. Some perform deterministic trajectory replay and seed-202 predictive sampling before scoring.
5. Saved metric CSV can be plotted separately. Other original producers calculate numerical results before plotting.

`P` below is `code/perturbation`; `A` is `code/analysis`. Brace lists in paths denote separate concrete files, not literal filenames. Execution commands are documented for the owner only; use an isolated working copy and environment from `REPRODUCIBILITY.md`.

### mixed_npieff_comparison_2018_2019

- Supplied outputs: `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/comparison_metrics.csv`, `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/comparison_metrics_by_year.csv`, `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/comparison_diagnostics.csv`, `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/comparison_status.csv`, `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/prediction_plot_data.csv`, `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/baseline_npieff_vs_humidity_no_npieff_2018_2019.png`.
- Producer: `code/analysis/2013-2017 training/Baseline vs Humidity anomaly.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data/japan_rsv.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Helpers: `code/analysis/npieff_support/npieff_prediction_helpers.R`, `code/analysis/npieff_support/npieff_workflow.R`, `code/analysis/npieff_support/humidity_no_npieff_support.R`.
- Frozen dependency set: `code/perturbation/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/hokkaido/baseline/fit_cache.rds` (B), `code/perturbation/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/hokkaido/baseline/result.rds` (B), `code/perturbation/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/shikoku/baseline/fit_cache.rds` (B), `code/perturbation/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/shikoku/baseline/result.rds` (B), `code/perturbation/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/ryukyu/baseline/fit_cache.rds` (B), `code/perturbation/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019/ryukyu/baseline/result.rds` (B), `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_result.rds` (A).
- Availability: Baseline fit/result B required; humidity pairs and baseline attempt metadata A.
- Computation: Original producer requires archived fits/results and replays saved humidity draws (expected counts and seed-202 NB2 noise), then computes metrics/plots. Existing CSV and PNG inspectable without archive.
- Seed/settings: 202 for predictive observation sampling; fitted seeds/settings inside caches authoritative.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/Baseline vs Humidity anomaly.R"
```

### metrics_chart_2018_2019

- Supplied outputs: `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/forecast_metrics_bars_2018_2019.png`, `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/forecast_metrics_bars_2018_2019.pdf`.
- Producer: `code/analysis/2013-2017 training/plot_forecast_metrics.R`.
- Inputs: `code/analysis/2013-2017 training/plots_baseline_npieff_vs_humidity_no_npieff/comparison_metrics.csv`.
- Availability: Saved CSV and chart producer A; plotting packages/fonts required.
- Computation: Primary-only figure regeneration from saved numerical table; ggplot2, Malgun Gothic font and Cairo device. No fitted artifact required.
- Seed/settings: none.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/plot_forecast_metrics.R"
```

### humidity_three_island_plot_2018_2019

- Supplied outputs: `code/analysis/2013-2017 training/plots_humidity_seasonal_no_npieff/humidity_seasonal_no_npieff_three_islands_2018_2019.png`.
- Producer: `code/analysis/2013-2017 training/plot_hokkaido_shikoku_ryukyu_2018_2019.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data/japan_rsv.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Helpers: `code/analysis/npieff_support/npieff_workflow.R`, `code/analysis/npieff_support/humidity_no_npieff_support.R`.
- Frozen dependency set: `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_result.rds` (A).
- Availability: All current humidity fit/result inputs included in A; static dependency completeness, runtime not verified.
- Computation: No MCMC, but hn_collect replays saved trajectories and seed-202 observation noise; all current humidity pairs included in Package A.
- Seed/settings: 202.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/plot_hokkaido_shikoku_ryukyu_2018_2019.R"
```

### humidity_island_results_2018_2019

- Supplied outputs: `code/perturbation/analysis/humidity_seasonal_{hokkaido,shikoku,ryukyu}/metrics.csv`, `code/perturbation/analysis/humidity_seasonal_{hokkaido,shikoku,ryukyu}/peaks_2018_2019.csv`, `code/perturbation/analysis/humidity_seasonal_{hokkaido,shikoku,ryukyu}/prediction_2018_2019.csv`, `code/perturbation/analysis/humidity_seasonal_{hokkaido,shikoku,ryukyu}/humidity_seasonal_2018_2019.png`.
- Producer: `code/analysis/2013-2017 training/humidity_seasonal_{hokkaido,sikoku,ryukyu}.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data/japan_rsv.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Helpers: `code/perturbation/script/setup.R`, `code/perturbation/script/script_data.R`.
- Frozen dependency set: `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_result.rds` (A).
- Availability: All current humidity fit/result inputs included in A; static dependency completeness, runtime not verified.
- Computation: Reuses matching saved posterior cache; if absent performs costly fitting. Embedded Stan code preserved. Always reruns downstream predictions/evaluation.
- Seed/settings: Hokkaido fitting 303; Shikoku/Ryukyu fitting 101; predictive observations 202.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/humidity_seasonal_hokkaido.R"
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/humidity_seasonal_sikoku.R"
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/humidity_seasonal_ryukyu.R"
```

### humidity_peak_shift_three_island_2018_2019

- Supplied outputs: `code/perturbation/analysis/humidity_seasonal_three_islands_train2013_2017/peak_shift_summary_three_islands_2018_2019.csv`, `code/perturbation/analysis/humidity_seasonal_three_islands_train2013_2017/peak_shift_checks_three_islands_2018_2019.csv`, `code/perturbation/analysis/humidity_seasonal_three_islands_train2013_2017/peak_shift_distribution_three_islands_2018_2019.csv`, `code/perturbation/analysis/humidity_seasonal_three_islands_train2013_2017/peak_shift_three_islands_2018_2019.png`, `code/perturbation/analysis/humidity_seasonal_three_islands_train2013_2017/peak_shift_three_islands_2018_2019_zoom.png`.
- Producer: `code/analysis/2013-2017 training/to define weak shift using humidity 2013-2017.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Frozen dependency set: `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_fit_cache.rds` (A), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_result.rds` (A).
- Availability: All current humidity fit/result inputs included in A; static dependency completeness, runtime not verified.
- Computation: Saved-draw deterministic trajectories A=actual humidity vs B=training same-week humidity; same states, continuous 104 weeks. Original producer computes numerical results before plotting. Zoom is a distinct view, not resolution duplicate.
- Seed/settings: no new random observation noise; saved fitted draws determine result.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/to define weak shift using humidity 2013-2017.R"
```

### humidity_three_island_plot_2016_2019

- Supplied outputs: `code/analysis/2013-2015 training/plots_hokkaido_shikoku_ryukyu/hokkaido_shikoku_ryukyu_prediction_2016_2019.png`.
- Producer: `code/analysis/2013-2015 training/plot_hokkaido_shikoku_ryukyu_2016_2019.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Frozen dependency set: `code/perturbation/analysis/humidity_seasonal_hokkaido_train2013_2015/humidity_seasonal_fit_cache.rds` (B), `code/perturbation/analysis/humidity_seasonal_hokkaido_train2013_2015/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_shikoku_train2013_2015/humidity_seasonal_fit_cache.rds` (B), `code/perturbation/analysis/humidity_seasonal_shikoku_train2013_2015/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_ryukyu_train2013_2015/humidity_seasonal_fit_cache.rds` (B), `code/perturbation/analysis/humidity_seasonal_ryukyu_train2013_2015/humidity_seasonal_result.rds` (B).
- Availability: See explicit dependency tiers; no runtime verification.
- Computation: Reads saved predictions and checks fit/input identity; no new trajectory generation. Requires archived RDS. Distinct training-window alternative.
- Seed/settings: none for plotting.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2015 training/plot_hokkaido_shikoku_ryukyu_2016_2019.R"
```

### humidity_five_island_plot_2016_2019

- Supplied outputs: `code/perturbation/analysis/humidity_seasonal_all_train2013_2015/humidity_seasonal_all_islands_2016_2019.png`.
- Producer: `code/analysis/2013-2015 training/humidity_seasonal_all plots.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data/japan_rsv.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Frozen dependency set: `code/perturbation/analysis/humidity_seasonal_hokkaido_train2013_2015/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_honshu_train2013_2015/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_shikoku_train2013_2015/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_kyushu_train2013_2015/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_ryukyu_train2013_2015/humidity_seasonal_result.rds` (B).
- Availability: See explicit dependency tiers; no runtime verification.
- Computation: Saved-curve plotting, five-island historical alternative. Peak-shift CSVs/figures in same output folder have no producer established in A and must be separately marked unresolved.
- Seed/settings: See source and configs/RECORDED_SETTINGS.md; unknown historical values not invented.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2015 training/humidity_seasonal_all plots.R"
```

### humidity_five_island_plot_2018_2019

- Supplied outputs: `code/perturbation/analysis/humidity_seasonal_all_islands.png`.
- Producer: `code/analysis/2013-2017 training/humidity_seasonal_all plots.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data/japan_rsv.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Frozen dependency set: `code/perturbation/analysis/humidity_seasonal_hokkaido/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_standalone/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_shikoku/humidity_seasonal_result.rds` (A), `code/perturbation/analysis/humidity_seasonal_kyushu/humidity_seasonal_result.rds` (B), `code/perturbation/analysis/humidity_seasonal_ryukyu/humidity_seasonal_result.rds` (A).
- Availability: See explicit dependency tiers; no runtime verification.
- Computation: Historical five-island saved curves; standalone is Honshu. Current source marks Honshu/Kyushu scripts unused, so acceptance unresolved.
- Seed/settings: See source and configs/RECORDED_SETTINGS.md; unknown historical values not invented.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/2013-2017 training/humidity_seasonal_all plots.R"
```

### humidity_annual_description

- Supplied outputs: `figure1_reproduced.png (not found within supplied A root; output directory lies outside supplied folders)`.
- Producer: `code/analysis/plot_specific_humidity_before_after_2016.R`.
- Inputs: `code/perturbation/data/japan_pop.csv`, `code/perturbation/data_processed/data_processed_humidity.csv`.
- Availability: See explicit dependency tiers; no runtime verification.
- Computation: No fit needed, computes annual/period means from retained humidity CSV. Output path must be relocated into package; no output generated in cleanup.
- Seed/settings: See source and configs/RECORDED_SETTINGS.md; unknown historical values not invented.

```powershell
Rscript --vanilla code/run_handoff.R "code/analysis/plot_specific_humidity_before_after_2016.R"
```

## Original manuscript results

See `docs/PROVENANCE_original.md` for every original retained PDF, missing main figure, exact producer, input data, original Stan variant, fitted/trajectory dependency and owner-run command. Original main Figures 2–4 and several standalone fits/trajectories are absent; no complete manuscript reconstruction is claimed.

## Historical and unresolved result families

| Family | Evidence retained | Limitation / policy |
|---|---|---|
| `analysis/humidity_seasonal_*_train2013_2016/` | Predictions/metrics/peaks/figures and B fitted objects | 2013–2016 producing scripts not found; cannot substitute 2013–2015 or 2013–2017 source. |
| `analysis/humidity_tuning_run1/` | Candidate settings, comparisons, predictions, failures, compact result objects; larger fit/manifest objects B | No standalone tuning orchestrator found. Retain historical outputs; no accepted winner inferred. |
| `analysis/five_island_tuning_run1/` | Job plan, score weights, candidate catalogue, ranked/best model tables, compact result/manifest objects; fit attempts B | Holdout and reconstruction differ; exact producer absent. No automatic best-model adoption. |
| `analysis/baseline_humidity_seasonal_three_islands/` | Small Hokkaido/Shikoku results A; larger results and fits B | Earlier baseline alternative, not current npieff comparison. Producer unresolved. |
| `analysis/humidity_seasonal_standalone/`, `humidity_seasonal_kyushu/` | Historical CSV/figures A and RDS B | Current corresponding source files marked `(사용x)`; historical output/version match not established. |
| `analysis/peak_shift_*`, `humidity_seasonal_all_train2013_2015/peak_shift_*`, `beta_w_training_comparison/` | Existing tables/figures | Precise standalone producer not identified. Session history contains candidate code but does not prove provenance. |
| `stanfit_sirs/honshu_*train2013_2017.rda`, humidity checkpoint | Exact original fits/checkpoint B plus histories | References to absent Desktop fitting scripts; no invented entrypoint. |
| Historical npieff humidity Shikoku pair | Cache/result B, attempt metadata A | Earlier paired-npieff variant, distinct from current no-npieff humidity. |
| Manuscript PDFs/TeX/Bib | Supplied documents A | Intended source mapping only; missing figure files prevent claiming a complete build. |
| R workspaces, histories, before-refit/write snapshots, differing editor source | Exact B hashes and paths | Unknown contents preserved; missing standalone objects may or may not be present. |

All unresolved individual outputs are explicitly marked in `RESULT_PROVENANCE.csv`. They are retained for review, not called accepted final results. No arbitrary parameters, numerical approximations, new scores or fictitious provenance were introduced.

## Recorded ambiguities

- npieff README describes superseded paired-npieff model and nonexistent/misplaced entrypoints.
- fit_baseline_three_islands_2013_2017.R currently runs Hokkaido only; historical other-island settings are in opaque fit caches.
- Baseline validation_only FALSE ignores option set by orchestrator; do not advertise option as safe non-sampling command.
- 2013-2015 Honshu script force_refit TRUE; output filenames say 2018_2019 despite code's 2016-2019 holdout; require owner inspection of dates.
- 2013-2015 to define weak shift using humidity.R depends on preexisting globals; not standalone and does not save listed combined-folder outputs.
- No explicit license found in A; Park et al. attribution alone does not establish redistribution permission.
- Historical ready/diagnostics CSVs are retained evidence, not runtime verification performed in cleanup.
- All .Rhistory files retained in-place in Tier B for scientific provenance; none deleted.
