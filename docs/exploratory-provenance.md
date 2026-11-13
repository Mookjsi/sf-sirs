# Exploratory climate figure provenance

The collaborator archive `final_현승.zip` is retained in ignored staging. Four owner-validated source CSV copies are committed under `data/climate/source/`. The archived source code and figures remain available in staging; the integrated calculations live under `extensions/humidity/exploratory/`. Run-specific numerical tables and figures are generated under ignored `artifacts/` paths.

| Archived figure | Original source | Integrated analysis and output | Status |
|---|---|---|---|
| Figure 1 | `paper_rsv_heatmap.py`; Park Figure 1A reference | `threshold_overlay`; `rsv_heatmap_data.csv` | Park-style plot data reproduced; week 53 is combined for plotting only |
| Figure 2 | `var.py` | `variance`; `all_prefecture_variances.csv`, `variance_heatmap.jpg` | Calculation runs from source cases |
| Figure 3 | `temperature_heatmap.py` | `temperature_heatmap`; `temperature_heatmap.jpg` | Prefecture-level outdoor temperature |
| Figure 4 | `humidity_heatmap.py` | `humidity_heatmap`; `humidity_heatmap.jpg` | Specific humidity converted to g/kg for display |
| Figure 5 | `five_region_humidity.py` | `humidity_profiles`; `humidity_profiles.jpg` | Five representative prefectures, not island aggregates |
| Figure 6 | `temperature_smoothing.py` | `temperature_smoothing`; `temperature_smooth.csv`, `temperature_smoothing.jpg` | Retrospective cubic GCV smoothing, no extrapolation |
| Figure 7 | `paper_rsv_heatmap.py` | `threshold_overlay`; `temperature_threshold_crossings.csv`, `rsv_temperature_thresholds.jpg` | Six exploratory outdoor temperature thresholds |
| Figure 8 | `Tokyo_plot.py` | `tokyo`; `tokyo_observations.csv`, `tokyo_peak_windows.jpg` | Corrected: archived script filtered Tottori under a Tokyo title |

The original Park heatmap implementation and MIT notice already exist in `references/parksw3/`; this module does not copy another version of the R code. The conceptual `plot5_right.py` is excluded from the empirical pipeline. The four input SHA-256 hashes are emitted by each run in `input_hashes.json`.

## Chapter 6 scope

The current repository includes separate R humidity extension code and six frozen 2013–2017 three-island humidity posterior objects. Those support the saved Hokkaido, Shikoku, and Ryukyu A/B analysis, subject to runtime verification. The full Chapter 6 `humidity` and `humidity_seasonal` model comparison across training windows and islands is not established by this Python exploratory module. Some baseline + npieff fitted caches remain archival dependencies; complete reproduction also needs the exact model inputs, fitted objects, and matching environment for each reported training window. No substitute Bayesian fitting was performed here.
