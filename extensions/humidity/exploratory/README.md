# Exploratory RSV and climate analysis

This Python module reproduces the collaborator's prefecture-level timing-variance and climate figures independently of the SF-SIRS M2 model. It uses four unchanged source CSVs in `data/climate/source/` and writes every generated table, figure, and run manifest to a fresh `artifacts/` directory.

```sh
python3 scripts/reproduce_exploratory.py --input-dir data/climate/source --out artifacts/exploratory-climate-new
python3 scripts/reproduce_exploratory.py --input-dir data/climate/source --out artifacts/exploratory-variance-new --analysis variance
```

Available selections are `variance`, `tokyo`, `temperature_smoothing`, `temperature_heatmap`, `humidity_heatmap`, `humidity_profiles`, and `threshold_overlay`. The latter generates smoothing data first. Install the recorded Python versions in `extensions/humidity/requirements.txt` in a separate environment. The source calculations reside in `exploratory/_source/`; the runner supplies explicit paths and a single output directory.

The Tokyo plot intentionally selects Tokyo. The collaborator's `Tokyo_plot.py` selected Tottori despite its Tokyo title, and archived Figure 8 shows a peak below 250 cases consistent with that error. The corrected plot is a changed figure, not a reproduction of that archived PDF.

Timing variance is the case-weighted population variance of reported-case week positions within a 41-week peak-centered window, not transmission-rate variance. Climate inputs are prefecture observations; the five profiles use Hokkaido, Tokyo, Kochi, Fukuoka, and Okinawa as representatives. Humidity is recorded in kg/kg and converted to g/kg for plots. The 12/14/16°C and 20/22/24°C crossings use smoothed outdoor temperature as exploratory proxies, not observed indoor operation. Climate observations end at 2020 W09 while RSV cases continue through 2023.

The original conceptual `plot5_right.py` is excluded because it is not an empirical SIRS result. The original archive remains outside Git in ignored staging. No SF-SIRS fitting or Park MCMC is performed by this module.
