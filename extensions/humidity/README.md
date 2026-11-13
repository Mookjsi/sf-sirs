# Humidity extension

This isolated module preserves the collaborator's Japan RSV specific-humidity analysis. It does not share fitted parameters or equations with the main SF-SIRS M2 pipeline. The retained 2013–2017 training and 2018–2019 conditional hindcast covers Hokkaido, Shikoku, and Ryukyu. Baseline + npieff versus seasonal humidity without npieff changes two model components, so it is not an isolated humidity effect.

The original `code/perturbation/` and `code/analysis/` hierarchies are retained because source scripts rely on them. Earlier training windows and exploratory variants remain in those directories for provenance and are not public entry points. `provenance/` and the original manifests retain path and checksum records. The six frozen humidity RDS files and collaborator data are included for repository tracking with the owner’s redistribution authorization. Baseline + npieff fit caches are absent and require the collaborator archive.

From the repository root, the owner may later run:

```sh
Rscript scripts/humidity/check_dependencies.R
Rscript scripts/humidity/verify.R
Rscript scripts/humidity/run.R "code/analysis/2013-2017 training/plot_forecast_metrics.R"
```

The plot command reads saved metrics and writes to `extensions/humidity/outputs/`. The peak-shift replay script is also supported and writes into the ignored `outputs/` directory. Other original scripts remain blocked by the launcher because some can refit or overwrite shipped outputs. `verify.R --replay PATH_TO_SUMMARY.csv` compares a future rerun with the saved three-island summary using a `1e-6` numeric tolerance. Package versions in `requirements-R.txt` are recorded observations, not a tested lockfile. No scientific runtime validation has been performed during this integration.

The repository owner has confirmed redistribution rights for the included collaborator code, supplied data, frozen RDS objects, and figures. The Park code under `code/perturbation/` carries its original MIT license.
