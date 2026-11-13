# Static availability report; no packages are installed or research code run.
args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
root <- normalizePath(file.path(dirname(script), '../..'), mustWork=TRUE)
module <- file.path(root, 'extensions/humidity')
required <- c('code/run_handoff.R', 'ARCHIVE_MANIFEST.csv', 'ARTIFACT_MANIFEST.csv', 'RESULT_PROVENANCE.csv',
 'code/perturbation/data/japan_pop.csv', 'code/perturbation/data_processed/data_processed_humidity.csv',
 'code/analysis/2013-2017 training/plot_forecast_metrics.R',
 'code/perturbation/analysis/humidity_seasonal_three_islands_train2013_2017/peak_shift_summary_three_islands_2018_2019.csv')
for (island in c('hokkaido','shikoku','ryukyu')) for (name in c('humidity_seasonal_fit_cache.rds','humidity_seasonal_result.rds'))
 required <- c(required, file.path('code/perturbation/analysis',paste0('humidity_seasonal_',island),name))
for (p in required) cat(if (file.exists(file.path(module,p))) 'READY' else 'MISSING_INPUT', p, '\n', sep='\t')
archive <- file.path(module,'code/perturbation/analysis/npieff_baseline_humidity_train2013_2017_forecast2018_2019')
cat('ARCHIVE_REQUIRED\tbaseline + npieff fitted caches; consult ARCHIVE_MANIFEST.csv\n')
req <- readLines(file.path(module,'requirements-R.txt'))
for (line in req[!grepl('^#|^$', req)]) {
 parts <- strsplit(line,'==',fixed=TRUE)[[1]]
 have <- if (requireNamespace(parts[1],quietly=TRUE)) as.character(packageVersion(parts[1])) else NA_character_
 cat(if (is.na(have)) 'MISSING_INPUT' else if (have==parts[2]) 'READY' else 'DEPENDENCY_MISMATCH',
     paste0('package:',line,'; installed=',have), '\n', sep='\t')
}
