# Static artifact and saved-result checks; --replay compares owner-produced CSVs only.
args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
module <- normalizePath(file.path(dirname(script),'../../extensions/humidity'),mustWork=TRUE)
lines <- readLines(file.path(module,'SHA256SUMS.txt'),warn=FALSE)
for (line in lines) {
 match <- regexec('^([0-9a-f]{64})[[:space:]]+[*]?(.+)$',line)
 fields <- regmatches(line,match)[[1]]
 if (length(fields)!=3L) next
 path <- file.path(module,fields[3])
 if (!file.exists(path)) { cat('MISSING_INPUT\t',fields[3],'\n',sep=''); next }
 actual <- tolower(unname(system2('sha256sum',path,stdout=TRUE)))
 actual <- substr(actual,1,64)
 if (actual!=fields[2]) stop('DEPENDENCY_MISMATCH: ',fields[3])
}
base <- file.path(module,'code/perturbation/analysis/humidity_seasonal_three_islands_train2013_2017')
summary <- read.csv(file.path(base,'peak_shift_summary_three_islands_2018_2019.csv'))
if (!setequal(unique(summary$island),c('Hokkaido','Shikoku','Ryukyu')) || !setequal(unique(summary$year),2018:2019)) stop('DEPENDENCY_MISMATCH: island/year scope')
if (anyDuplicated(summary[c('island','year')]) || any(!is.finite(summary$median_shift_weeks))) stop('DEPENDENCY_MISMATCH: summary keys or values')
distribution <- read.csv(file.path(base,'peak_shift_distribution_three_islands_2018_2019.csv'))
if (!setequal(unique(distribution$island),unique(summary$island)) || !setequal(unique(distribution$year),unique(summary$year))) stop('DEPENDENCY_MISMATCH: distribution scope')
if (anyDuplicated(distribution[c('island','year','shift_weeks')]) || any(!is.finite(distribution$proportion))) stop('DEPENDENCY_MISMATCH: distribution keys or values')
for (i in seq_len(nrow(summary))) {
 subset <- distribution[distribution$island==summary$island[i] & distribution$year==summary$year[i],]
 if (abs(sum(subset$proportion)-1)>1e-6 || sum(subset$n)!=summary$n_draws[i]) stop('DEPENDENCY_MISMATCH: distribution totals')
}
if (length(args)==2L && args[1]=='--replay') {
 replay <- read.csv(args[2]); reference <- summary
 if (!identical(names(replay),names(reference)) || !identical(replay[c('island','year')],reference[c('island','year')])) stop('DEPENDENCY_MISMATCH: replay schema or keys')
 numeric <- names(reference)[vapply(reference,is.numeric,logical(1))]
 for (key in numeric) if (any(abs(replay[[key]]-reference[[key]])>1e-6,na.rm=TRUE)) stop('DEPENDENCY_MISMATCH: replay ',key)
}
cat('READY\tavailable checksum records and summary scope\n')
