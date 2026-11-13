# Reproduce the authors' script_data.R aggregation using base R only.
# Defining this function does not read or write files.
preprocess_japan_data <- function(input_root, output_dir) {
  input_root <- normalizePath(input_root, mustWork=TRUE)
  paths <- c(
    cases=file.path(input_root,"data/japan_rsv.csv"),
    population=file.path(input_root,"data/japan_pop.csv"),
    humidity=file.path(input_root,"data_processed/data_processed_humidity.csv"),
    temperature=file.path(input_root,"data_processed/data_processed_temperature.csv"))
  if(any(!file.exists(paths))) stop("Missing inputs: ",paste(paths[!file.exists(paths)],collapse=", "))
  if(dir.exists(output_dir) || file.exists(output_dir)) stop("Output already exists: ",output_dir)
  read_input <- function(path) utils::read.csv(path,stringsAsFactors=FALSE)
  pop <- read_input(paths[["population"]])
  pop$prefecture <- gsub("-.*","",pop$prefecture)
  pop$pop <- as.numeric(gsub(",","",pop$pop))
  pop <- pop[order(pop$prefecture),]; rownames(pop) <- NULL
  stopifnot(!anyNA(pop),!anyDuplicated(pop$prefecture),all(pop$pop>0))
  population_island <- aggregate(pop~island,pop,sum)

  raw <- read_input(paths[["cases"]])
  stopifnot(!anyNA(raw),all(raw$cases>=0),all(raw$cases==round(raw$cases)),
            all(raw$week %in% 1:53),all(raw$prefecture %in% pop$prefecture))
  mapping <- data.frame(raw_row=seq_len(nrow(raw)),prefecture=raw$prefecture,
    year=raw$year,raw_week=raw$week,prepared_week=ifelse(raw$week==53,52,raw$week),
    island=pop$island[match(raw$prefecture,pop$prefecture)])
  mapping$original_island_index <- (mapping$year-2013L)*52L+mapping$prepared_week
  if("time" %in% names(raw)) mapping$raw_time <- raw$time
  raw$week[raw$week==53] <- 52
  keys <- c("prefecture","year","week")
  cases <- aggregate(cases~prefecture+year+week,raw,function(x) round(mean(x)))
  # Coordinates follow the first original row for each prefecture/year/week.
  coordinates <- raw[!duplicated(raw[keys]),c(keys,"long","lat")]
  cases <- merge(cases,coordinates,by=keys)
  cases <- cases[order(cases$prefecture,cases$year,cases$week),]
  cases <- merge(cases,pop,by="prefecture",sort=TRUE)
  cases_island <- aggregate(cases~island+year+week,cases,sum)
  cases_island <- cases_island[order(cases_island$island,cases_island$year,cases_island$week),]

  result <- list(japan_pop_prefecture=pop,japan_pop_island=population_island,
                 japan_rsv_prefecture=cases,japan_rsv_island=cases_island)
  unmatched_climate_prefectures <- list()
  for(variable in c("humidity","temperature")) {
    climate <- read_input(paths[[variable]])
    names(climate)[names(climate)=="key"] <- "prefecture"
    # Preserve the original inner join, including unmatched source spellings.
    unmatched_climate_prefectures[[variable]] <- setdiff(climate$prefecture,pop$prefecture)
    climate <- merge(climate,pop,by="prefecture",sort=TRUE)
    stopifnot(!anyNA(climate),all(is.finite(climate[[variable]])))
    island <- aggregate(climate[variable],climate[c("island","year","week")],mean)
    island <- island[order(island$island,island$year,island$week),]
    result[[paste0("japan_",variable,"_prefecture")]] <- climate
    result[[paste0("japan_",variable,"_island")]] <- island
  }
  stopifnot(all(table(cases_island$island)==572),
            identical(sort(unique(cases_island$year)),2013:2023),
            all(vapply(result,function(x) !anyNA(x),logical(1))))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  for(name in names(result)) utils::write.csv(result[[name]],
    file.path(output_dir,paste0(name,".csv")),row.names=FALSE)
  utils::write.csv(mapping,file.path(output_dir,"prepared_row_mapping.csv"),row.names=FALSE)
  provenance <- list(preprocessing_id="original-island-v1",inputs=paths,
    input_md5=tools::md5sum(paths),
    unmatched_climate_prefectures=unmatched_climate_prefectures,
    output_md5=tools::md5sum(file.path(output_dir,paste0(names(result),".csv"))),
    mapping_file="prepared_row_mapping.csv",
    mapping_md5=unname(tools::md5sum(file.path(output_dir,"prepared_row_mapping.csv"))),
    rules=c("week 53 -> 52; prefecture/year/week: round(mean(cases))",
            "island/year/week: sum(cases); island population: sum(pop)",
            "island climate: unweighted prefecture mean"))
  saveRDS(provenance,file.path(output_dir,"provenance.rds"))
  capture.output(dput(provenance),file=file.path(output_dir,"provenance.txt"))
  invisible(result)
}
