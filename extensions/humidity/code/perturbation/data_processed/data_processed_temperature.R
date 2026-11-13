# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("data_processed")

library(tidyr)
library(dplyr)
library(lubridate)

data_climate <- read.csv("../data_climate/era5.t2m.daily.Japan.subregion.201301-202002.csv")

data_processed_temperature <- data_climate %>%
  gather(key, value, -time) %>%
  mutate(
    year=epiyear(time),
    week=epiweek(time),
    week=ifelse(week==53, 52, week)
  ) %>%
  group_by(key, year, week) %>%
  summarize(
    temperature=mean(value)
  )

write.csv(data_processed_temperature, "data_processed_temperature.csv")
