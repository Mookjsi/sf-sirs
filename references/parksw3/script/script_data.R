library(tidyr)
library(dplyr)

japan_pop_prefecture <- read.csv("../data/japan_pop.csv") %>%
  mutate(
    prefecture=gsub("-.*", "", prefecture),
    pop=as.numeric(gsub(",", "", pop))
  ) %>%
  arrange(prefecture)

japan_pop_island <- japan_pop_prefecture %>%
  group_by(island) %>%
  summarize(
    pop=sum(pop)
  )

japan_rsv_prefecture <- read.csv("../data/japan_rsv.csv") %>%
  mutate(
    week=ifelse(week==53, 52, week)
  ) %>%
  group_by(prefecture, year, week) %>%
  summarize(
    cases=round(mean(cases)),
    long=long[1],
    lat=lat[1]
  ) %>%
  merge(japan_pop_prefecture)

japan_rsv_island <- japan_rsv_prefecture %>%
  group_by(island, year, week) %>%
  summarize(
    cases=sum(cases)
  ) %>%
  group_by(island) %>%
  mutate(
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

japan_humidity_prefecture <- read.csv("../data_processed/data_processed_humidity.csv") %>%
  rename(prefecture=key) %>%
  merge(japan_pop_prefecture)

japan_humidity_island <- japan_humidity_prefecture %>%
  group_by(island, year, week) %>%
  summarize(
    humidity=mean(humidity)
  )

japan_temperature_prefecture <- read.csv("../data_processed/data_processed_temperature.csv") %>%
  rename(prefecture=key) %>%
  merge(japan_pop_prefecture)

japan_temperature_island <- japan_temperature_prefecture %>%
  group_by(island, year, week) %>%
  summarize(
    temperature=mean(temperature)
  )
