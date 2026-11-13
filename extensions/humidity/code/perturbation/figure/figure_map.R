# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("figure")

library(tidyr)
library(dplyr)
library(mapdata)
library(ggplot2)
theme_set(theme_bw(base_family="sans"))
library(ggfortify)
source("../script/script_data.R")

jp <- ggplot2::map_data('world2', 'japan')

jp_filter <- jp %>%
  filter(subregion %in% c("Amami O Shima", "Iriomote Jima", "Ishigaki Shima", "Miyako Jima", "Tokuno Shima", "Okinawa", "Shikoku", "Honshu", "Kyushu",
                          "Hokkaido")) %>%
  mutate(
    subregion=factor(subregion,
                     levels=c("Amami O Shima", "Iriomote Jima", "Ishigaki Shima", "Miyako Jima", "Tokuno Shima", "Okinawa", "Shikoku", "Honshu", "Kyushu",
                              "Hokkaido"),
                     labels=c("Ryukyu", "Ryukyu", "Ryukyu", "Ryukyu", "Ryukyu", "Ryukyu", "Shikoku", "Honshu", "Kyushu",
                              "Hokkaido")),
    subregion=factor(subregion,
                     levels=c("Hokkaido",
                              "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

g1 <- ggplot(jp, aes(x = long, y = lat, group = group)) +
  geom_polygon() +
  geom_polygon(data=jp_filter, aes(fill=subregion)) +
  scale_x_continuous("Longitude") +
  scale_y_continuous("Latitude") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.2, 0.8),
    legend.title = element_blank()
  )

ggsave("figure_map.pdf", width=4, height = 4)
