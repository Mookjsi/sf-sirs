# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("figure")

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="sans"))
library(egg)
library(rstan)
source("../R/firstup.R")
source("../script/script_data.R")

island <- c("hokkaido", "honshu", "shikoku", "kyushu", "ryukyu")

fitlist <- betalist <- npilist <- Slist <- vector('list', length(island))

for (i in 1:length(island)) {
  print(i)
  ii <- island[[i]]
  
  file <- paste0("../stanfit_sirs/stanfit_sirs_", ii, "_npi.rda")
  
  load(file)
  
  ff <- get(paste0("stanfit_sirs_", ii, "_npi"))
  
  ss <- summary(ff)
  
  fitlist[[i]] <- data.frame(
    week=rep(1:52, 11),
    year=rep(2013:2023, each=52),
    est=ss$summary[grepl("C\\[", rownames(ss$summary)),6],
    lwr=ss$summary[grepl("C\\[", rownames(ss$summary)),4],
    upr=ss$summary[grepl("C\\[", rownames(ss$summary)),8],
    island=firstup(ii)
  )
  
  npilist[[i]] <- data.frame(
    week=rep(1:52, 11),
    year=rep(2013:2023, each=52),
    est=ss$summary[grepl("npieff\\[", rownames(ss$summary)),6],
    lwr=ss$summary[grepl("npieff\\[", rownames(ss$summary)),4],
    upr=ss$summary[grepl("npieff\\[", rownames(ss$summary)),8],
    island=firstup(ii)
  )
  
  Slist[[i]] <- data.frame(
    week=rep(1:52, 11),
    year=rep(2013:2023, each=52),
    est=ss$summary[grepl("S\\[", rownames(ss$summary)),6],
    lwr=ss$summary[grepl("S\\[", rownames(ss$summary)),4],
    upr=ss$summary[grepl("S\\[", rownames(ss$summary)),8],
    island=firstup(ii)
  )
  
  betalist[[i]] <- data.frame(
    week=1:52,
    est=ss$summary[grepl("beta\\[", rownames(ss$summary)),6],
    lwr=ss$summary[grepl("beta\\[", rownames(ss$summary)),4],
    upr=ss$summary[grepl("beta\\[", rownames(ss$summary)),8],
    island=firstup(ii)
  )  
}

betadata <- betalist %>%
  bind_rows %>%
  mutate(
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

japan_humidity_island2 <- japan_humidity_island %>%
  # filter(year < 2020) %>%
  group_by(island, week) %>%
  summarise(humidity=mean(humidity))

betadata2 <- betadata %>%
  merge(japan_humidity_island2)

g1 <- ggplot(betadata2) +
  geom_point(aes(humidity*1000, est, col=island, shape=island)) +
  geom_smooth(aes(humidity*1000, est), col="black",
              method="gam") +
  scale_x_continuous("Mean specific humidity (g/kg)") +
  scale_y_continuous("Transmission rate\n(1/week)", limits=c(5, 10.5)) +
  scale_color_viridis_d(end=0.95) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1),
    legend.title = element_blank()
  )

ggsave("figure_joint_climate.pdf", g1, width=6, height=4)
