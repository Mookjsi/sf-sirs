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

island <- "honshu"

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

fitdata <- fitlist %>%
  bind_rows %>%
  merge(japan_rsv_island) %>%
  mutate(
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

npidata <- npilist %>%
  bind_rows %>%
  mutate(
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

Sdata <- Slist %>%
  bind_rows %>%
  merge(japan_pop_island) %>%
  mutate(
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

betadata <- betalist %>%
  bind_rows %>%
  mutate(
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

g1 <- ggplot(fitdata) +
  geom_vline(xintercept = 2013:2023, lty=3, col="gray50", lwd=0.5) +
  # geom_ribbon(data=fitdata %>% filter(type=="Without imported cases"),
  #             aes(year+(week-1)/52, ymin=lwr, ymax=upr), fill="#40B0A6", alpha=0.5) +
  geom_ribbon(aes(year+(week-1)/52, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.5) +
  geom_point(aes(year+(week-1)/52, cases), shape=1) +
  # geom_line(data=fitdata %>% filter(type=="Without imported cases"),
  #           aes(year+(week-1)/52, est), col="#40B0A6", lwd=0.7)  +
  geom_line(aes(year+(week-1)/52, est), col="#E02938", lwd=0.7)  +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Cases") +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g2 <- ggplot(betadata) +
  # geom_ribbon(data=fitdata %>% filter(type=="Without imported cases"),
  #             aes(year+(week-1)/52, ymin=lwr, ymax=upr), fill="#40B0A6", alpha=0.5) +
  geom_ribbon(aes(week, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.5) +
  # geom_line(data=fitdata %>% filter(type=="Without imported cases"),
  #           aes(year+(week-1)/52, est), col="#40B0A6", lwd=0.7)  +
  geom_line(aes(week, est), col="#E02938", lwd=0.7)  +
  scale_x_continuous("Week", expand=c(0, 0)) +
  scale_y_continuous("Transmission rate\n(1/week)", limits=c(4.4, 12)) +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

japan_humidity_island2 <- japan_humidity_island %>%
  # filter(year < 2020) %>%
  group_by(island, week) %>%
  summarise(humidity=mean(humidity))

japan_temperature_island2 <- japan_temperature_island %>%
  # filter(year < 2020) %>%
  group_by(island, week) %>%
  summarise(temperature=mean(temperature))

betadata2 <- betadata %>%
  merge(japan_humidity_island2)

betadata_t2 <- betadata %>%
  merge(japan_temperature_island2)

betadata_comb <- betadata %>%
  merge(japan_temperature_island2) %>%
  merge(japan_humidity_island2)

g3 <- ggplot(betadata2) +
  geom_point(aes(humidity*1000, est)) +
  geom_smooth(aes(humidity*1000, est), col="blue", fill="blue") +
  geom_smooth(aes(humidity*1000, est), col="#E02938", fill="#E02938",
              method="lm", formula=y~x+I(x^2)) +
  scale_x_continuous("Mean specific humidity (g/kg)") +
  scale_y_continuous("Transmission rate\n(1/week)", limits=c(5, 10.5)) +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g3a <- ggplot(betadata_t2) +
  geom_point(aes(temperature, est)) +
  geom_smooth(aes(temperature, est), col="blue", fill="blue") +
  geom_smooth(aes(temperature, est), col="#E02938", fill="#E02938",
              method="lm", formula=y~x+I(x^2)) +
  scale_x_continuous("Mean temperature (°C)") +
  scale_y_continuous("Transmission rate\n(1/week)", limits=c(5, 10.5)) +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g3b <- ggplot(betadata_comb) +
  geom_point(aes(temperature, humidity)) +
  scale_x_continuous("Mean temperature (°C)") +
  scale_y_continuous("Mean specific humidity (g/kg)") +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

lapply(split(betadata_comb, betadata_comb$island), function(x) {
  cor(x$temperature, x$humidity)
}) 

lapply(split(betadata_comb, betadata_comb$island), function(x) {
  summary(lm(est~1+I(humidity*1000)+I((humidity*1000)^2)+temperature+I(temperature^2), data=x))
}) 

lapply(split(betadata_comb, betadata_comb$island), function(x) {
  summary(lm(est~1+I(humidity*1000)+I((humidity*1000)^2), data=x))
}) 

lapply(split(betadata_comb, betadata_comb$island), function(x) {
  summary(lm(est~1+temperature+I(temperature^2), data=x))
}) 

g4 <- ggplot(npidata) +
  geom_vline(xintercept = 2013:2023, lty=3, col="gray50", lwd=0.5) +
  geom_hline(yintercept = 1, lty=2) +
  # geom_ribbon(data=fitdata %>% filter(type=="Without imported cases"),
  #             aes(year+(week-1)/52, ymin=lwr, ymax=upr), fill="#40B0A6", alpha=0.5) +
  geom_ribbon(aes(year+(week-1)/52, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.5) +
  # geom_line(data=fitdata %>% filter(type=="Without imported cases"),
  #           aes(year+(week-1)/52, est), col="#40B0A6", lwd=0.7)  +
  geom_line(aes(year+(week-1)/52, est), col="#E02938", lwd=0.7)  +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Relative changes\nin transmission", expand=c(0, 0)) +
  coord_cartesian(ylim=c(0, 1.6)) +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g5 <- ggplot(Sdata) +
  geom_vline(xintercept = 2013:2023, lty=3, col="gray50", lwd=0.5) +
  # geom_ribbon(data=fitdata %>% filter(type=="Without imported cases"),
  #             aes(year+(week-1)/52, ymin=lwr, ymax=upr), fill="#40B0A6", alpha=0.5) +
  geom_ribbon(aes(year+(week-1)/52, ymin=lwr/pop, ymax=upr/pop), fill="#E02938", alpha=0.5) +
  # geom_line(data=fitdata %>% filter(type=="Without imported cases"),
  #           aes(year+(week-1)/52, est), col="#40B0A6", lwd=0.7)  +
  geom_line(aes(year+(week-1)/52, est/pop), col="#E02938", lwd=0.7)  +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Proportion\nsusceptible", limits=c(0, 0.26), expand=c(0, 0)) +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

gcomb <- ggarrange(g1, g2, g3, g4, g5, nrow=5,
                   labels=LETTERS[1:5])

ggsave("figure_comb_sirs_npi.pdf", gcomb, width=12, height=10)
ggsave("figure_comb_sirs_npi_temp.pdf", g3a, width=12, height=3)
ggsave("figure_comb_sirs_npi_temp_hum.pdf", g3b, width=12, height=3)
