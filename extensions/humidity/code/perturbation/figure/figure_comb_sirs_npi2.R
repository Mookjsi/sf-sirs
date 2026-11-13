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

island <- c("honshu", "shikoku", "kyushu")

fitlist <- betalist <- beta2list <- npilist <- Slist <- vector('list', length(island))

for (i in 1:length(island)) {
  print(i)
  ii <- island[[i]]
  
  file <- paste0("../stanfit_sirs2/stanfit_sirs_", ii, "_npi2.rda")
  
  load(file)
  
  ff <- get(paste0("stanfit_sirs_", ii, "_npi2"))
  
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
  
  beta2list[[i]] <- data.frame(
    week=1:52,
    est=ss$summary[grepl("beta2\\[", rownames(ss$summary)),6],
    lwr=ss$summary[grepl("beta2\\[", rownames(ss$summary)),4],
    upr=ss$summary[grepl("beta2\\[", rownames(ss$summary)),8],
    island=firstup(ii)
  )  
}

fitdata <- fitlist %>%
  bind_rows %>%
  merge(japan_rsv_island) %>%
  mutate(
    island=factor(island, levels=c("Honshu", "Shikoku", "Kyushu"))
  )

npidata <- npilist %>%
  bind_rows %>%
  mutate(
    island=factor(island, levels=c("Honshu", "Shikoku", "Kyushu"))
  )

Sdata <- Slist %>%
  bind_rows %>%
  merge(japan_pop_island) %>%
  mutate(
    island=factor(island, levels=c("Honshu", "Shikoku", "Kyushu"))
  )

betadata <- betalist %>%
  bind_rows %>%
  mutate(
    island=factor(island, levels=c("Honshu", "Shikoku", "Kyushu"))
  )

beta2data <- beta2list %>%
  bind_rows %>%
  mutate(
    island=factor(island, levels=c("Honshu", "Shikoku", "Kyushu"))
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
  geom_ribbon(data=beta2data, aes(week, ymin=lwr, ymax=upr), fill="blue", alpha=0.5) +
  geom_line(data=beta2data, aes(week, est), col="blue", lwd=0.7)  +
  scale_x_continuous("Week", expand=c(0, 0)) +
  scale_y_continuous("Transmission rate\n(1/week)") +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g3 <- ggplot(Sdata) +
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

Sdata_summ <- Sdata %>%
  group_by(island) %>%
  filter(week==26)

g4 <- ggplot(Sdata_summ) +
  geom_vline(xintercept = 2013:2023, lty=3, col="gray50", lwd=0.5) +
  # geom_ribbon(data=fitdata %>% filter(type=="Without imported cases"),
  #             aes(year+(week-1)/52, ymin=lwr, ymax=upr), fill="#40B0A6", alpha=0.5) +
  geom_errorbar(aes(year+(week-1)/52, ymin=lwr/pop, ymax=upr/pop), col="#E02938", width=0, lwd=0.7) +
  # geom_line(data=fitdata %>% filter(type=="Without imported cases"),
  #           aes(year+(week-1)/52, est), col="#40B0A6", lwd=0.7)  +
  geom_point(aes(year+(week-1)/52, est/pop), col="#E02938", lwd=0.7)  +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2020)) +
  scale_y_continuous("Proportion\nsusceptible", limits=c(0.08, 0.11), expand=c(0, 0)) +
  facet_wrap(~island, scale="free_y", nrow=1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

gcomb <- ggarrange(g1, g2, g3, g4, ncol=1,
                   labels=LETTERS[1:4])

ggsave("figure_comb_sirs_npi2.pdf", gcomb, width=12, height=10)
