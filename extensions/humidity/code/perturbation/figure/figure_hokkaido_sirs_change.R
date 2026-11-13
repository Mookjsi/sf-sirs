# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("figure")

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="sans"))
library(egg)
library(rstan)
library(mgcv)
source("../R/firstup.R")
source("../script/script_data.R")

load("../simulate_sirs/simulate_sirs_hokkaido.rda")

island <- c("hokkaido")

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
    S=ss$summary[grepl("S\\[", rownames(ss$summary)),6],
    I=ss$summary[grepl("I\\[", rownames(ss$summary)),6],
    island=firstup(ii)
  )
}

fitdata <- fitlist %>%
  bind_rows %>%
  merge(japan_pop_island)

fitdata_filter <- fitdata %>% 
  filter(week==25, year < 2020) %>%
  mutate(
    yy=gsub("20", "", year)
  )

allsim <- bind_rows(
  simulate_sirs_hokkaido
)

allsim_summ <- allsim %>%
  group_by(island, I0, S0) %>%
  mutate(
    p=cases/sum(cases)
  ) %>%
  summarize(
    cog=sum(cases*1:52)/sum(cases)+26,
    intensity=-1/sum(p * log(p))
  )

g1 <- ggplot(allsim_summ) +
  geom_raster(aes(I0, S0, fill=cog)) +
  geom_vline(xintercept=unique(allsim$I0)[9], lty=2, color="white") +
  geom_point(data=fitdata_filter, aes(I/pop, S/pop), size=6, shape=21, fill="white", stroke=0.8, alpha=0.7) +
  geom_text(data=fitdata_filter, aes(I/pop, S/pop, label=yy), size=3) +
  scale_x_log10("I(0)", expand=c(0, 0)) +
  scale_y_continuous("S(0)", expand=c(0, 0)) +
  scale_fill_viridis_c("Center of gravity\n(weeks)",
                       breaks=c(40, 46, 52, 58, 64),
                       labels=c(40, 46, 52, 6, 12)) +
  facet_wrap(~island, nrow=1) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
  )

allsim_filter <- allsim %>%
  filter(I0==unique(allsim$I0)[12],
         S0 %in% unique(allsim$S0)[1:10*2]) %>%
  mutate(
    week=factor(week,
                levels=c(26:52, 1:25))
  )

allsim_filter %>%
  group_by(island, S0) %>%
  filter(cases==max(cases)) %>%
  as.data.frame

g3 <- ggplot(allsim_filter) +
  geom_line(aes(week, cases, group=S0, col=S0)) +
  scale_x_discrete("Week",
                     breaks=c("26", "32", "38", "44", "50", "4", "10", "16", "22"),
                   expand=c(0, 0)) +
  scale_y_continuous("Cases", limits=c(0, NA), expand=c(0, 0)) +
  scale_colour_gradient("S(0)", low = "gray", high = "black") +
  facet_wrap(~island, nrow=1, scale="free_y") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

gcomb1 <- ggarrange(g1, g3, nrow=2, labels=c("A", "B"))

ggsave("figure_hokkaido_sirs_change.pdf", gcomb1, width=6, height=8)
