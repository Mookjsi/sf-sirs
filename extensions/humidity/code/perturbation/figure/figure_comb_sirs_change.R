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

load("../simulate_sirs/simulate_sirs_honshu.rda")
load("../simulate_sirs/simulate_sirs_shikoku.rda")
load("../simulate_sirs/simulate_sirs_kyushu.rda")

island <- c("honshu", "shikoku", "kyushu")

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
  merge(japan_pop_island) %>%
  mutate(
    island=factor(island, levels=c("Honshu", "Shikoku", "Kyushu"))
  )

fitdata_filter <- fitdata %>% 
  filter(week==25, year < 2020) %>%
  mutate(
    yy=gsub("20", "", year)
  )

allsim <- bind_rows(
  simulate_sirs_honshu,
  simulate_sirs_kyushu,
  simulate_sirs_shikoku
) %>%
  mutate(
    island=factor(island, levels=c("Honshu", "Shikoku", "Kyushu"))
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

g2 <- ggplot(allsim_summ) +
  geom_raster(aes(I0, S0, fill=intensity)) +
  geom_vline(xintercept=unique(allsim$I0)[9], lty=2, color="white") +
  geom_point(data=fitdata_filter, aes(I/pop, S/pop), size=6, shape=21, fill="white", stroke=0.8, alpha=0.7) +
  geom_text(data=fitdata_filter, aes(I/pop, S/pop, label=yy), size=3) +
  scale_x_log10("I(0)", expand=c(0, 0)) +
  scale_y_continuous("S(0)", expand=c(0, 0)) +
  scale_fill_viridis_c("Epidemic\nintensity") +
  facet_wrap(~island, nrow=1) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
  )

allsim_filter <- allsim %>%
  filter(I0==unique(allsim$I0)[9],
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

load("../stanfit_sirs/stanfit_sirs_honshu_npi.rda")
load("../stanfit_sirs/stanfit_sirs_kyushu_npi.rda")
load("../simulate_sirs/simulate_sirs_interpolate.rda")

ss_honshu <- summary(stanfit_sirs_honshu_npi)
ss_kyushu <- summary(stanfit_sirs_kyushu_npi)

beta_honshu <- ss_honshu$summary[grepl("beta\\[",rownames(ss_honshu$summary)),6]
beta_kyushu <- ss_kyushu$summary[grepl("beta\\[",rownames(ss_kyushu$summary)),6]

week_gam <- 1:52
week_gam_pred <- seq(1, 52, by=0.1)

gfit_honshu <- gam(beta_honshu~s(week_gam, bs="cc"))
gfit_kyushu <- gam(beta_kyushu~s(week_gam, bs="cc"))

beta_smooth_honshu <- predict(gfit_honshu, newdata=data.frame(week_gam=week_gam_pred))
beta_smooth_kyushu <- predict(gfit_kyushu, newdata=data.frame(week_gam=week_gam_pred))

seas_smooth_honshu <- beta_smooth_honshu/mean(beta_smooth_honshu)-1
seas_smooth_kyushu <- beta_smooth_kyushu/mean(beta_smooth_kyushu)-1

amp_honshu <- (max(seas_smooth_honshu)-min(seas_smooth_honshu))/2
amp_kyushu <- (max(seas_smooth_kyushu)-min(seas_smooth_kyushu))/2

g4 <- ggplot() +
  geom_line(aes(week_gam_pred, seas_smooth_honshu)) +
  annotate("text", x=-Inf, y=Inf, label="Honshu", hjust=0, vjust=1.2, family="serif") +
  scale_x_continuous("Week", expand=c(0, 0)) +
  scale_y_continuous("Seasonality", limits=c(-0.3, 0.33), expand=c(0, 0)) +
  theme(
    panel.grid = element_blank()
  )

g5 <- ggplot() +
  geom_line(aes(week_gam_pred, seas_smooth_honshu), col="gray") +
  geom_line(aes(week_gam_pred, seas_smooth_honshu/amp_honshu*amp_kyushu)) +
  annotate("text", x=-Inf, y=Inf, label="Honshu, increased amplitude", hjust=0, vjust=1.2, family="serif") +
  scale_x_continuous("Week", expand=c(0, 0)) +
  scale_y_continuous("Seasonality", limits=c(-0.3, 0.33), expand=c(0, 0)) +
  theme(
    panel.grid = element_blank()
  )

g6 <- ggplot() +
  geom_line(aes(week_gam_pred, seas_smooth_honshu), col="gray") +
  geom_line(aes(week_gam_pred, seas_smooth_kyushu/amp_kyushu*amp_honshu)) +
  annotate("text", x=-Inf, y=Inf, label="Kyushu, decreased amplitude", hjust=0, vjust=1.2, family="serif") +
  scale_x_continuous("Week", expand=c(0, 0)) +
  scale_y_continuous("Seasonality", limits=c(-0.3, 0.33), expand=c(0, 0)) +
  theme(
    panel.grid = element_blank()
  )

g7 <- ggplot() +
  geom_line(aes(week_gam_pred, seas_smooth_honshu), col="gray") +
  geom_line(aes(week_gam_pred, seas_smooth_kyushu)) +
  annotate("text", x=-Inf, y=Inf, label="Kyushu", hjust=0, vjust=1.2, family="serif") +
  scale_x_continuous("Week", expand=c(0, 0)) +
  scale_y_continuous("Seasonality", limits=c(-0.3, 0.33), expand=c(0, 0)) +
  theme(
    panel.grid = element_blank()
  )

g8 <- ggplot(simulate_sirs_interpolate) +
  geom_raster(aes(amp, interp, fill=shift_cases)) +
  annotate("text", x=0.15, y=0, label="C", hjust=0, vjust=1,
           size=5, fontface = "bold") +
  annotate("text", x=0.247, y=0, label="D", hjust=1, vjust=1,
           size=5, fontface = "bold") +
  annotate("text", x=0.15, y=1, label="E", hjust=0, vjust=0,
           size=5, fontface = "bold") +
  annotate("text", x=0.247, y=1, label="F", hjust=1, vjust=0,
           size=5, fontface = "bold") +
  scale_x_continuous("Seasonal transmission amplitude", expand=c(0, 0), position="top") +
  scale_y_reverse("Interpolation coefficient", expand=c(0, 0),
                  breaks=c(0, 0.25, 0.5, 0.75, 1),
                  labels=c("Honshu\n0.00", 0.25, "0.50", 0.75, "Kyushu\n1.00")) +
  scale_fill_steps("Differences in peak timing (weeks)",
                   limits=c(12, 18),
                   low="white",
                   high="#FF495C") +
  # scale_fill_viridis_c() +
  theme(
    legend.position = "bottom"
  )

g8a <- ggarrange(g8, labels=c("G"))

gtmp1 <- ggarrange(g4, g6, labels=c("C", "E"))
gtmp2 <- ggarrange(g5, g7, labels=c("D", "F"))

gcomb2 <- arrangeGrob(gtmp1, g8a, gtmp2,
                     widths=c(1, 1.5, 1))

gfinal <- arrangeGrob(gcomb1, gcomb2, nrow=2, heights=c(1.2, 1))

ggsave("figure_comb_sirs_change.pdf", gfinal, width=10, height=8)
ggsave("figure_comb_sirs_change.png", gfinal, width=10, height=8)
