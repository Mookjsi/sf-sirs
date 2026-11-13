# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("figure")

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="sans"))
library(egg)
library(rstan)
library(readxl)
source("../R/firstup.R")
source("../script/script_data.R")

# https://www.cfa.go.jp/assets/contents/node/basic_page/field_ref_resources/f699fe5b-bf3d-46b1-8028-c5f450718d1a/7803b525/20230901_policies_hoiku_torimatome_r5_02.pdf
childcare <- read_xlsx("../data/japan_childcare.xlsx")

childcare2 <- childcare %>%
  gather(key, value, -year) %>%
  group_by(year) %>%
  summarize(
    total=sum(value)
  )

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
    S=ss$summary[grepl("S\\[", rownames(ss$summary)),6],
    Slwr=ss$summary[grepl("S\\[", rownames(ss$summary)),4],
    Supr=ss$summary[grepl("S\\[", rownames(ss$summary)),8],
    I=ss$summary[grepl("I\\[", rownames(ss$summary)),6],
    island=firstup(ii)
  )
}

fitdata <- fitlist %>%
  bind_rows %>%
  merge(japan_pop_island) %>%
  mutate(
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

fitdata_filter <- fitdata %>% 
  filter(week==26, year < 2020,
         island %in% c("Honshu", "Shikoku", "Kyushu")) %>%
  mutate(
    yy=gsub("20", "", year)
  ) %>%
  merge(childcare2)

g1 <- ggplot(fitdata_filter) +
  geom_ribbon(aes(year, ymin=Slwr/pop, ymax=Supr/pop), alpha=0.2) +
  geom_line(aes(year, S/pop)) +
  geom_point(aes(year, S/pop)) +
  geom_line(aes(year, total/2.6e7), col="#E02938") +
  geom_point(aes(year, total/2.6e7), col="#E02938", shape="triangle") +
  scale_x_continuous("Year") +
  scale_y_continuous("Proportion susceptible",
                     sec.axis = sec_axis(~.*26, "Number of children\nattending childcare facilities\n(million)")) +
  facet_wrap(~island, nrow=1) +
  theme(
    strip.background = element_blank(),
    axis.line.y.right = element_line(color="#E02938"),
    axis.ticks.y.right = element_line(color="#E02938"),
    axis.text.y.right = element_text(color="#E02938"),
    axis.title.y.right = element_text(color="#E02938")
  )

lapply(split(fitdata_filter, fitdata_filter$island)[2:4], function(x) {
  cor.test(x$S, x$total)
})

g2 <- ggplot(fitdata_filter) +
  geom_smooth(aes(total/1e6, S/pop), method="lm", col="#E02938", fill="#E02938", fullrange=TRUE) +
  geom_errorbar(aes(total/1e6, ymin=Slwr/pop, ymax=Supr/pop), width=0) +
  geom_point(aes(total/1e6, S/pop)) +
  scale_x_continuous("Number of children\nattending childcare facilities (million)",
                     limits=c(2.15, 2.75),
                     expand=c(0, 0)) +
  scale_y_continuous("Proportion susceptible") +
  facet_wrap(~island, nrow=1) +
  theme(
    strip.background = element_blank()
  )

gcomb <- ggarrange(g1, g2, nrow=2,
                   labels=c("A", "B"))

ggsave("figure_childcare.pdf", gcomb, width=8, height=6)
