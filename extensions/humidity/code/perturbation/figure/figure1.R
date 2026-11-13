# Project environment and working directory (no model changes).
source(file.path(Sys.getenv("PERTURBATION_ROOT",
    unset = Sys.getenv("PERTURBATION_ROOT")), "script", "setup.R"))
perturbation_setup("figure")

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="sans"))
library(egg)
library(gridExtra)
source("../script/script_data.R")

japan_rsv_max <- japan_rsv_prefecture %>%
  filter(year < 2020) %>%
  group_by(prefecture) %>%
  summarize(
    max=max(cases)
  )

japan_rsv_prefecture2 <- japan_rsv_prefecture %>%
  merge(japan_rsv_max) %>%
  mutate(
    relcases=pmin(cases/max, 2)
  ) %>%
  merge(japan_pop_prefecture)

g1 <- ggplot(japan_rsv_prefecture2) +
  geom_raster(aes(year+(week-1)/52, reorder(prefecture, lat), fill=relcases)) +
  geom_vline(xintercept = 2013:2024, lty=2, col="white") +
  geom_segment(
    # Change the 0 and 25's to be appropriate to your data
    x = 2016.7, xend = 2016.7, y = 52, yend= Inf,
    arrow = arrow(length = unit(8, "pt")),
    col="#E02938",
    lwd=1
  ) +
  scale_x_continuous("Year", breaks=2013:2024, expand=c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c("Relative cases", 
                       breaks=c(0, 1, 2),
                       labels=c("0", "1", ">2")) +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(13.5, 5.5, 5.5, 5.5)
  )

japan_rsv_summ <- lapply(split(japan_rsv_prefecture, japan_rsv_prefecture$prefecture), function(x) {
  island <- x$island[1]
  
  x <- x %>%
    mutate(
      trueyear=year
    )
  
  if (island!="Ryukyu") {
    x <- x %>%
      mutate(
        week=week-26,
        year=ifelse(week <= 0, year-1, year),
        week=week%%52,
        week=ifelse(week==0, 52, week)
      )
  }
  
  out <- x %>%
    filter(year > 2012, year < 2020) %>%
    group_by(prefecture, island, year) %>%
    mutate(
      prop=(cases)/sum(cases)
    ) %>%
    summarize(
      cog=sum(week*cases)/sum(cases),
      intensity=1/(-sum(prop*log(prop), na.rm=TRUE)),
      size=sum(cases),
      trough=min(cases[trueyear<2020])
    ) %>%
    mutate(
      cog=ifelse(island!="Ryukyu", cog+26, cog)
    )
}) %>%
  bind_rows %>%
  merge(
    japan_pop_prefecture
  ) %>%
  mutate(
    season=paste0(year,"-",year+1),
    island=factor(island, levels=c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu"))
  )

g2 <- ggplot(japan_rsv_summ %>% filter(prefecture != "Okinawa")) +
  geom_boxplot(aes(season, cog), show.legend = FALSE, fill="#CCECF8") +
  scale_x_discrete("Season") +
  scale_y_continuous("Center of gravity (weeks)") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.text.x = element_text(hjust=1, angle=45)
  )

dd1 <- japan_rsv_summ %>%
  filter(prefecture != "Okinawa",
         season=="2015-2016")

dd2 <- japan_rsv_summ %>%
  filter(prefecture != "Okinawa",
         season=="2016-2017")

mean(dd1$cog-dd2$cog)

g3 <- ggplot(japan_rsv_summ) +
  geom_boxplot(aes(season, trough/pop), show.legend = FALSE, fill="#CCECF8") +
  scale_x_discrete("Season") +
  scale_y_continuous("Epidemic trough") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.text.x = element_text(hjust=1, angle=45)
  )

g4 <- ggplot(japan_rsv_island) +
  geom_vline(xintercept = 2013:2024, lty=3, col="gray80", lwd=0.5) +
  # annotate("rect", xmin=2016+0.4, xmax=2017+0.4, min=-Inf, ymax=Inf, alpha=0.2) +
  # geom_line(data=data_processed_travel3, aes(year+(month-1)/12, travel/pop*10), col="#E02938", lty=1) +
  geom_line(aes(year+(week-1)/52, cases), col="#224B95") +
  scale_x_continuous("Year", limits=c(2013, 2024.2), breaks=c(2013, 2015, 2017, 2019, 2021, 2023), expand=c(0, 0)) +
  scale_y_continuous("Cases", expand=c(0, 0)) +
  facet_wrap(~island, nrow=1, scale="free") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(size=1)
  )

g1a <- ggarrange(g1, labels="A")
g23 <- ggarrange(g2, g3, nrow=2, labels=c("B", "C"), byrow=FALSE)
g4a <- ggarrange(g4, labels=c("D"))

gcomb1 <- arrangeGrob(g1a, g23, nrow=1, widths=c(3, 1))
gcomb2 <- arrangeGrob(gcomb1, g4a, nrow=2, heights=c(3, 1))

ggsave("figure1.pdf", gcomb2, width=14, height=8)
ggsave("figure1.png", gcomb2, width=14, height=8)

g5 <- ggplot(japan_rsv_summ) +
  geom_boxplot(aes(season, cog), show.legend = FALSE, fill="#CCECF8") +
  scale_x_discrete("Season") +
  scale_y_continuous("Center of gravity (weeks)") +
  facet_wrap(~island, scale="free", nrow=1) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.text.x = element_text(hjust=1, angle=45),
    axis.title.x = element_blank(),
    strip.background = element_blank()
  )

ggsave("figure1_cog.pdf", g5, width=10, height=3)
