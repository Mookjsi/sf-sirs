# Reproduce figure1 (2).png from the Park et al. humidity CSV.
# Run this whole file with Source. No earlier scripts or MCMC are needed.
# Objects left in R: humidity_weekly, humidity_annual, humidity_period, g_humidity.

# 1. Paths and package ---------------------------------------------------------
project_root <- Sys.getenv("PERTURBATION_ROOT")
output_dir <- paste0(Sys.getenv("HANDOFF_ROOT"), "/generated/rsv_humidity_analysis_summary_assets")
output_file <- file.path(output_dir, "figure1_reproduced.png")
save_plot <- TRUE

# Use the packages already installed for this project, when available.
project_library <- file.path(project_root, ".R-library")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop('Install ggplot2 first: install.packages("ggplot2")', call. = FALSE)
}
library(ggplot2)

# 2. Read and match all 47 prefectures -----------------------------------------
humidity_file <- file.path(project_root, "data_processed", "data_processed_humidity.csv")
population_file <- file.path(project_root, "data", "japan_pop.csv")
stopifnot(file.exists(humidity_file), file.exists(population_file))

# encoding="UTF-8" avoids truncating names on a Windows native locale.
q_raw <- read.csv(humidity_file, encoding = "UTF-8", stringsAsFactors = FALSE)
prefecture_map <- read.csv(population_file, encoding = "UTF-8", stringsAsFactors = FALSE)
stopifnot(all(c("key", "year", "week", "humidity") %in% names(q_raw)),
          all(c("prefecture", "island") %in% names(prefecture_map)))

prefecture_map$prefecture <- sub("-.*", "", prefecture_map$prefecture)
q_raw$prefecture <- trimws(q_raw$key)
q_raw$prefecture[q_raw$prefecture == "Hy\u014dgo"] <- "Hyogo"
q_raw$prefecture[q_raw$prefecture == "Naoasaki"] <- "Nagasaki"

island_order <- c("Hokkaido", "Honshu", "Shikoku", "Kyushu", "Ryukyu")
expected_counts <- c(Hokkaido = 1L, Honshu = 34L, Shikoku = 4L,
                     Kyushu = 7L, Ryukyu = 1L)
map_counts <- table(factor(prefecture_map$island, levels = island_order))
stopifnot(nrow(prefecture_map) == 47L,
          !anyDuplicated(prefecture_map$prefecture),
          !anyNA(prefecture_map$island),
          all(prefecture_map$island %in% island_order),
          identical(as.integer(map_counts), unname(expected_counts)))

q_raw <- q_raw[q_raw$year %in% 2013:2019, ]
q_raw$island <- prefecture_map$island[match(q_raw$prefecture, prefecture_map$prefecture)]
stopifnot(setequal(unique(q_raw$prefecture), prefecture_map$prefecture),
          !anyNA(q_raw[c("prefecture", "island", "year", "week", "humidity")]),
          all(is.finite(q_raw$humidity)), all(q_raw$humidity > 0),
          all(q_raw$week %in% 1:52),
          !anyDuplicated(q_raw[c("prefecture", "year", "week")]))
coverage <- aggregate(humidity ~ island + year + week, q_raw, length)
stopifnot(nrow(coverage) == 5L * 7L * 52L,
          all(coverage$humidity == expected_counts[coverage$island]))

# 3. Weekly island means -> annual means -> period means -----------------------
# All prefectures within an island receive equal weight (not population weights).
# The source CSV is in kg/kg; multiply by 1000 to display g/kg.
q_raw$q_gkg <- q_raw$humidity * 1000
humidity_weekly <- aggregate(q_gkg ~ island + year + week, q_raw, mean)
humidity_annual <- aggregate(q_gkg ~ island + year, humidity_weekly, mean)
annual_weeks <- aggregate(q_gkg ~ island + year, humidity_weekly, length)
stopifnot(nrow(humidity_annual) == 35L, all(annual_weeks$q_gkg == 52L))

period_order <- c("2013-2015", "2016-2019")
humidity_annual$island <- factor(humidity_annual$island, levels = island_order)
humidity_annual$period <- factor(
  ifelse(humidity_annual$year <= 2015, "2013-2015", "2016-2019"),
  levels = period_order
)
humidity_annual <- humidity_annual[order(humidity_annual$island, humidity_annual$year), ]
rownames(humidity_annual) <- NULL

# Each year has 52 weeks, so averaging annual means equals averaging all weeks.
humidity_period <- aggregate(q_gkg ~ island + period, humidity_annual, mean)
humidity_period$x_start <- ifelse(humidity_period$period == "2013-2015", 2013, 2015.5)
humidity_period$x_end <- ifelse(humidity_period$period == "2013-2015", 2015.5, 2019)

# 4. Draw the five panels ------------------------------------------------------
g_humidity <- ggplot(humidity_annual, aes(x = year, y = q_gkg)) +
  geom_vline(xintercept = 2015.5, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_line(aes(group = island), colour = "grey40", linewidth = 0.5) +
  geom_segment(
    data = humidity_period,
    aes(x = x_start, xend = x_end, y = q_gkg, yend = q_gkg, colour = period),
    inherit.aes = FALSE, linetype = "dotted", linewidth = 0.7
  ) +
  geom_point(aes(colour = period), size = 2.8) +
  facet_wrap(~ island, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("2013-2015" = "#D62728", "2016-2019" = "#0072B2"),
                      breaks = period_order, name = "Period") +
  scale_x_continuous(breaks = 2013:2019, minor_breaks = NULL) +
  labs(
    title = "Specific humidity before and after 2016",
    subtitle = "Annual means, 2013-2019",
    x = "Year", y = "Specific humidity (g/kg)",
    caption = paste0(
      "Points: mean of 52 weeks; dotted horizontal lines: period means.\n",
      "Vertical dashed line separates 2015 and 2016; y-axis scales differ across islands."
    )
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", hjust = 0),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0),
    plot.title.position = "panel"
  )

print(g_humidity)

# 5. Save at the same pixel dimensions as the supplied figure ------------------
# The supplied image is preserved; the new PNG has a different filename.
if (save_plot) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename = output_file, width = 1600, height = 1760,
                units = "px", res = 160,
                type = if (.Platform$OS.type == "windows") "windows" else "cairo",
                bg = "white")
  tryCatch(print(g_humidity), finally = grDevices::dev.off())
  message("Saved: ", normalizePath(output_file, winslash = "/"))
}

# To display it again later in this R session: print(g_humidity)
# To inspect the underlying values: print(humidity_annual)
