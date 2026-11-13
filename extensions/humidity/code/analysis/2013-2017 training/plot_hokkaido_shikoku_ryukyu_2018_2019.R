# Three restored seasonal-humidity models WITHOUT npieff. No MCMC here.
g_three_islands_2013_2017 <- NULL
project_root <- Sys.getenv("PERTURBATION_ROOT")
support_dir <- paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/npieff_support")
humidity_script_dir <- paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training")
.libPaths(c(file.path(project_root, ".R-library"), .libPaths()))
library(ggplot2)
source(file.path(support_dir, "npieff_workflow.R"), encoding = "UTF-8")
source(file.path(support_dir, "humidity_no_npieff_support.R"), encoding = "UTF-8")
islands <- c("Hokkaido", "Shikoku", "Ryukyu")
checked <- lapply(islands, function(island) hn_collect(project_root, island, humidity_script_dir))
prediction_three_islands <- do.call(rbind, lapply(checked, function(x) x$curves))
diagnostics_three_islands <- do.call(rbind, lapply(checked, function(x) x$diagnostics))
prediction_three_islands$island <- factor(prediction_three_islands$island, levels = islands)
g_three_islands_2013_2017 <- ggplot(prediction_three_islands, aes(x = step)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#90CAF9", alpha = 0.45) +
  geom_line(aes(y = observed), colour = "black", linewidth = .35) +
  geom_line(aes(y = pred_median), colour = "#0072B2", linewidth = .55) +
  geom_vline(xintercept = 52.5, linetype = "dotted") +
  facet_grid(island ~ ., scales = "free_y") +
  scale_x_continuous(breaks = c(1,27,53,79,104),
                    labels = c("2018 W1","W27","2019 W1","W27","W52")) +
  labs(title = "Seasonal humidity without npieff: three islands",
       subtitle = "Training: 2013-2017; conditional hindcast: 2018-2019",
       x = NULL, y = "Weekly reported cases",
       caption = "Black: observed; blue: predictive median; shading: 95% predictive interval.") +
  theme_bw()
print(diagnostics_three_islands, row.names = FALSE, digits = 4)
print(g_three_islands_2013_2017)
output_dir <- getOption("humidity_no_npieff.plot_output_dir",
  paste0(Sys.getenv("HANDOFF_ANALYSIS_ROOT"), "/2013-2017 training/plots_humidity_seasonal_no_npieff"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "humidity_seasonal_no_npieff_three_islands_2018_2019.png"),
       plot = g_three_islands_2013_2017, width = 12, height = 10, units = "in", dpi = 200)
message("Saved no-npieff humidity figure in: ", output_dir)

