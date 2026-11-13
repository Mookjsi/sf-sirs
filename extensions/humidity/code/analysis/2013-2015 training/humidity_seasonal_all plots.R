project_root <- Sys.getenv("PERTURBATION_ROOT")
Sys.setenv(PERTURBATION_ROOT = project_root)
source(file.path(project_root, "script", "setup.R"), encoding = "UTF-8")
perturbation_setup("analysis")
library(dplyr)
library(ggplot2)

train_years <- 2013:2015
test_years <- 2016:2019
combined_output_dir <- file.path(project_root, "analysis", "humidity_seasonal_all_train2013_2015")

# 1. 각 섬의 저장 폴더
island_folders <- c(
  Hokkaido = "humidity_seasonal_hokkaido_train2013_2015",
  Honshu   = "humidity_seasonal_honshu_train2013_2015",
  Shikoku  = "humidity_seasonal_shikoku_train2013_2015",
  Kyushu   = "humidity_seasonal_kyushu_train2013_2015",
  Ryukyu   = "humidity_seasonal_ryukyu_train2013_2015"
)

result_files <- file.path(
  project_root,
  "analysis",
  island_folders,
  "humidity_seasonal_result.rds"
)
names(result_files) <- names(island_folders)

# 2. 저장 파일 확인
if (!all(file.exists(result_files))) {
  stop(
    "Missing results: ",
    paste(
      names(result_files)[!file.exists(result_files)],
      collapse = ", "
    )
  )
}

# 3. 기존 예측 곡선 불러오기
curves <- lapply(result_files, function(path) {
  
  res <- readRDS(path)

  stopifnot(
    res$standata$N == 52L * length(train_years),
    identical(as.integer(res$train$year), rep(train_years, each = 52L)),
    identical(as.integer(res$test$year), rep(test_years, each = 52L)),
    identical(as.integer(res$test$week), rep(1:52, length(test_years))),
    identical(as.integer(res$result$curve$year), rep(test_years, each = 52L)),
    identical(as.integer(res$result$curve$step), seq_len(52L * length(test_years)))
  )
  
  if (!isTRUE(res$diagnostics$summary$diagnostics_ok)) {
    stop("MCMC diagnostics failed: ", path)
  }
  
  res$result$curve
})

prediction_all <- bind_rows(curves, .id = "island") %>%
  mutate(
    island = factor(island, levels = names(island_folders))
  ) %>%
  arrange(island, step)

# 4. 다섯 섬의 예측 그림
g_seasonal_all <- ggplot(
  prediction_all,
  aes(x = step)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#90CAF9",
    alpha = 0.45
  ) +
  geom_line(
    aes(y = observed),
    colour = "black",
    linewidth = 0.45
  ) +
  geom_line(
    aes(y = pred_median),
    colour = "#0072B2",
    linewidth = 0.65
  ) +
  geom_vline(
    xintercept = c(52.5, 104.5, 156.5),
    linetype = "dotted"
  ) +
  facet_wrap(
    ~island,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = c(1, 53, 105, 157, 208),
    labels = c(
      "2016 W1", "2017 W1",
      "2018 W1", "2019 W1", "2019 W52"
    )
  ) +
  labs(
    title = "humidity_seasonal: five islands",
    subtitle = "Conditional hindcast, 2016-2019 (trained on 2013-2015)",
    x = NULL,
    y = "Weekly reported cases",
    caption = paste(
      "Black: observed; blue: predictive median;",
      "shading: 95% predictive interval."
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(g_seasonal_all)

# 5. 큰 이미지로 저장
dir.create(combined_output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(
    combined_output_dir,
    "humidity_seasonal_all_islands_2016_2019.png"
  ),
  g_seasonal_all,
  width = 16,
  height = 15,
  units = "in",
  dpi = 200,
  bg = "white"
)
message("Saved 2016-2019 prediction plot in: ", combined_output_dir)
print(g_seasonal_all)
