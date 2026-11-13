library(dplyr)
library(ggplot2)

project_root <- Sys.getenv("PERTURBATION_ROOT")

# 1. 각 섬의 저장 폴더
island_folders <- c(
  Hokkaido = "humidity_seasonal_hokkaido",
  Honshu   = "humidity_seasonal_standalone",
  Shikoku  = "humidity_seasonal_shikoku",
  Kyushu   = "humidity_seasonal_kyushu",
  Ryukyu   = "humidity_seasonal_ryukyu"
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
    xintercept = 52.5,
    linetype = "dotted"
  ) +
  facet_wrap(
    ~island,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = c(1, 27, 53, 79, 104),
    labels = c(
      "2018 W1", "2018 W27",
      "2019 W1", "2019 W27", "2019 W52"
    )
  ) +
  labs(
    title = "humidity_seasonal: five islands",
    subtitle = "Conditional hindcast, 2018-2019",
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
ggsave(
  file.path(
    project_root,
    "analysis",
    "humidity_seasonal_all_islands.png"
  ),
  g_seasonal_all,
  width = 12,
  height = 15,
  units = "in",
  dpi = 200,
  bg = "white"
)