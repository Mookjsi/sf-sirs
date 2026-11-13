# 1. 각 사후표본에서 정점 주차를 찾는 함수

library(dplyr)
library(ggplot2)

# B: 2013–2015년 같은 주차의 평균 비습 적용
z_B_seasonal <- rep(0, nrow(honshu_test_seasonal))

# 계절 반응에 필요한 세 입력 만들기
phase_B <- 2 * pi * honshu_test_seasonal$week / 52

X_B_seasonal <- cbind(
  humidity     = z_B_seasonal,
  humidity_sin = z_B_seasonal * sin(phase_B),
  humidity_cos = z_B_seasonal * cos(phase_B)
)

# 동일한 사후표본으로 B의 기대 신고 수 계산
C_humidity_seasonal_B <- hs_predict(
  draws = draws_humidity_seasonal,
  week = honshu_test_seasonal$week,
  standata = standata_humidity_seasonal,
  X_test = X_B_seasonal
)

stopifnot(
  identical(
    dim(C_humidity_seasonal),
    dim(C_humidity_seasonal_B)
  ),
  all(is.finite(C_humidity_seasonal_B)),
  all(C_humidity_seasonal_B >= 0)
)

# 1. 각 사후표본에서 정점 주차를 찾는 함수
peak_info <- function(M, weeks) {
  
  k <- apply(M, 1, which.max)
  maximum <- M[cbind(seq_len(nrow(M)), k)]
  
  data.frame(
    peak_week = weeks[k],
    
    # 정점이 해당 연도의 첫 주 또는 마지막 주인지
    boundary = k == 1L | k == ncol(M),
    
    # 동일한 최대값을 가진 주가 여러 개인지
    tied = rowSums(M == maximum) > 1L
  )
}

# 2. 연도별로 A와 B의 같은 사후표본끼리 비교
peak_comparison <- bind_rows(
  lapply(sort(unique(honshu_test_seasonal$year)), function(yr) {
    
    idx <- which(honshu_test_seasonal$year == yr)
    weeks <- honshu_test_seasonal$week[idx]
    
    stopifnot(identical(as.integer(weeks), 1:52))
    
    A <- peak_info(
      C_humidity_seasonal[, idx, drop = FALSE],
      weeks
    )
    
    B <- peak_info(
      C_humidity_seasonal_B[, idx, drop = FALSE],
      weeks
    )
    
    data.frame(
      year = yr,
      draw = seq_len(nrow(A)),
      peak_week_A = A$peak_week,
      peak_week_B = B$peak_week,
      shift_weeks = A$peak_week - B$peak_week,
      boundary = A$boundary | B$boundary,
      tied = A$tied | B$tied
    )
  })
)

# 3. 사후표본별 차이를 먼저 구한 뒤 요약
peak_summary <- peak_comparison %>%
  group_by(year) %>%
  summarise(
    median_shift_weeks = median(shift_weeks),
    lower_95 = quantile(shift_weeks, 0.025),
    upper_95 = quantile(shift_weeks, 0.975),
    prob_earlier = mean(shift_weeks < 0),
    prob_same = mean(shift_weeks == 0),
    prob_later = mean(shift_weeks > 0),
    .groups = "drop"
  )

print(peak_summary, width = Inf)

# 4. 연도 경계와 정점 동률 여부 확인
peak_checks <- peak_comparison %>%
  group_by(year) %>%
  summarise(
    boundary_draws = sum(boundary),
    tied_draws = sum(tied),
    .groups = "drop"
  )

print(peak_checks)

print(peak_summary, width = Inf)

library(dplyr)
library(ggplot2)

# 연도별로 정점 차이의 비율 계산
peak_distribution <- peak_comparison %>%
  count(year, shift_weeks, name = "n") %>%
  group_by(year) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    direction = factor(
      case_when(
        shift_weeks < 0 ~ "Earlier",
        shift_weeks == 0 ~ "Same week",
        shift_weeks > 0 ~ "Later"
      ),
      levels = c("Earlier", "Same week", "Later")
    )
  )

# 정점 시기 차이의 분포
g_peak_shift <- ggplot(
  peak_distribution,
  aes(x = shift_weeks, y = proportion, fill = direction)
) +
  geom_col(width = 0.85) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  facet_wrap(~year, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Earlier" = "#0072B2",
      "Same week" = "#999999",
      "Later" = "#D55E00"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(
      min(peak_distribution$shift_weeks),
      max(peak_distribution$shift_weeks),
      by = 1
    )
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Peak timing difference: humidity_seasonal",
    subtitle = "A: actual humidity; B: 2013-2015 weekly-mean humidity",
    x = "Peak timing difference, A - B (weeks)",
    y = "Posterior proportion",
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(g_peak_shift)
