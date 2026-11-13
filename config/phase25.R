phase25_config <- function() list(
  version="phase25-1",
  phase24_report="artifacts/phase24-report-002/result.rds",
  phase24_adjudicate="artifacts/phase24-adjudicate-002/result.rds",
  phase24_dir="artifacts/phase24-report-002",
  replay_tolerance=1e-7,
  padding_weeks=13L,
  minimum_bisections=5L,
  transition_width=1e-6,
  stable_rounds=2L,
  maximum_bisections=20L,
  identity=list(match_distance=3,max_log_height_change=log(2),
    height_weight=2,context_weight=1,margin=.5),
  early_max_week=19L,late_min_week=20L)
