phase28_config <- function() list(
  version="phase28-1",region="Honshu",
  park_commit="475c7dcdadc12987033d332b101f025625074301",
  park_stan="references/parksw3/stanmodel/sirs_npi.stan",
  park_driver="references/parksw3/stanfit_sirs/stanfit_sirs_honshu_npi.R",
  park_preprocessing="references/parksw3/script/script_data.R",
  park_license="references/parksw3/LICENSE",
  original_cases="data/processed/original-island-v1/japan_rsv_island.csv",
  original_population="data/processed/original-island-v1/japan_pop_island.csv",
  park_seed=102L,park_chains=4L,park_cores=4L,park_iter=4000L,
  park_treedepth=15L,park_adapt_delta=0.9,
  rhat_max=1.01,ess_min=400,bfmi_min=0.3,
  m2_fit="artifacts/phase26-m2-fit-001/result.rds",
  m2_evaluation="artifacts/phase26-m2-evaluate-001/result.rds")
