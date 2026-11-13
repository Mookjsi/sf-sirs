# Recorded settings, not newly validated

All values below are transcribed from supplied source. Actual historical per-fit settings and object fingerprints were not deserialized; the exact frozen object is authoritative. Initial values are existing inference initializers, not substitutes for fitted simulation states.

| Entry point / variant | Recorded settings |
|---|---|
| Current baseline driver | Hokkaido only; baseline; force_refit FALSE; 4 chains/cores; iter 16000; warmup 8000; seed 101; adapt_delta .99; max_treedepth 15 |
| Current no-npieff Hokkaido | 4 chains/cores; iter 12000; warmup 8000; seed 303; adapt_delta .95; max_treedepth 15 |
| Current no-npieff Shikoku/Ryukyu | 4 chains/cores; iter 4000; warmup 2000; seed 101; adapt_delta .95; max_treedepth 15 |
| Predictive observations | Seed 202; NB2 observations are distinct from deterministic expected-count trajectories |
| Current baseline prior horizon | N_train 260; N_total 364; Nnonpi Hokkaido 208, Shikoku 364, Ryukyu 312 |
| Humidity features | Training same-week reference and anomaly scale; K=3: z, z*sin(2*pi*week/52), z*cos(2*pi*week/52) |
| Original Honshu npi / npi2 | seed 102; adapt_delta .9; Nnonpi 364 |
| Original Kyushu npi / npi2 | seeds 101 / 102; adapt_delta .95 / .99; Nnonpi 364 |
| Original Ryukyu npi | seed 102; adapt_delta .99; Nnonpi 312 |
| Original Hokkaido npi / raw alternative | seed 102; adapt_delta .95; Nnonpi 208 / 364 |
| Original Shikoku npi / npi2 | seeds 104 / 102; adapt_delta .99 / .9; Nnonpi 364; tweak / tweak2 Stan |

Original npi2 time indexing uses `(year + week/52 > 2016.5)+1`. Baseline `np_run` offsets seeds when enumerating islands/models; do not reinterpret driver seed as proof of every stored fit's seed. Detailed chain/iteration settings for original scripts remain in the retained source. Tuning candidate settings, score weights, job plans, run manifests and diagnostics remain under their existing `code/perturbation/analysis/*tuning_run1/` paths. Alternative training windows and unused variants remain separate.

No value was retuned and no model was refit during cleanup.
