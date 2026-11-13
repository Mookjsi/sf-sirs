# License and data-source notes

This is a local integration handoff for owner review. A file being retained here does not establish that it may be publicly redistributed. No online license review was performed.

## Software attribution

The original perturbation project includes an MIT License, copyright (c) 2025 Sang Woo Park. Preserve that LICENSE and copyright notice with copies or substantial portions of its software. The existing README identifies the research project as *Interplay between climate, childhood mixing, and population-level susceptibility explains a sudden shift in RSV seasonality in Japan* and cites Zenodo DOI 10.5281/zenodo.17396084. The manuscript also records DOI 10.5281/zenodo.17396085. These identifiers are retained as source provenance; their relationship and remote contents were not checked.

The separately supplied Desktop/analysis folder has no independent license statement identified. Its scripts and additions must remain identifiable as the collaborator contribution. The owner should establish authorship and redistribution authority before public integration. Do not assume the upstream MIT notice automatically licenses every collaborator addition, imported dataset, figure, or third-party component.

## Dataset attribution and review flags

| Retained input | Source stated in supplied manuscript | Redistribution status |
|---|---|---|
| data/japan_rsv.csv | NIID Infectious Diseases Weekly Report sentinel surveillance tables, beginning 2013 through end 2023 | Source attribution recorded; terms and extraction rights not independently established |
| data/japan_pop.csv | Statistics Bureau of Japan; 2022 prefecture population; https://www.e-stat.go.jp/en | Terms not independently established |
| data/japan_childcare.xlsx | Children and Families Agency | Terms not independently established |
| data_processed/data_processed_humidity.csv | ECMWF ERA5 specific humidity, aggregated from hourly approximately 31 km grids to daily/prefecture and weekly/island summaries | ERA5/Copernicus attribution and redistribution requirements require review |
| data_processed/data_processed_temperature.csv | ECMWF ERA5 surface air temperature, same aggregation description | ERA5/Copernicus attribution and redistribution requirements require review |

The local manuscript `doc/perturbation.tex` lines 320–338 is the source for these descriptions. It records the original RSV URL as `https://www.niid.go.jp/niid/en/survaillance-data-table-english.html` and childcare URL literally as `cfa.jp.gov`. These are unverified original citations; the latter must be checked rather than silently treated as a verified address. This handoff does not claim that the acquisition and redistribution terms were satisfied or that all original raw inputs are present.

## Figures, results, documents, and dependencies

Numerical result files and figures preserve the supplied record; their scientific origin and acceptance status are separate questions from copyright authorization. The owner should check authorship and publication permissions for manuscript PDFs, TeX sources, figures, and derivative analyses before public release. The source MIT notice is preserved but no new blanket license is asserted over the handoff contents.

Installed packages under .R-library are excluded from the primary package. Their third-party source and binary licenses are not replaced by the research project's MIT license. Reinstall dependencies under their own applicable terms; the package version inventory records the observed local environment and is not a lockfile or proof of portability.

Before public release, resolve data-provider attribution, collaborator licensing, and any manuscript/figure publication rights. These are owner review flags, not requests to execute or publish anything during cleanup.
