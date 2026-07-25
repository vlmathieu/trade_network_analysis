# Conda Environments

One isolated environment per tool stack. Snakemake builds these automatically on the
first `snakemake --sdm conda` run and activates the right one for each rule (via the
`conda:` directive in `workflow/rules/*.smk`). Pinned package versions live in
`config/requirements.txt`.

| File | Stack | Rules that use it |
|------|-------|-------------------|
| `comtradeapicall.yaml` | `comtradeapicall` client | `get_uncomtrade` |
| `wbgapi.yaml` | World Bank `wbgapi` client | `wb_data` |
| `polars.yaml` | polars data wrangling | Stage-1 transforms (`correspondence_FAO_HS`, `deflate_uncomtrade`, `join_fao_code`, `filter_data`, `aggregate_eu`) + Stage-2 `network_objects` |
| `network_metrics.yaml` | networkx + scipy | metric rules `network_composition`, `market_concentration`, `network_contribution`, `contributor_profiles` |
| `network_connectivity.yaml` | connectivity/graph deps | `network_connectivity` (kept separate — heavier graph algorithms) |
| `r_plots.yaml` | R + tidyverse/ggplot2 | all `plot_*` figure rules except trade flows |
| `r_trade_flows.yaml` | R + chord/network-graph deps | `plot_trade_flows` (chord + node-link diagrams) |

> Why split? Isolating environments keeps dependency conflicts out (e.g. the R
> plotting stack never has to coexist with polars), and lets Snakemake rebuild only
> what changed.
