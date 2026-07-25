# Rules

One [Snakemake rule](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
per `.smk` file, loaded by `workflow/Snakefile`. Each declares a step's **inputs,
outputs, params, log, benchmark, and conda environment**, then calls the matching
script in `workflow/scripts/`. Rules are the "wiring"; the scripts are the logic.

Run one with `snakemake <rule> --sdm conda` (drop the `.smk`).

### Stage 1 — build the trade dataset (polars, sequential)

| Rule | Does |
|------|------|
| `correspondence_FAO_HS` | Load the FAO↔HS correspondence table. |
| `get_uncomtrade` | Download raw bilateral trade records from UN Comtrade (needs `comtrade_apikey`). |
| `wb_data` | Download World Bank unit-value indices for deflation. |
| `deflate_uncomtrade` | Convert trade values to constant USD 2015. |
| `join_fao_code` | Tag each record with its FAO product via the correspondence table. |
| `filter_data` | Keep the target divisions/columns, drop non-country reporters. |
| `aggregate_eu` | Produce the `agg_eu` level (EU as one node, intra-EU flows dropped). |

### Stage 2 — build network objects

| Rule | Does |
|------|------|
| `network_objects` | Emit `mirror_flows.csv` + `edge_lists.pkl` (directed edges, mirror sides kept separate). |

### Stage 3 — metrics (compute) + figures (plots)

Each metric has a Python compute rule and a matching R plot rule; both run at both
aggregation levels.

| Metric | Compute rule | Plot rule |
|--------|--------------|-----------|
| Network composition (trader typing) | `network_composition` | `plot_network_composition` |
| Market concentration | `market_concentration` | `plot_market_concentration` |
| Network connectivity | `network_connectivity` | `plot_network_connectivity` |
| Network contribution | `network_contribution` | `plot_network_contribution` |
| Contributor profiles | `contributor_profiles` | `plot_contributor_profiles` |

Plot-only rules (no separate metric file):

| Rule | Produces |
|------|----------|
| `plot_prices_figures` | Price descriptive stats + flow-level ridgeline figures. |
| `plot_trade_flows` | Chord diagram + node-link trade-network graph. |
