# Scripts

The actual logic behind each rule. **Python (`*.py`)** does data wrangling and metric
computation; **R (`*.R`)** does the figures. Scripts are invoked by the matching rule
in `workflow/rules/` — run them through Snakemake, not by hand (they expect
Snakemake-provided inputs, outputs, and params).

## Python — data & metrics

### Stage 1 — build the trade dataset (polars)

| Script | Role |
|--------|------|
| `correspondence_FAO_HS.py` | Load/prepare the FAO↔HS correspondence table. |
| `get_uncomtrade_data.py` | Download raw bilateral trade from UN Comtrade. |
| `wb_data.py` | Download World Bank unit-value indices. |
| `deflate_uncomtrade.py` | Deflate values to constant USD 2015. |
| `join_fao_code.py` | Attach FAO product labels to each record. |
| `filter_data.py` | Filter divisions/columns; guarded sums keep all-null groups null (data-quality signal preserved). |
| `aggregate_eu.py` | Merge EU members into one node for the `agg_eu` level. |

### Stage 2

| Script | Role |
|--------|------|
| `network_objects.py` | Build directed edge lists + mirror-flows table; the two mirror sides (`_exp`/`_imp`) are kept separate, never reconciled. |

### Stage 3 — metrics

| Script | Computes |
|--------|----------|
| `network_composition.py` | Classifies each country main-exporter / main-importer / balanced (80% own-reported value-share threshold). |
| `market_concentration.py` | Concentration indices (two rows per year, one per edge weight). |
| `network_connectivity.py` | Connectivity / graph-structure metrics. |
| `network_contribution.py` | Per-country contribution shares (sum to 200% by design — each edge touches 2 nodes). |
| `contributor_profiles.py` | Trader-profile time series for the dominant contributors. |

## R — figures

Each `plot_*.R` reads a Stage-3 CSV and writes `.png` + `.svg` to `results/`.

| Script | Figure |
|--------|--------|
| `plot_network_composition.R` | Trader-type composition over time. |
| `plot_market_concentration.R` | Market-concentration trends. |
| `plot_network_connectivity.R` | Connectivity metrics. |
| `plot_network_contribution.R` | Country contribution (mean line = exporter/importer reports; band = the discrepancy). |
| `plot_contributor_profiles.R` | Dominant-contributor profiles. |
| `plot_prices_figures.R` | Price descriptive stats + flow-level ridgelines. |
| `plot_trade_flows.R` | Chord diagram + node-link trade-network graph. |

> **Shared plot helpers live in `utils.R`** (palettes `PAL_EXP`/`PAL_IMP`/`PAL_BALANCED`/`PAL_TOTAL`, `theme_trade`, `benchmark_lines`, `save_figure`), sourced by every `plot_*.R`.
> ⚠ `utils.R` is **not yet committed to git** — a fresh clone will be missing it and all plot rules will fail. Commit it before pushing so the pipeline is reproducible for a new user.
