# Configuration

Settings and dependency pins for the trade-network-analysis pipeline. Edit these
files to change *what* the workflow computes without touching any code.

---

## `config.yaml`

Central parameter file read by `workflow/Snakefile` and every rule. A new user
can retune the whole analysis from here.

**Key parameters**

| Parameter | Meaning |
|-----------|---------|
| `years.start` / `years.stop` | Time span. ⚠ `stop` is **exclusive by 2** — `stop: 2025` processes up to **2023**. |
| `HS_version` | Harmonized System revisions to map (classifications change over time). |
| `flowCode` | Trade flows kept: `M` (imports), `X` (exports). |
| `fao_divisions` | Fine FAO product codes downloaded. |
| `fao_divisions_agg` | Aggregated divisions analysed: `01` Roundwood, `05` Sawnwood, `07` Wood-based panels. |
| `wb_series` | World Bank unit-value-index series used to deflate values to constant USD 2015. |
| `excluded_iso` | Regex patterns dropping non-country reporters (World, "…, nes", Free Zones, etc.). |
| `col_keep` | Columns retained from raw Comtrade to lighten computation. |
| `eu.desc` / `eu.iso` | EU members merged into one "EU" node in the `agg_eu` level (intra-EU flows dropped). |
| `agg_lvl` | The two aggregation levels: `country_lvl` (every country) and `agg_eu` (EU as one node). |
| `weight` | Edge weights computed: `primary_value` and `net_wgt`. |
| `threshold_main_contributors` | Min share (0.01) to count as a main contributor. |
| `reference_years` | Reference window (2000–2020) shaded across figures. |
| `network_contribution.top_frac` | **Upper-tail probability, not a floor** — keeps countries above the `1 - top_frac` quantile over the last `time_span` years (0.03 = 97th pctile). |
| `chord` / `trade_network` | Year, top-N, and node threshold for the `plot_trade_flows` diagrams. |
| `figure_ext` | Output formats for all figures (`png`, `svg`). |

---

## `requirements.txt`

Pinned versions of every Python and R package the pipeline depends on (polars,
networkx, scipy, wbgapi; r-tidyverse, r-ggplot2, r-patchwork, …). Reference
manifest — the actual per-rule environments are built from `workflow/envs/*.yaml`
via `snakemake --sdm conda`.
