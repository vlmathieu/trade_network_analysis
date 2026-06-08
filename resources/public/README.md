# Public Datasets

This folder contains the two raw input datasets for the trade network analysis pipeline.

---

## `uncomtrade_data.parquet.gzip`

**Source:** [UN Comtrade](https://comtradeplus.un.org/) — downloaded via the `comtradeapicall` Python library (`workflow/scripts/get_uncomtrade_data.py`)

**Size:** ~41 MB (compressed), ~2.05 million records

**Coverage**
| Dimension | Scope |
|-----------|-------|
| Time | 1996–2023, annual |
| Geography | Global bilateral flows (all UN Comtrade reporters and partners) |
| Products | Wood and wood products — FAO (2022) divisions 01 (Wood in the rough), 05 (Sawnwood), 07 (Wood-based panels) |
| Flow types | Exports (`X`) and Imports (`M`) |

> FAO, 2022. *Classification of forest products 2022*. Rome. https://doi.org/10.4060/cb8216en — correspondence tables between FAO divisions and HS codes are provided in this reference and in `resources/correspondence_FAO_HS.json`.

**Key columns**

| Column | Description |
|--------|-------------|
| `reporterISO` | ISO 3-letter code of the reporting country |
| `reporterDesc` | Name of the reporting country |
| `partnerISO` | ISO 3-letter code of the partner country |
| `partnerDesc` | Name of the partner country |
| `flowCode` | Trade flow direction: `X` (export) or `M` (import) |
| `flowDesc` | Trade flow description: `Export` or `Import` |
| `cmdCode` | HS product code |
| `classificationCode` | HS revision (H0, H1, H2, …) |
| `period` | Reference year (string), used as join key in processing scripts |
| `netWgt` | Net weight in kilograms |
| `primaryValue` | Trade value in current USD |

**Pipeline role**

This is the primary raw data source. It flows through three processing steps before analysis:
1. `workflow/scripts/deflate_uncomtrade.py` — adds `primaryValue_deflated` using World Bank UVI indices
2. `workflow/scripts/join_fao_code.py` — joins FAO↔HS code correspondence table, adding `FAO Code`, `FAO Code Agg`, `FAO Product`, and `FAO Product Agg` columns
3. `workflow/scripts/filter_data.py` — selects and filters columns for network analysis

---

## `wb_series_data.csv`

**Source:** [World Bank](https://data.worldbank.org/) — downloaded via the `wbgapi` Python library (`workflow/scripts/wb_data.py`)

**Size:** ~148 KB, 7,715 rows

**Coverage**
| Dimension | Scope |
|-----------|-------|
| Time | 1996–2024 |
| Geography | ~215 economies (ISO 3-letter codes, including regional aggregates) |

**Columns**

| Column | Description |
|--------|-------------|
| `economy` | ISO 3-letter economy code (e.g., `CHN`, `AFE`) |
| `time` | Year in `YR####` format (e.g., `YR2015`) |
| `TM.UVI.MRCH.XD.WD` | Import Unit Value Index, merchandise trade (base year 2015 = 100) |
| `TX.UVI.MRCH.XD.WD` | Export Unit Value Index, merchandise trade (base year 2015 = 100) |

Missing values are common in earlier years and for economies with limited trade reporting.

**Pipeline role**

Used in `workflow/scripts/deflate_uncomtrade.py` to convert UN Comtrade trade values from current USD to constant 2015 USD, enabling consistent comparison across years.
