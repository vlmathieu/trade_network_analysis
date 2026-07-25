# In-house Resources

Data curated by the authors for this project (not downloaded from an external API).

---

## `correspondence_FAO_HS.json`

Lookup table mapping **FAO forest-product classes** to the **Harmonized System (HS)
commodity codes** used by UN Comtrade, across every HS revision (1996–2022). The
`join_fao_code` rule (`workflow/scripts/join_fao_code.py`) uses it to label each raw
trade record with its FAO product before aggregation.

**Format:** JSON array of 357 records, one per (FAO product × HS-code) row.

**Fields per record**

| Field | Description |
|-------|-------------|
| `FAO Product Agg` | Aggregated FAO product name, e.g. `WOOD IN THE ROUGH (ROUNDWOOD)`. |
| `FAO Code Agg` | Aggregated FAO division code — `01`, `05`, `07` (the three products analysed). |
| `FAO 1982` | Code under the older FAO 1982 classification. |
| `FAO Product` / `FAO Code` | Finer FAO product name / code (may be empty at the aggregated level). |
| `HS 1996` … `HS 2022` | HS commodity code for each revision year — the join key against Comtrade. |

> Source classification: FAO, 2022. *Classification of forest products 2022*. Rome.
> https://doi.org/10.4060/cb8216en
