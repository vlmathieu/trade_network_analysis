from snakemake.script import snakemake
import logging
import polars as pl
import re

def format_col_names(col_names: list) -> dict:
    '''
    Function that returns a mapping between column names and their formated 
    version as a dictionnary.

    Parameters
    ----------
    col_names : list of strings
        The list of column names to format.

    Returns
    -------
    col_mapping : dictionnary
        Mapping between column names and their formated version.
    '''

    # Insert an underscore before a uppercase
    col_names_formated = [re.sub(r"([a-z])([A-Z])", "\\1_\\2", s) for s in col_names]

    # Replace space by underscore
    col_names_formated = [re.sub(" ", "_", s) for s in col_names_formated]

    # Put every column name to lowercase
    col_names_formated = [s.lower() for s in col_names_formated]

    # Create mapping between column names and their formated version
    col_mapping = dict(zip(col_names, col_names_formated))

    return col_mapping

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load data
# Read the parquet handoff from join_fao_code. Parquet preserves dtypes, so no
# schema_overrides are needed: 'period' is already String (the deflate_uncomtrade
# join on 'period' requires it) and 'FAO Code Agg' keeps the zero-padded codes
# '01'/'05'/'07' as String, so the str(...) comparisons below still work.
merged_data = pl.read_parquet(snakemake.input[0])

# Filter data for network analysis
input_data = (
    merged_data
        .select(snakemake.params['col_keep'])
        # Format column names for homogenity
        .rename(format_col_names(snakemake.params['col_keep']))
        .filter(
            # Keep data in specified time range
            pl.col('period') >= str(snakemake.params['year_start']),
            pl.col('period') <= str(snakemake.params['year_stop']-2),

            # Remove non-country reporters and partners
            (~pl.col('reporter_iso')
             .str.contains('|'.join(snakemake.params['excluded_iso']))),
            (~pl.col('partner_iso')
             .str.contains('|'.join(snakemake.params['excluded_iso']))),

            # Delete "auto-" imports or exports (reporter = partner)
            pl.col('reporter_desc') != pl.col('partner_desc'),
        )
        # Drop potential duplicates
        .unique()
)

# Drop outliers = raw trade flows under the fifth percentile for weight (kg) and
# value (USD), computed PER aggregated FAO product (fao_product_agg) so the very
# different magnitude scales of divisions 01/05/07 are not pooled into one global
# threshold. Done on the RAW flows, before aggregation to FAO product.
thresholds = (
    input_data
    .group_by('fao_product_agg')
    .agg(
        # quantile ignores nulls
        pl.col('net_wgt').quantile(0.05).alias('min_weight'),
        pl.col('primary_value').quantile(0.05).alias('min_value'),
    )
    .sort('fao_product_agg')
)
logging.info(f"\nPer-product 5th-percentile thresholds:\n {thresholds}\n")

input_data = (
    input_data
    .join(thresholds, on='fao_product_agg', how='left')
    .filter(
        # Drop trade flow with net weight (kg) under its product's fifth percentile.
        # The is_null() pass keeps unreported weights (e.g. the 2000-2006 quantity
        # gaps), so only reported-but-implausibly-small flows are removed.
        ((pl.col('net_wgt') > pl.col('min_weight')) |
         (pl.col('net_wgt').is_null())),

        # Drop trade flow with value (USD) under its product's fifth percentile
        # (same null-passing guard as the weight condition).
        ((pl.col('primary_value') > pl.col('min_value')) |
         (pl.col('primary_value').is_null())),
    )
    .drop(['min_weight', 'min_value'])
)

# Sum weight and values by FAO division product.
# Guarded sum (same pattern as aggregate_eu.py): a group whose values are ALL
# null aggregates to null, not 0. This keeps the released data honest — a missing
# net-weight report (e.g. the 2000-2006 quantity gaps) stays null instead of being
# turned into a fake reported 0. Groups with at least one reported value sum
# exactly as a plain sum would, so downstream results are unchanged (every metric
# consumer already treats null/None as 0).
sum_cols = ['net_wgt', 'primary_value', 'primary_value_deflated']

groupby_cols = [_ for _ in input_data.columns if _ not in sum_cols]

input_data = (
    input_data
    .group_by(groupby_cols)
    .agg(pl.when(pl.col(sum_cols).count() > 0).then(pl.sum(sum_cols)))
 )
logging.info(f"\nFinal data:\n {input_data}\n")

# Final descriptive statistics of the saved dataset, per FAO product
# (5th / 10th / 50th / 90th / 95th percentiles for weight, value and deflated value)
with pl.Config(tbl_rows=-1, tbl_cols=-1):
    for cmd in sorted(input_data.select('fao_code_agg').unique().to_series().to_list()):
        desc = (
            input_data
            .filter(pl.col('fao_code_agg') == cmd)
            .select(['net_wgt', 'primary_value', 'primary_value_deflated'])
            .describe(percentiles=[0.05, 0.10, 0.50, 0.90, 0.95])
        )
        logging.info(f"\nFinal descriptive stats - FAO product {cmd}:\n {desc}\n")

# Save input data
input_data.write_parquet(
    snakemake.output[0],
    compression='gzip'
    )
