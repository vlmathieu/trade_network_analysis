from snakemake.script import snakemake
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

# Load data
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

            # Keep imports and exports only
            pl.col('flow_code').is_in(snakemake.params['flow_to_keep']),

            # Keep FAO product division specified
            pl.col('fao_code').is_in(snakemake.params['fao_divisions']),

            # Remove non-country reporters and partners
            (~pl.col('reporter_iso')
             .str.contains('|'.join(snakemake.params['excluded_iso']))),
            (~pl.col('partner_iso')
             .str.contains('|'.join(snakemake.params['excluded_iso']))),

            # Delete "auto-" imports or exports (reporter = partner)
            pl.col('reporter_desc') != pl.col('partner_desc'),

            # # Remove primaryValue_deflated null values (~8% of the dataset)
            # ~pl.col('primary_value_deflated').is_null()
        )
        # Drop potential duplicates
        .unique()
)

# Drop outliers = values under fifth percentile for weight (kg) and value (USD)
stats_desc = (
    input_data
    .select(['net_wgt', 'primary_value'])
    .describe(percentiles=[0.05])
)

min_weight, min_value = (
    stats_desc
    .filter(pl.col('statistic') == '5%')
    .select(['net_wgt', 'primary_value'])
)

input_data = (
    input_data.filter(
        # Drop trade flow with net weight (kg) under fifth percentile
        pl.col('net_wgt') > min_weight.item(),

        # Drop trade flow with value (USD) under fifth percentile
        pl.col('primary_value') > min_value.item()
    )
)

# Sum weight and values by FAO division product
sum_cols = ['net_wgt', 'primary_value', 'primary_value_deflated']

groupby_cols = [_ for _ in input_data.columns if _ not in sum_cols]

input_data = (
    input_data
    .group_by(groupby_cols)
    .agg(pl.sum(sum_cols))
 )

# Save input data
input_data.write_parquet(
    snakemake.output[0],
    compression='gzip'
    )
