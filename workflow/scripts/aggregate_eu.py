from snakemake.script import snakemake
import logging
import polars as pl
import polars.selectors as cs

def agg_eu(trade_data: pl.dataframe.frame.DataFrame, 
           eu_iso: list):
    '''
    Function that aggregates trade between EU and each trading partner and 
    deletes intra-EU trade.

    Parameters
    ----------
    trade_data : polars dataframe
        The trade data on which aggregation must be undergone.
    eu_iso : list
        The list of all ISO 3 codes for EU member countries.

    Returns
    -------
    trade_data_agg_eu : polars dataframe
        The trade data with aggregated trade flows between EU and each trading 
        partner and without intra-EU trade flows.
    '''

    # Convert each European Union country desc and iso
    # desc becomes "European Union"
    # iso becomes "EU"
    trade_data_eu = (
        trade_data
        .with_columns(
            reporter_desc = pl.when(pl.col('reporter_iso').is_in(eu_iso))
                            .then(pl.lit('European Union'))
                            .otherwise(pl.col('reporter_desc')),
            partner_desc = pl.when(pl.col('partner_iso').is_in(eu_iso))
                            .then(pl.lit('European Union'))
                            .otherwise(pl.col('partner_desc')),
            reporter_iso = pl.when(pl.col('reporter_iso').is_in(eu_iso))
                            .then(pl.lit('EU'))
                            .otherwise(pl.col('reporter_iso')),
            partner_iso = pl.when(pl.col('partner_iso').is_in(eu_iso))
                            .then(pl.lit('EU'))
                            .otherwise(pl.col('partner_iso'))
        )
        .filter(pl.col('reporter_desc') != pl.col('partner_desc'))
    )

    # Sum netwgt and primary value by all other columns, i.e., 
    # {'period', 'reporter_iso', 'reporter_desc', 'flow_code', 'partner_iso', 'partner_desc', 'fao_code_agg', 'fao_product_agg'}
    # ~ aggregation for duplicated trade flows EU <> partner | reporter <> EU
    sum_cols = (
        trade_data_eu
        .select(cs.starts_with('net_wgt') | cs.starts_with('primary_value'))
        .columns
    )
    groupby_cols = [_ for _ in trade_data_eu.columns if _ not in sum_cols]

    trade_data_agg_eu = (
        trade_data_eu
        .group_by(groupby_cols)
        .agg(
            pl.when(pl.col(sum_cols).count() > 0).then(pl.sum(sum_cols))
        )
    )

    return trade_data_agg_eu

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load UNComtrade data with new deflated trade value column and merged data
input_data = pl.read_parquet(snakemake.input[0])

# Collect EU countries desc
eu_iso = snakemake.params['eu_iso']

# Aggregate trade between EU and each trading partner and delete intra-EU trade
input_data_agg_eu = agg_eu(trade_data=input_data, eu_iso=eu_iso)
logging.info(f"\nColumns:\n {input_data_agg_eu.columns}\n")
logging.info(f"\nDataframe overview:\n {input_data_agg_eu}\n")
logging.info(f"\nAggregate data descriptive stats:\n {input_data_agg_eu.describe()}\n")


# Save aggregated data
input_data_agg_eu.write_parquet(
        snakemake.output[0],
        compression='gzip'
        )