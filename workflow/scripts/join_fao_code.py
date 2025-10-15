from snakemake.script import snakemake
import logging
import polars as pl
import polars.selectors as cs

def join_uncomtrade_FAO(FAO_HS: pl.dataframe.frame.DataFrame, 
                        uncomtrade_data: pl.dataframe.frame.DataFrame):
    '''
    Function that joins the corresponding FAO forest product codes to the HS 
    codes in UNComtrade data.

    Parameters
    ----------
    FAO_HS : polars dataframe
        The corresponding table between HS codes and FAO forest products 
        classification.
    uncomtrade_data : polars dataframe
        The UNComtrade data (commodities in HS codes).

    Returns
    -------
    uncomtrade_data_join : polars dataframe
        The UN Comtrade data with joined FAO forest products codes corresponding
        to HS codes.
    '''

    # Collect HS classification version code in uncomtrade data (H0, H1, ...)
    HS_vers_code = (sorted(
        uncomtrade_data.select('classificationCode')
                        .unique()
                        .to_series()
                        .to_list()
        )
    )

    # Collect HS classification version year in FAO_HS corresponding table
    # N.B. 'HS 1996' is duplicated to match H0 and H1 in HS_vers_code
    HS_vers_year = (
        [code for code in FAO_HS.columns if '1996' in code] +
        [code for code in FAO_HS.columns if 'HS' in code]
    )

    # Collect FAO columns to join to uncomtrade data
    FAO_col = [col for col in FAO_HS.columns if 'HS' not in col]

    # Join FAO codes to uncomtrade data by batch of HS classification version
    uncomtrade_joined_batch = (
        [
            (uncomtrade_data
             .filter(pl.col('classificationCode') == HS_code)
             .join(FAO_HS.unique(subset=FAO_col+[HS_year]),
                   left_on='cmdCode',
                   right_on=HS_year,
                   how='left')
                   .drop(cs.starts_with('HS'))
            )
            for HS_code, HS_year in zip(HS_vers_code,HS_vers_year)
        ]
    )

    # Concatenate all batch
    uncomtrade_data_join = pl.concat(
        [df for df in uncomtrade_joined_batch if df.shape != (0,0)],
        how='vertical_relaxed'
    )

    return uncomtrade_data_join

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load correspondance between HS codes and FAO forest product classification
FAO_HS = pl.read_json(snakemake.input[0])

# Load UNComtrade data with new deflated trade value column
deflate_uncomtrade = pl.read_parquet(snakemake.input[1])

# Join UNComtrade data with FAO_HS corresponding table
merged_data = join_uncomtrade_FAO(FAO_HS, deflate_uncomtrade)
logging.info(f"\nMerged dataframe head: \n {merged_data.head(5)}\n")
logging.info(f"\nMerged dataframe size (rows, columns): {merged_data.shape}\n")

# Save merged data
merged_data.write_parquet(
    snakemake.output[0],
    compression='gzip'
    )
