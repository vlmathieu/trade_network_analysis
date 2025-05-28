from snakemake.script import snakemake
import polars as pl

def reshape_uvi_data(uvi_data: pl.dataframe.frame.DataFrame):
    '''
    Functions that reshapes World Bank unit value index (UVI) data for deflation
    operation on uncomtrade data.

    Parameters
    ----------
    uvi_data : polars dataframe
        World Bank data on import/export unit values indexes.

    Returns
    -------
    uvi_data_reshaped : polars dataframe
        Reshaped World Bank data on import/export unit values indexes.

    '''
    
    # Convert time to values in years
    uvi_data = uvi_data.with_columns(pl.col('time').str.replace(r'YR', ''))

    # Split data between export and import
    uvi_data_exp = (uvi_data.select(['economy','time','TX.UVI.MRCH.XD.WD'])
                            .rename({'TX.UVI.MRCH.XD.WD': 'index'})
                            .with_columns(
                                flowDesc = pl.lit('Export')
                            ))
    
    uvi_data_imp = (uvi_data.select(['economy','time','TM.UVI.MRCH.XD.WD'])
                            .rename({'TM.UVI.MRCH.XD.WD': 'index'})
                            .with_columns(
                                flowDesc = pl.lit('Import')
                            ))
    
    # Merge uvi data on export and on import
    uvi_data_reshaped = (pl.concat([uvi_data_exp,uvi_data_imp], how="vertical_relaxed")
                           .rename({'economy':'reporterISO','time':'period'})
                        )
    
    return uvi_data_reshaped

def deflate(uncomtrade_data: pl.dataframe.frame.DataFrame, 
            uvi_data: pl.dataframe.frame.DataFrame):
    '''
    Function that deflate trade values using World Bank import/export unit 
    values indexes.

    Parameters
    ----------
    comtrade_data : polars dataframe
        UNComtrade trade data on which deflation operation is conducted.
    uvi_data : dataframe
        Dataframe containing unit values indexes through time.

    Returns
    -------
    data_deflated : polars dataframe
        Dataframe of the input trade data with deflated primary value as an 
        additional column (nan in case of missing data).

    '''
    
    uvi_data_reshaped = reshape_uvi_data(uvi_data)
    
    data_deflated = uncomtrade_data.join(uvi_data_reshaped,
                                         on=['reporterISO','period','flowDesc'],
                                         how='left')
    
    data_deflated = (
        data_deflated.with_columns(
            primaryValue_deflated = (pl.col('primaryValue')/pl.col('index')*100)
            ).drop('index')
    )
        
    return data_deflated

# Load data and select unit value index data
uvi_data = (
    pl.read_csv(snakemake.input[0], separator=';')
    .select(['economy', 'time', 'TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'])
)

# Load trade data
uncomtrade_data = pl.read_parquet(snakemake.input[1])

# Deflate primary column in trade date in current USD into constant USD
deflate_uncomtrade = deflate(uncomtrade_data, uvi_data)

# Save processed data
deflate_uncomtrade.write_parquet(
        snakemake.output[0],
        compression='gzip'
        )
