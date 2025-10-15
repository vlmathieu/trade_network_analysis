from snakemake.script import snakemake
import logging
import polars as pl
import wbgapi as wb

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

wb_list = (
    [pl.DataFrame(
        wb.data.DataFrame(serie, 
                          time=range(
                              snakemake.params['year_start'],
                              snakemake.params['year_stop']), 
                          columns='series')
                          .reset_index()
    ) for serie in snakemake.params['wb_series']]
)

wb_series_data = (
    pl.concat(wb_list, how='align')
)
logging.info(f"\nWB series:\n {wb_series_data}\n")

check = (
    sorted(snakemake.params['wb_series']) == 
    sorted([_ for _ in wb_series_data.columns if _ not in {'economy', 'time'}])
)

if check:
    logging.info('Data have been checked.\n')
    wb_series_data.write_csv(snakemake.output[0], separator=';')
else:
    logging.info('Issues found in data download.\n')
    logging.info(sorted(snakemake.params['wb_series']))
    logging.info(sorted([_ for _ in wb_series_data.columns 
                         if _ not in {'economy', 'time'}]))
