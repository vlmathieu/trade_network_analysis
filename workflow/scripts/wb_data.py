from snakemake.script import snakemake
import polars as pl
import wbgapi as wb

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

wb_countries_data = (
    pl.DataFrame(
        wb.economy.DataFrame(skipAggs=True).reset_index()
    )
)

check = (
    sorted(snakemake.params['wb_series']) == 
    sorted([_ for _ in wb_series_data.columns if _ not in {'economy', 'time'}])
)

if check:
    print('Data have been checked.\n')
    wb_series_data.write_csv(snakemake.output[0], separator=';')
    wb_countries_data.write_csv(snakemake.output[1], separator=';')
else:
    print('Issues found in data download.\n')
    print(sorted(snakemake.params['wb_series']))
    print(sorted([_ for _ in wb_series_data.columns 
                  if _ not in {'economy', 'time'}]))
