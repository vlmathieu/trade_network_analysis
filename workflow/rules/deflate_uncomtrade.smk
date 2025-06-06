rule deflate_uncomtrade:
    input:
        'resources/public/wb_series_data.csv',
        'resources/public/uncomtrade_data.parquet.gzip'
    output:
        'results/global/deflate_uncomtrade_data.parquet.gzip'
    log:
        'workflow/logs/deflate_uncomtrade.log'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/deflate_uncomtrade.py'