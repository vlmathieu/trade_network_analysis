rule merge_data:
    input:
        'resources/raw_data/inhouse/correspondence_FAO_HS.json',
        'resources/raw_data/public/wb_series_data.csv',
        'resources/raw_data/public/wb_countries_data.csv',
        'results/processed_data/global/deflate_uncomtrade_data.parquet.gzip'
    output:
        'results/processed_data/global/merged_data.parquet.gzip'
    params:
        wb_series_drop      = ['TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'],
        wb_countries_keep   = ['id', 'longitude', 'latitude','capitalCity'],
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/merge_data.py'