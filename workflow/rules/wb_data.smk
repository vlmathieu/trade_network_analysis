rule wb_data:
    output:
        'resources/raw_data/public/wb_series_data.csv',
        'resources/raw_data/public/wb_countries_data.csv'
    params:
        year_start  = config['years']['start'],
        year_stop   = config['years']['stop'],
        wb_series   = config['wb_series']
    threads: 1
    conda:
        '../envs/wbgapi.yaml'
    script:
        '../scripts/wb_data.py'