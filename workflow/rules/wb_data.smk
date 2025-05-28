rule wb_data:
    output:
        'resources/public/wb_series_data.csv',
    params:
        year_start  = config['years']['start'],
        year_stop   = config['years']['stop'],
        wb_series   = config['wb_series']
    threads: 1
    conda:
        '../envs/wbgapi.yaml'
    script:
        '../scripts/wb_data.py'