rule get_uncomtrade:
    input:
        'resources/inhouse/correspondence_FAO_HS.json'
    output:
        'resources/public/uncomtrade_data.parquet.gzip'
    params:
        year_start  = config['years']['start'],
        year_stop   = config['years']['stop'],
        HS_version  = config['HS_version'],
        flowCode    = list(str(flow) for flow in config['flowCode']),
        apikey      = os.environ['comtrade_apikey']
    threads: 1
    conda:
        '../envs/comtradeapicall.yaml'
    script:
        '../scripts/get_uncomtrade_data.py'