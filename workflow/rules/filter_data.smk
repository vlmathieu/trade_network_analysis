rule filter_data:
    input:
        'results/global/merged_data.parquet.gzip'
    output:
        'results/network_analysis/input/input_data.parquet.gzip'
    params:
        year_start      = config['years']['start'],
        year_stop       = config['years']['stop'],
        excluded_iso    = config['excluded_iso'],
        col_keep        = config['col_keep']
    threads: 2
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/filter_data.py'