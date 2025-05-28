rule format_data_network_analysis:
    input:
        'results/processed_data/global/merged_data.parquet.gzip'
    output:
        'results/processed_data/network_analysis/input/input_data.parquet.gzip'
    params:
        year_start      = config['years']['start'],
        year_stop       = config['years']['stop'],
        flow_to_keep    = config['flow_to_keep'],
        fao_divisions   = config['fao_divisions'],
        excluded_iso    = config['excluded_iso'],
        col_keep        = config['col_keep']
    threads: 2
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/format_data_network_analysis.py'