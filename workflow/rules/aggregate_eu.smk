rule aggregate_eu:
    input:
        'results/network_analysis/country_lvl/input/input_data.parquet.gzip'
    output:
        'results/network_analysis/agg_eu/input/input_data.parquet.gzip'
    log:
        'workflow/logs/aggregate_eu.log'
    params:
        eu_iso  = config['eu']['iso']
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/aggregate_eu.py'