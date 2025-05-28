rule network_objects:
    input:
        'results/processed_data/network_analysis/input/input_data.parquet.gzip'
    output:
        'results/processed_data/network_analysis/intermediary/mirror_flows.csv',
        'results/processed_data/network_analysis/intermediary/edge_lists.pkl'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/network_objects.py'