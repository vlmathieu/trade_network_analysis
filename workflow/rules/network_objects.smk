rule network_objects:
    input:
        'results/network_analysis/input/input_data.parquet.gzip'
    output:
        'results/network_analysis/intermediary/mirror_flows.csv',
        'results/network_analysis/intermediary/edge_lists.pkl'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/network_objects.py'