rule network_objects:
    input:
        'results/network_analysis/{agg_lvl}/input/input_data.parquet.gzip'
    output:
        mirror_flows = 'results/network_analysis/{agg_lvl}/intermediary/mirror_flows.csv',
        edge_lists   = 'results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl'
    log:
        'workflow/logs/network_objects_{agg_lvl}.log'
    benchmark:
        'workflow/benchmarks/network_objects_{agg_lvl}.tsv'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/network_objects.py'