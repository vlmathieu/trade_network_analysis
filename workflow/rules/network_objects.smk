rule network_objects:
    input:
        expand('results/network_analysis/{agg_lvl}/input/input_data.parquet.gzip',
                agg_lvl  = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/intermediary/mirror_flows.csv',
                agg_lvl  = config['agg_lvl']),
        expand('results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl',
                agg_lvl  = config['agg_lvl'])
    log:
        'workflow/logs/network_objects.log'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/network_objects.py'