rule network_connectivity:
    input:
        'results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/{agg_lvl}/output/network_connectivity.csv'
    log:
        'workflow/logs/network_connectivity_{agg_lvl}.log'
    benchmark:
        'workflow/benchmarks/network_connectivity_{agg_lvl}.tsv'
    threads: 1
    conda:
        '../envs/network_connectivity.yaml'
    script: 
        '../scripts/network_connectivity.py'