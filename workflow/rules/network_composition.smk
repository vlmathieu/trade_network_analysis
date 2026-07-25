rule network_composition:
    input:
        'results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/{agg_lvl}/output/network_composition.csv'
    log:
        'workflow/logs/network_composition_{agg_lvl}.log'
    benchmark:
        'workflow/benchmarks/network_composition_{agg_lvl}.tsv'
    threads: 1
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_composition.py'