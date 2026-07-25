rule network_contribution:
    input:
        'results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/{agg_lvl}/output/network_contribution.csv'
    params:
        weight = config['weight']
    log:
        'workflow/logs/network_contribution_{agg_lvl}.log'
    benchmark:
        'workflow/benchmarks/network_contribution_{agg_lvl}.tsv'
    threads: 1
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_contribution.py'