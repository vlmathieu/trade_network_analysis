rule network_composition:
    input:
        'results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl'
    output:
        nominal  = 'results/network_analysis/{agg_lvl}/output/network_composition.csv',
        deflated = 'results/network_analysis/{agg_lvl}/output/network_composition_deflated.csv'
    params:
        value_col_deflated = config['deflation']['weight']
    log:
        'workflow/logs/network_composition_{agg_lvl}.log'
    benchmark:
        'workflow/benchmarks/network_composition_{agg_lvl}.tsv'
    threads: 1
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_composition.py'