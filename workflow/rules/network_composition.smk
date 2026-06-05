rule network_composition:
    input:
        expand('results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl',
                agg_lvl  = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/output/network_composition.csv',
                agg_lvl  = config['agg_lvl'])
    log:
        'workflow/logs/network_composition.log'
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_composition.py'