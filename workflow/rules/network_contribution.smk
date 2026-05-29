rule network_contribution:
    input:
        expand('results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl',
                agg_lvl  = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/output/network_contribution.csv',
                agg_lvl  = config['agg_lvl'])
    params:
        weight      = config['weight']
    log:
        'workflow/logs/network_contribution.log'
    threads: 4
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_contribution.py'