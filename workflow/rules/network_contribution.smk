rule network_contribution:
    input:
        # 'results/network_analysis/intermediary/edge_lists.pkl'
        expand('results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl',
                agg_lvl  = config['agg_lvl'])
    output:
        # 'results/network_analysis/output/network_contribution.csv'
        expand('results/network_analysis/{agg_lvl}/output/network_contribution.csv',
                agg_lvl  = config['agg_lvl'])
    params:
        weight      = config['weight']
    log:
        'workflow/logs/network_contribution.log'
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_contribution.py'