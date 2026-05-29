rule market_concentration:
    input:
        expand('results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl',
                agg_lvl  = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/output/market_concentration.csv',
                agg_lvl  = config['agg_lvl'])
    params:
        weight      = config['weight']
    log:
        'workflow/logs/market_concentration.log'
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/market_concentration.py'