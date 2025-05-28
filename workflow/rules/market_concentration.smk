rule market_concentration:
    input:
        'results/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/output/market_concentration.csv'
    params:
        weight      = config['weight']
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/market_concentration.py'