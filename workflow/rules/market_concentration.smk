rule market_concentration:
    input:
        'results/processed_data/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/processed_data/network_analysis/output/market_concentration.csv'
    params:
        weight      = config['weight']
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/market_concentration.py'