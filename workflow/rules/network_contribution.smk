rule network_contribution:
    input:
        'results/processed_data/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/processed_data/network_analysis/output/network_contribution.csv'
    params:
        weight      = config['weight']
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_contribution.py'