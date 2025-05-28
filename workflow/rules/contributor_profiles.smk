rule contributor_profiles:
    input:
        'results/processed_data/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/processed_data/network_analysis/output/contributor_profiles.csv'
    params:
        threshold   = config['threshold_main_contributors'],
        weight      = config['weight']
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/contributor_profiles.py'