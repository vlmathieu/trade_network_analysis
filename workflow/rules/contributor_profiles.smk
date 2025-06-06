rule contributor_profiles:
    input:
        'results/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/output/contributor_profiles.csv'
    params:
        threshold   = config['threshold_main_contributors'],
        weight      = config['weight']
    log:
        "workflow/logs/contributor_profiles.log"
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/contributor_profiles.py'