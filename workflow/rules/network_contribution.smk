rule network_contribution:
    input:
        'results/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/output/network_contribution.csv'
    params:
        weight      = config['weight']
    log:
        'workflow/logs/network_contribution.log'
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_contribution.py'