rule network_connectivity:
    input:
        'results/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/output/network_connectivity.csv'
    log:
        'workflow/logs/network_connectivity.log'
    threads: 2
    conda:
        '../envs/network_connectivity.yaml'
    script: 
        '../scripts/network_connectivity.py'