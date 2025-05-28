rule network_connectivity:
    input:
        'results/processed_data/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/processed_data/network_analysis/output/network_connectivity.csv'
    threads: 2
    conda:
        '../envs/network_connectivity.yaml'
    script: 
        '../scripts/network_connectivity.py'