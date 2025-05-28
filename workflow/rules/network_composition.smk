rule network_composition:
    input:
        'results/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/output/network_composition.csv'
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_composition.py'