rule network_connectivity:
    input:
        expand('results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl',
                agg_lvl  = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/output/network_connectivity.csv',
                agg_lvl  = config['agg_lvl'])
    log:
        'workflow/logs/network_connectivity.log'
    threads: 2
    conda:
        '../envs/network_connectivity.yaml'
    script: 
        '../scripts/network_connectivity.py'