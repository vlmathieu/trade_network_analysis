rule contributor_profiles:
    input:
        expand('results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl',
                agg_lvl = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/output/contributor_profiles.csv',
                agg_lvl = config['agg_lvl'])
    log:
        "workflow/logs/contributor_profiles.log"
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/contributor_profiles.py'