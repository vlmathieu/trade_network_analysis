rule contributor_profiles:
    input:
        'results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/{agg_lvl}/output/contributor_profiles.csv'
    log:
        "workflow/logs/contributor_profiles_{agg_lvl}.log"
    benchmark:
        'workflow/benchmarks/contributor_profiles_{agg_lvl}.tsv'
    threads: 1
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/contributor_profiles.py'