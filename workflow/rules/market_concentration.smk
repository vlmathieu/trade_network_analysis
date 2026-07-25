rule market_concentration:
    input:
        'results/network_analysis/{agg_lvl}/intermediary/edge_lists.pkl'
    output:
        'results/network_analysis/{agg_lvl}/output/market_concentration.csv'
    params:
        weight = config['weight']
    log:
        'workflow/logs/market_concentration_{agg_lvl}.log'
    benchmark:
        'workflow/benchmarks/market_concentration_{agg_lvl}.tsv'
    threads: 1
    conda:
        '../envs/network_metrics.yaml'
    script:
        '../scripts/market_concentration.py'