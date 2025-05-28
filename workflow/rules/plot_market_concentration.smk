rule plot_market_concentration:
    input:
        'results/network_analysis/output/market_concentration.csv'
    output:
        expand('results/network_analysis/plot/{fao_div}/market_concentration.{ext}',
               fao_div = config['fao_divisions'],
               ext = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions'],
        ext             = ['png', 'svg']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_market_concentration.R'