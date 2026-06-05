rule plot_market_concentration:
    input:
        expand('results/network_analysis/{agg_lvl}/output/market_concentration.csv',
                agg_lvl  = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/market_concentration.{ext}',
                agg_lvl  = config['agg_lvl'],
                fao_div  = config['fao_divisions_agg'],
                ext      = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions_agg'],
        wgt             = config['weight'],
        ext             = ['png', 'svg']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_market_concentration.R'
