rule plot_network_connectivity:
    input:
        expand('results/network_analysis/{agg_lvl}/output/network_connectivity.csv',
                agg_lvl = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/network_connectivity.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions_agg'],
        ext             = ['png', 'svg']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_connectivity.R'