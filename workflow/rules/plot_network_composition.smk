rule plot_network_composition:
    input:
        composition = expand('results/network_analysis/{agg_lvl}/output/network_composition.csv',
                agg_lvl  = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/network_composition.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = ['png', 'svg']),
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/network_mirrored_desc_stat_{metric}.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                metric  = ['value', 'weight'],
                ext     = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions_agg'],
        ext             = ['png', 'svg']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_composition.R'