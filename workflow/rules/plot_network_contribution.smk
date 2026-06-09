rule plot_network_contribution:
    input:
        expand('results/network_analysis/{agg_lvl}/output/network_contribution.csv',
                agg_lvl = config['agg_lvl']),
        expand('results/network_analysis/{agg_lvl}/output/network_composition.csv',
                agg_lvl = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/network_contribution.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = ['png', 'svg'])
    params:
        fao_divisions = config['fao_divisions_agg'],
        top_frac      = 0.03,
        time_span     = 5,
        ext           = ['png', 'svg']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_network_contribution.R'
