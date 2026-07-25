rule plot_prices_figures:
    input:
        composition = expand('results/network_analysis/{agg_lvl}/output/network_composition.csv',
                agg_lvl  = config['agg_lvl']),
        mirror_flows = 'results/network_analysis/country_lvl/intermediary/mirror_flows.csv'
    output:
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/network_price_desc_stat.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = config['figure_ext']),
        expand('results/network_analysis/country_lvl/plot/01/network_price_{fig}.{ext}',
                fig     = ['partners',
                           'ridgeline_countries_imp',
                           'ridgeline_countries_exp',
                           'ridgeline_years_imp',
                           'ridgeline_years_exp'],
                ext     = config['figure_ext'])
    params:
        fao_divisions   = config['fao_divisions_agg'],
        ext             = config['figure_ext']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_prices_figures.R'
