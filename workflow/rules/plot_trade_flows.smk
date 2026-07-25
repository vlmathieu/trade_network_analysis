rule plot_trade_flows:
    input:
        expand('results/network_analysis/{agg_lvl}/intermediary/mirror_flows.csv',
                agg_lvl = config['agg_lvl']),
        expand('results/network_analysis/{agg_lvl}/output/contributor_profiles.csv',
                agg_lvl = config['agg_lvl'])
    output:
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/chord_diagram.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = config['figure_ext']),
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/chord_diagram_fob.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = config['figure_ext']),
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/trade_network.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = config['figure_ext'])
    params:
        fao_divisions   = config['fao_divisions_agg'],
        year_start      = config['reference_years']['start'],
        year_end        = config['reference_years']['end'],
        chord_year      = config['chord']['year'],
        threshold       = config['trade_network']['node_threshold'],
        chord_n         = config['chord']['n_top'],
        ext             = config['figure_ext']
    threads: 1
    conda:
        '../envs/r_trade_flows.yaml'
    script:
        '../scripts/plot_trade_flows.R'
