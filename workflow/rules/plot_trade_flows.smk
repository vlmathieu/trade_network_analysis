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
                ext     = ['png', 'svg']),
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/chord_diagram_fob.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = ['png', 'svg']),
        expand('results/network_analysis/{agg_lvl}/plot/{fao_div}/trade_network.{ext}',
                agg_lvl = config['agg_lvl'],
                fao_div = config['fao_divisions_agg'],
                ext     = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions_agg'],
        year_start      = 2000,
        year_end        = 2020,
        chord_year      = 2020,
        threshold       = config['threshold_main_contributors'],
        chord_n         = 10,
        ext             = ['png', 'svg']
    threads: 1
    conda:
        '../envs/r_trade_flows.yaml'
    script:
        '../scripts/plot_trade_flows.R'
