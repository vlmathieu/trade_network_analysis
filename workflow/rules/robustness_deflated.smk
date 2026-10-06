# Robustness check on inflation-adjusted value (supplementary material).
# Each rule reuses the core plot script with a different weight and output
# stem; the core figures are untouched. Figures are paired nominal vs deflated
# over the full span, with the under-covered years (config deflation.flag_years)
# shaded on the deflated panel.

DEFL_LVL = config['deflation']['agg_lvl']
DEFL_DIV = config['deflation']['fao_div']
DEFL_OUT = f'results/network_analysis/{DEFL_LVL}/output'
DEFL_PLT = f'results/network_analysis/{DEFL_LVL}/plot/{DEFL_DIV}'


rule plot_market_concentration_deflated:
    input:
        f'{DEFL_OUT}/market_concentration.csv'
    output:
        expand(f'{DEFL_PLT}/market_concentration_deflated.{{ext}}',
               ext = config['figure_ext'])
    params:
        fao_divisions = [DEFL_DIV],
        wgt           = ['primary_value', config['deflation']['weight']],
        flag_years    = config['deflation']['flag_years'],
        out_stem      = 'market_concentration_deflated',
        ext           = config['figure_ext']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_market_concentration.R'


rule plot_network_contribution_deflated:
    input:
        f'{DEFL_OUT}/network_contribution.csv',
        f'{DEFL_OUT}/network_composition.csv'
    output:
        expand(f'{DEFL_PLT}/network_contribution_deflated.{{ext}}',
               ext = config['figure_ext'])
    params:
        fao_divisions = [DEFL_DIV],
        weights       = ['primary_value', config['deflation']['weight']],
        flag_years    = config['deflation']['flag_years'],
        out_stem      = 'network_contribution_deflated',
        top_frac      = config['network_contribution']['top_frac'],
        time_span     = config['network_contribution']['time_span'],
        ext           = config['figure_ext']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_network_contribution.R'


rule plot_network_composition_deflated:
    input:
        composition = f'{DEFL_OUT}/network_composition_deflated.csv'
    output:
        expand(f'{DEFL_PLT}/network_composition_deflated.{{ext}}',
               ext = config['figure_ext'])
    params:
        fao_divisions = [DEFL_DIV],
        flag_years    = config['deflation']['flag_years'],
        out_stem      = 'network_composition_deflated',
        mirror        = False,
        ext           = config['figure_ext']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_network_composition.R'


rule plot_contributor_profiles_deflated:
    input:
        f'{DEFL_OUT}/contributor_profiles.csv'
    output:
        expand(f'{DEFL_PLT}/contributor_profiles_deflated.{{ext}}',
               ext = config['figure_ext'])
    params:
        fao_divisions = [DEFL_DIV],
        year_start    = config['reference_years']['start'],
        year_end      = config['reference_years']['end'],
        size          = config['deflation']['weight'],
        threshold     = config['threshold_main_contributors'],
        out_stem      = 'contributor_profiles_deflated',
        size_breaks   = [1000, 3000, 6000],
        zoom_a_y      = [23, 90],
        ext           = config['figure_ext']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_contributor_profiles.R'


rule plot_chord_diagram_deflated:
    input:
        f'results/network_analysis/{DEFL_LVL}/intermediary/mirror_flows.csv',
        f'{DEFL_OUT}/contributor_profiles.csv'
    output:
        expand(f'{DEFL_PLT}/chord_diagram_deflated.{{ext}}',
               ext = config['figure_ext'])
    params:
        fao_divisions = [DEFL_DIV],
        year_start    = config['reference_years']['start'],
        year_end      = config['reference_years']['end'],
        chord_year    = config['chord']['year'],
        chord_n       = config['chord']['n_top'],
        threshold     = config['trade_network']['node_threshold'],
        panels        = 'deflated',
        out_stem      = 'chord_diagram_deflated',
        ext           = config['figure_ext']
    threads: 1
    conda:
        '../envs/r_trade_flows.yaml'
    script:
        '../scripts/plot_trade_flows.R'
