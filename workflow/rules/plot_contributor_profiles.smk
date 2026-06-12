rule plot_contributor_profiles:
    input:
        'results/network_analysis/agg_eu/output/contributor_profiles.csv'
    output:
        expand('results/network_analysis/agg_eu/plot/01/contributor_profiles.{ext}',
               ext = ['png', 'svg'])
    params:
        fao_divisions = ['01'],
        year_start    = 2000,
        year_end      = 2020,
        ext           = ['png', 'svg'],
        size          = 'primary_value',
        threshold     = config['threshold_main_contributors']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_contributor_profiles.R'
