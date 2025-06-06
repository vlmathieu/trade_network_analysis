rule plot_contributor_profiles:
    input:
        'results/network_analysis/output/contributor_profiles.csv'
    output:
        expand('results/network_analysis/plot/{fao_div}/contributor_profiles/profiles_{year}.{ext}',
               fao_div = config['fao_divisions'],
               year = list(range(config['years']['start'], 
                                 config['years']['stop']-1)),
               ext = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions'],
        years           = list(range(config['years']['start'], 
                                     config['years']['stop']-1)),
        ext             = ['png', 'svg'],
        size            = "primary_value"
    log:
        expand('workflow/logs/plot_contributor_profiles/{fao_div}_{year}_{ext}.log',
               fao_div = config['fao_divisions'],
               year = list(range(config['years']['start'], 
                                     config['years']['stop']-1)),
               ext = ['png', 'svg'])
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_contributor_profiles.R'