rule plot_network_contribution:
    input:
        'results/network_analysis/output/network_contribution.csv',
        'results/network_analysis/output/contributor_profiles.csv'
    output:
        expand('results/network_analysis/plot/{fao_div}/network_contribution.{ext}',
               fao_div = config['fao_divisions'],
               ext = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions'],
        ext             = ['png', 'svg']
    log:
        expand('workflow/logsnetwork_contribution/{fao_div}_{ext}.log',
               fao_div = config['fao_divisions'],
               ext = ['png', 'svg'])
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_contribution.R'