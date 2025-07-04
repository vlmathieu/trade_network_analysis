rule plot_network_connectivity:
    input:
        'results/network_analysis/output/network_connectivity.csv'
    output:
        expand('results/network_analysis/plot/{fao_div}/network_connectivity.{ext}',
               fao_div = config['fao_divisions'],
               ext = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions'],
        ext             = ['png', 'svg']
    log:
        expand('workflow/logs/network_connectivity/{fao_div}_{ext}.log',
               fao_div = config['fao_divisions'],
               ext = ['png', 'svg'])
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_connectivity.R'