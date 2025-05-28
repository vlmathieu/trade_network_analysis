rule plot_network_connectivity:
    input:
        'results/processed_data/network_analysis/output/network_connectivity.csv'
    output:
        expand('results/processed_data/network_analysis/plot/{fao_div}/network_connectivity.{ext}',
               fao_div = config['fao_divisions'],
               ext = ['png', 'svg'])
    params:
        fao_divisions   = config['fao_divisions'],
        ext             = ['png', 'svg']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_connectivity.R'