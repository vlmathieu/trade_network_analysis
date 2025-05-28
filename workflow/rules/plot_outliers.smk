rule plot_outliers:
    input:
        data='resources/raw_data/uncomtrade_4403.parquet.gzip'
    output:
        plot='results/plot/outliers/outliers_4403.png'
    script:
        '../scripts/plot_outliers.py'
