rule merge_data:
    input:
        'resources/inhouse/correspondence_FAO_HS.json',
        'results/global/deflate_uncomtrade_data.parquet.gzip'
    output:
        'results/global/merged_data.parquet.gzip'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/merge_data.py'