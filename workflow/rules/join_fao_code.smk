rule join_fao_code:
    input:
        'resources/inhouse/correspondence_FAO_HS.json',
        'results/global/deflate_uncomtrade_data.parquet.gzip'
    output:
        'results/global/merged_data.parquet.gzip'
    log:
        'workflow/logs/join_fao_code.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/join_fao_code.py'