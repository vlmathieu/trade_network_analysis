rule correspondence_FAO_HS:
    output:
        'resources/inhouse/correspondence_FAO_HS.json'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/correspondence_FAO_HS.py'