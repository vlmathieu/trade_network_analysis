from snakemake.script import snakemake
import polars as pl
import matplotlib.pyplot as plt

data = pl.read_parquet(snakemake.input['data'])

fig, ax = plt.subplots()

ax.scatter(x=data['netWgt'], y=data['primaryValue'])

ax.set(xlim=(0, 8),
       ylim=(0, 300))

plt.savefig(snakemake.output['plot'])
