import polars as pl
import polars.selectors as cs
import pickle
import re
import numpy as np
import networkx as nx
from scipy.stats import kurtosis
from scipy.stats import skew
import pandas as pd
import itertools

# Parameters
HS_version = [1996, 2002, 2007, 2012, 2017, 2022]
year_start = 2020
year_stop = 2024
excluded_iso = ['XX', '_X', '\\d']
flow_to_keep = ['M', 'X']
col_keep = ['period', 
           'reporterISO', 
           'reporterDesc', 
           'flowCode', 
           'partnerISO', 
           'partnerDesc',
           'FAO Code',
           'FAO Product',  
           'netWgt', 
           'primaryValue', 
           'primaryValue_deflated']

# Load data
FAO_HS = pl.read_json('/Users/valentinmathieu/Desktop/wd/trade_network_analysis/resources/inhouse/correspondence_FAO_HS.json')

comtrade_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_network_analysis/resources/public/uncomtrade_data.parquet.gzip')

wb_data = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_network_analysis/resources/public/wb_series_data.csv',
                      separator=';')

deflate_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/global/deflate_uncomtrade_data.parquet.gzip')

merged_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/global/merged_data.parquet.gzip')

input_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/input/input_data.parquet.gzip')

mirror_flows = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/intermediary/mirror_flows.csv',
    separator=';'
)

path = '/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/intermediary/edge_lists.pkl'
with open(path, 'rb') as f:
    net_dict = pickle.load(f)

contributor_profiles = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/output/contributor_profiles.csv',
    separator=';'
)

market_concentration = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/output/market_concentration.csv',
    separator=';'
)

sorted(market_concentration.select('period').unique().to_series().to_list())
market_concentration.filter(pl.col('cmd') == 12).sort('period').select('hhi_imp')

network_composition = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/output/network_composition.csv',
    separator=';'
)

network_connectivity = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/output/network_connectivity.csv',
    separator=';'
)

network_contribution = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/output/network_contribution.csv',
    separator=';'
)

# Draft

list(net_dict.keys())
dict(net_dict[(12, 2022)])
edge_lists = [{k: v} for (k, v) in net_dict.items()]
unit_edge_list_dict = edge_lists[-1]
[[keys, edge_list]] = unit_edge_list_dict.items()
type(net_dict[(12, 2022)])

keys = (12, 1996)
edge_list = net_dict[(12, 1996)]
cmd, period = keys
to_omit = "Japan"
weight = 'primary_value'

######
network_contribution.filter(pl.col('country') == 'Japan')
mirror_flows.columns
tot_val_exp = (mirror_flows.filter(pl.col('period') == 1996,
                                   pl.col('fao_code') == 12)
                           .select(pl.sum('primary_value_exp'))
               )
tot_val_imp = (mirror_flows.filter(pl.col('period') == 1996,
                                   pl.col('fao_code') == 12)
                           .select(pl.sum('primary_value_imp'))
               )
tot_val_exp_wo_Japan = (mirror_flows.filter(pl.col('period') == 1996,
                                   pl.col('fao_code') == 12,
                                   ((pl.col('exporter_desc') != 'Japan') &
                                   (pl.col('importer_desc') != 'Japan')))
                           .select(pl.sum('primary_value_exp'))
               )
tot_val_imp_wo_Japan = (mirror_flows.filter(pl.col('period') == 1996,
                                   pl.col('fao_code') == 12,
                                   ((pl.col('exporter_desc') != 'Japan') &
                                   (pl.col('importer_desc') != 'Japan')))
                           .select(pl.sum('primary_value_imp'))
               )
(tot_val_imp - tot_val_imp_wo_Japan) / tot_val_imp
(tot_val_exp - tot_val_exp_wo_Japan) / tot_val_exp

edge_list = net_dict[(12, 1996)]
# Build directed network based on edge_list
net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)

# Replace None weights by 0
for _,_,d in net.edges(data=True):
    for key in d:
        if d[key] is None:
            d[key] = 0

# Remove country on which contribution is calculated from network
to_omit = 'Japan'
net_omit = nx.from_edgelist(edge_list, create_using=nx.DiGraph)
net_omit.remove_node(to_omit)

# Replace None weights by 0
for _,_,d in net_omit.edges(data=True):
    for key in d:
        if d[key] is None:
            d[key] = 0

# Build unweighted degree lists for exp=out_degree | imp=in_degree
# Without omission
degree_unweighted = [
    net.out_degree(x) for x in net.nodes()
]
# With omission
degree_unweighted_omit = [
    net_omit.out_degree(x) for x in net_omit.nodes()
]

# Build weighted degree lists for exporters=out_degree | importers=in_degree
# Without omission
degree_weighted_exp = [
    net.out_degree(x, weight = f'{weight}_exp') for x in net.nodes()
]
degree_weighted_imp = [
    net.in_degree(x, weight = f'{weight}_imp') for x in net.nodes()
]
# With omission
degree_weighted_exp_omit = [
    net_omit.out_degree(x, weight = f'{weight}_exp') 
    for x in net_omit.nodes()
]
degree_weighted_imp_omit = [
    net_omit.in_degree(x, weight = f'{weight}_imp') 
    for x in net_omit.nodes()
]

(sum(degree_weighted_exp) - sum(degree_weighted_exp_omit)) / sum(degree_weighted_exp)
(sum(degree_weighted_imp) - sum(degree_weighted_imp_omit)) / sum(degree_weighted_imp)

unit_network_contribution = pl.from_dict(
        {
            "period": 1996,
            "cmd": 12,
            # Country on which contribution is calculated
            "country": to_omit,
            # Assign total number of edges in network
            # Without omission
            "nb_edges": sum(degree_unweighted),
            # With omission
            "nb_edges_country": (sum(degree_unweighted) 
                                    - sum(degree_unweighted_omit)),
            # Assign total circulating value for exports and imports
            # Without omission
            "traded_value_exp": sum(degree_weighted_exp),
            "traded_value_imp": sum(degree_weighted_imp),
            # With omission
            "traded_value_exp_country": (sum(degree_weighted_exp) -
                                            sum(degree_weighted_exp_omit)),
            "traded_value_imp_country": (sum(degree_weighted_imp) -
                                            sum(degree_weighted_imp_omit))
        }
)

# Compute country contribution to total traded value and nb. of edges
unit_network_contribution = unit_network_contribution.with_columns(
    (pl.col("nb_edges_country") / 
    pl.col("nb_edges")).alias("contrib_nb_edges"),
    (pl.col("traded_value_exp_country") / 
    pl.col("traded_value_exp")).alias("contrib_trade_value_exp"),
    (pl.col("traded_value_imp_country") / 
    pl.col("traded_value_imp")).alias("contrib_trade_value_imp")
)

network_contribution.filter(pl.col('country') == 'Japan',
                            pl.col('cmd') == 12,
                            pl.col('period') == 1996)