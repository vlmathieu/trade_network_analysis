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
FAO_HS = pl.read_json('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/inhouse/correspondence_FAO_HS.json')

comtrade_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/public/uncomtrade_data.parquet.gzip')

wb_countries = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/public/wb_countries_data.csv',
                      separator=';')
wb_data = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/public/wb_series_data.csv',
                      separator=';')

deflate_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/deflate_uncomtrade_data.parquet.gzip')

merged_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/merged_data.parquet.gzip')

input_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/input/input_data.parquet.gzip')

mirror_flows = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/intermediary/mirror_flows.csv',
    separator=';'
)

path = '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/intermediary/edge_lists.pkl'
with open(path, 'rb') as f:
    net_dict = pickle.load(f)

contributor_profiles = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/contributor_profiles.csv',
    separator=';'
)

market_concentration = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/market_concentration.csv',
    separator=';'
)

sorted(market_concentration.select('period').unique().to_series().to_list())
market_concentration.filter(pl.col('cmd') == 12).sort('period').select('hhi_imp')

network_composition = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/network_composition.csv',
    separator=';'
)

network_connectivity = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/network_connectivity.csv',
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

cmd, period = keys

# Build directed network based on edge_list
net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)

# Replace None weights by 0
for _,_,d in net.edges(data=True):
    for key in d:
        if d[key] is None:
            d[key] = 0

net.degree('Japan', weight = f'{weight}_imp')
net.degree('Japan')
# Remove country on which contribution is calculated from network
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

(sum(degree_weighted_exp) - sum(degree_weighted_exp_omit)) / sum(degree_weighted_exp) * 100
(sum(degree_weighted_imp) - sum(degree_weighted_imp_omit)) / sum(degree_weighted_imp) * 100

# Build directed network based on edge_list
net = nx.from_edgelist(edge_list, create_using=nx.Graph)

# Replace None weights by 0
for _,_,d in net.edges(data=True):
    for key in d:
        if d[key] is None:
            d[key] = 0

# Remove country on which contribution is calculated from network
net_omit = nx.from_edgelist(edge_list, create_using=nx.Graph)
net_omit.remove_node(to_omit)

# Replace None weights by 0
for _,_,d in net_omit.edges(data=True):
    for key in d:
        if d[key] is None:
            d[key] = 0

# Build weighted degree lists for exporters=out_degree | importers=in_degree
# Without omission
degree_weighted_exp = [
    net.degree(x, weight = f'{weight}_exp') for x in net.nodes()
]
degree_weighted_imp = [
    net.degree(x, weight = f'{weight}_imp') for x in net.nodes()
]
# With omission
degree_weighted_exp_omit = [
    net_omit.degree(x, weight = f'{weight}_exp') 
    for x in net_omit.nodes()
]
degree_weighted_imp_omit = [
    net_omit.degree(x, weight = f'{weight}_imp') 
    for x in net_omit.nodes()
]

(sum(degree_weighted_exp) - sum(degree_weighted_exp_omit)) / sum(degree_weighted_exp) * 100
(sum(degree_weighted_imp) - sum(degree_weighted_imp_omit)) / sum(degree_weighted_imp) * 100

########
merged_data.columns
subset_m = merged_data.filter(pl.col('cmdCode').str.contains('4403'),
                              pl.col('period') == '2020',
                              pl.col('reporterDesc') == 'Japan',
                              pl.col('flowCode') == 'M')
subset_m.select(pl.col('FAO Code')).unique()

subset_4403_m = subset_m.filter(pl.col('cmdCode') == '4403')
subset_HS6_m = subset_m.filter(pl.col('cmdCode') != '4403')

subset_4403_m.select(pl.sum('primaryValue'))
subset_HS6_m.select(pl.sum('primaryValue'))

subset_m_fao = merged_data.filter(pl.col('FAO Code') == '012',
                                  pl.col('period') == '2020',
                                  pl.col('reporterDesc') == 'Japan',
                                  pl.col('flowCode') == 'M')
subset_m_fao.select(pl.sum('primaryValue'))

lst_HS6 = subset_HS6_m.select(pl.col('cmdCode')).unique().to_series().to_list()
lst_fao = subset_m_fao.select(pl.col('cmdCode')).unique().to_series().to_list()
lst_HS6.sort()
lst_fao.sort()
lst_HS6
lst_fao

################
########
input_data.columns
subset_i = input_data.filter(pl.col('fao_code') == "012",
                              pl.col('period') == '1996',
                              pl.col('reporter_desc') == 'Japan',
                              pl.col('flow_code') == 'M')
subset_i.select(pl.col('fao_code')).unique()
subset_i.select(pl.sum('primary_value'))
input_data.select('reporter_desc').unique().to_series().to_list()
input_data.filter(pl.col('reporter_desc') == 'Japan').select('period').unique().to_series().to_list()
merged_data.filter(pl.col('reporterDesc') == 'Japan').select('period').unique().to_series().to_list()

############
mirror_flows.columns
subset_mi = mirror_flows.filter(pl.col('fao_code') == 12,
                              pl.col('period') == 2020,
                              pl.col('importer_desc') == 'Japan')
subset_mi.select(pl.col('fao_code')).unique()
subset_mi.select(pl.sum('primary_value_exp'))

val = (mirror_flows.filter(pl.col('fao_code') == 12,
                     pl.col('period') == 1996)
             .select(pl.sum('primary_value_imp')))
val_omit = (mirror_flows.filter(pl.col('fao_code') == 12,
                     pl.col('period') == 1996,
                     pl.col('importer_desc') != 'Japan',
                     pl.col('exporter_desc') != 'Japan')
             .select(pl.sum('primary_value_imp')))
(val - val_omit) / val * 100


#################
def format_col_names(col_names: list) -> dict:
    '''
    Function that returns a mapping between column names and their formated 
    version as a dictionnary.

    Parameters
    ----------
    col_names : list of strings
        The list of column names to format.

    Returns
    -------
    col_mapping : dictionnary
        Mapping between column names and their formated version.
    '''

    # Insert an underscore before a uppercase
    col_names_formated = [re.sub(r"([a-z])([A-Z])", "\\1_\\2", s) for s in col_names]

    # Replace space by underscore
    col_names_formated = [re.sub(" ", "_", s) for s in col_names_formated]

    # Put every column name to lowercase
    col_names_formated = [s.lower() for s in col_names_formated]

    # Create mapping between column names and their formated version
    col_mapping = dict(zip(col_names, col_names_formated))

    return col_mapping

# Load data
merged_data

# Parameters
HS_version = [1996, 2002, 2007, 2012, 2017, 2022]
year_start = 2020
year_stop = 2024
fao_divisions = ['011', '012']
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

# Filter data for network analysis
input_data = (
    merged_data
        .select(col_keep)
        # Format column names for homogenity
        .rename(format_col_names(col_keep))
        .filter(
            # Keep data in specified time range
            pl.col('period') >= str(year_start),
            pl.col('period') <= str(year_stop - 2),

            # Keep imports and exports only
            pl.col('flow_code').is_in(flow_to_keep),

            # Keep FAO product division specified
            pl.col('fao_code').is_in(fao_divisions),

            # Remove non-country reporters and partners
            (~pl.col('reporter_iso')
             .str.contains('|'.join(excluded_iso))),
            (~pl.col('partner_iso')
             .str.contains('|'.join(excluded_iso))),

            # Delete "auto-" imports or exports (reporter = partner)
            pl.col('reporter_desc') != pl.col('partner_desc'),

            # # Remove primaryValue_deflated null values (~8% of the dataset)
            # ~pl.col('primary_value_deflated').is_null()
        )
        # Drop potential duplicates
        .unique()
)

# Drop outliers = values under fifth percentile for weight (kg) and value (USD)
stats_desc = (
    input_data
    .select(['net_wgt', 'primary_value'])
    .describe(percentiles=[0.05])
)

min_weight, min_value = (
    stats_desc
    .filter(pl.col('statistic') == '5%')
    .select(['net_wgt', 'primary_value'])
)

input_data = (
    input_data.filter(
        # Drop trade flow with net weight (kg) under fifth percentile
        pl.col('net_wgt') > min_weight.item(),

        # Drop trade flow with value (USD) under fifth percentile
        pl.col('primary_value') > min_value.item()
    )
)

# Sum weight and values by FAO division product
sum_cols = ['net_wgt', 'primary_value', 'primary_value_deflated']

groupby_cols = [_ for _ in input_data.columns if _ not in sum_cols]

input_data = (
    input_data
    .group_by(groupby_cols)
    .agg(pl.sum(sum_cols))
 )


subset_i = input_data.filter(pl.col('fao_code') == "012",
                              pl.col('period') == '1996',
                              pl.col('reporter_desc') == 'Japan',
                              pl.col('flow_code') == 'M')
subset_i.select(pl.col('fao_code')).unique()

subset_m_fao = merged_data.filter(pl.col('FAO Code') == '012',
                                  pl.col('period') == '1996',
                                  pl.col('reporterDesc') == 'Japan',
                                  pl.col('flowCode') == 'M',
                                  pl.col('partnerDesc') == 'World')

subset_i.select(pl.sum('primary_value'))
subset_m_fao.select(pl.sum('primaryValue'))

##############
year_start = 1996
year_stop = 2024
fao_divisions = ['011', '012']
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

test = (
    merged_data
        .select(col_keep)
        # Format column names for homogenity
        .rename(format_col_names(col_keep))
        .filter(
            # Keep data in specified time range
            pl.col('period') >= str(year_start),
            pl.col('period') <= str(year_stop - 2),

            # Keep imports and exports only
            pl.col('flow_code').is_in(flow_to_keep),

            # Keep FAO product division specified
            pl.col('fao_code').is_in(fao_divisions),

            # # Remove non-country reporters and partners
            # (~pl.col('reporter_iso')
            #  .str.contains('|'.join(excluded_iso))),
            # (~pl.col('partner_iso')
            #  .str.contains('|'.join(excluded_iso))),

            # Delete "auto-" imports or exports (reporter = partner)
            pl.col('reporter_desc') != pl.col('partner_desc'),

            # # Remove primaryValue_deflated null values (~8% of the dataset)
            # ~pl.col('primary_value_deflated').is_null()
        )
        # Drop potential duplicates
        .unique()
)
test.columns
(test.filter(pl.col('fao_code') == "012",
             pl.col('period') == '2020',
             pl.col('reporter_desc') == 'Japan',
             pl.col('flow_code') == 'M',
            #  pl.col('partner_desc') != 'World',
            #  pl.col('partner_desc') != 'Other Asia, nes',
             (~pl.col('partner_iso')
              .str.contains('|'.join(excluded_iso)))
            )
     .select('partner_desc')
     .unique()
     .to_series()
     .to_list())

(test.filter(pl.col('fao_code') == "012",
             pl.col('period') == '1996',
             pl.col('reporter_desc') == 'Japan',
             pl.col('flow_code') == 'M',
             pl.col('partner_desc') != 'World',
             pl.col('partner_desc') != 'Other Asia, nes',
            #  (~pl.col('partner_iso')
            #   .str.contains('|'.join(excluded_iso)))
             )
     .select(pl.sum('primary_value')))

merged_data.columns
exclusion = (
    merged_data.filter(
        pl.col('partnerISO').str.contains('|'.join(excluded_iso)),
        pl.col('reporterISO').str.contains('|'.join(excluded_iso))
        )
               .select('partnerDesc')
               .unique()
               .to_series()
               .to_list()
)

(test.filter(pl.col('fao_code') == "012",
             pl.col('period') == '1996',
             pl.col('reporter_desc') == 'Japan',
             pl.col('flow_code') == 'M',
             pl.col('partner_desc') != 'World',
             pl.col('partner_desc') != 'Other Asia, nes',
            #  (~pl.col('partner_desc').is_in(exclusion)),
            #  (~pl.col('partner_iso')
            #   .str.contains('|'.join(excluded_iso)))
             )
     .select(pl.sum('primary_value')))








######################
# Params
wb_series_drop      = ['TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'],
wb_countries_keep   = ['id', 'longitude', 'latitude','capitalCity']

# Load data
FAO_HS = pl.read_json('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/inhouse/correspondence_FAO_HS.json')

wb_countries = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/public/wb_countries_data.csv',
                      separator=';').select(['id', 'longitude', 'latitude','capitalCity'])

wb_series = (pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/public/wb_series_data.csv', separator=';')
             .with_columns(pl.col('time').str.replace(r'YR', ''))
             .drop(['TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD']))

deflate_uncomtrade = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/deflate_uncomtrade_data.parquet.gzip')

def join_uncomtrade_FAO(FAO_HS: pl.dataframe.frame.DataFrame, 
                        uncomtrade_data: pl.dataframe.frame.DataFrame):
    '''
    Function that joins the corresponding FAO forest product codes to the HS 
    codes in UNComtrade data.

    Parameters
    ----------
    FAO_HS : polars dataframe
        The corresponding table between HS codes and FAO forest products 
        classification.
    uncomtrade_data : polars dataframe
        The UNComtrade data (commodities in HS codes).

    Returns
    -------
    uncomtrade_data_join : polars dataframe
        The UN Comtrade data with joined FAO forest products codes corresponding
        to HS codes.
    '''

    # Collect HS classification version code in uncomtrade data (H0, H1, ...)
    HS_vers_code = (sorted(
        uncomtrade_data.select('classificationCode')
                        .unique()
                        .to_series()
                        .to_list()
        )
    )

    # Collect HS classification version year in FAO_HS corresponding table
    # N.B. 'HS 1996' is duplicated to match H0 and H1 in HS_vers_code
    HS_vers_year = (
        [code for code in FAO_HS.columns if '1996' in code] +
        [code for code in FAO_HS.columns if 'HS' in code]
    )

    # Collect FAO columns to join to uncomtrade data
    FAO_col = [col for col in FAO_HS.columns if 'HS' not in col]

    # Join FAO codes to uncomtrade data by batch of HS classification version
    uncomtrade_joined_batch = (
        [
            (uncomtrade_data
             .filter(pl.col('classificationCode') == HS_code)
             .join(FAO_HS.unique(subset=FAO_col+[HS_year]),
                   left_on='cmdCode',
                   right_on=HS_year,
                   how='left')
                   .drop(cs.starts_with('HS'))
            )
            for HS_code, HS_year in zip(HS_vers_code,HS_vers_year)
        ]
    )

    # Concatenate all batch
    uncomtrade_data_join = pl.concat(
        [df for df in uncomtrade_joined_batch if df.shape != (0,0)],
        how='vertical_relaxed'
    )

    return uncomtrade_data_join


# Join UNComtrade data with World Bank data
merged_data = (
    deflate_uncomtrade
    .join(wb_series,
          left_on=['reporterISO', 'period'],
          right_on=['economy', 'time'],
          how='left')
    .join(wb_countries,
          left_on=['reporterISO'],
          right_on=['id'],
          how='left')
)

# Join UNComtrade data with FAO_HS corresponding table
merged_data = join_uncomtrade_FAO(FAO_HS, merged_data)

