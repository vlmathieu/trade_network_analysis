from snakemake.script import snakemake
import polars as pl
import pickle
import numpy as np
import networkx as nx

def unit_market_concentration(
        unit_edge_list_dict: dict, 
        weight: str = 'primary_value') -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of total traded value and market 
    concentration indices for exports and imports based on an edge list 
    describing a trade network for a given year and traded product. The weight 
    to be taken into account when calculating the total traded value and market 
    concentration indices is given by the weight parameter. The network refers 
    to the trade of one year and one product and is directed and weighted.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (cmd, period) of the commodity
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network total traded value and market 
        concentration indices are calculated.
    weight : string
        The weight on which network total traded value and market concentration 
        indices are calculated. The default value is "primary_value".

    Returns
    -------
    unit_market_concentration : polars data frame
        A polars data frame of the total traded value and market concentration
        indices for exports and imports for a given year and traded product.
        
    '''

    # Extract keys and edge list from dict
    [[keys, edge_list]] = unit_edge_list_dict.items()

    # Extract commodity code (=cmd) and year (=period)
    cmd, period = keys

    # Build directed network based on edge_list
    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)
    
    # Replace None weights by 0
    for _,_,d in net.edges(data=True):
        for key in d:
            if d[key] is None:
                d[key] = 0
    
    # Build weighted degree lists for exporters=out_degree | importers=in_degree
    degree_weighted_exp = [
        net.out_degree(x, weight = f'{weight}_exp') for x in net.nodes() 
        # Consider only non-None edges = > 0 weights
        if net.out_degree(x, weight = f'{weight}_exp') > 0
    ]
    degree_weighted_imp = [
        net.in_degree(x, weight = f'{weight}_imp') for x in net.nodes() 
        # Consider only non-None edges = > 0 weights
        if net.in_degree(x, weight = f'{weight}_imp') > 0
    ]

    # Compute total traded value for exports and imports
    traded_value_exp = sum(degree_weighted_exp)
    traded_value_imp = sum(degree_weighted_imp)

    # Build dictionnary of tot traded value and market concentration indices
    unit_market_concentration = pl.from_dict(
            {
                "period": period,
                "cmd": cmd,
                # Assign total circulating value for exports and imports
                "traded_value_exp": traded_value_exp,
                "traded_value_imp": traded_value_imp,
                # Compute Herfindahl-Hirschmann index for exports and imports
                "hhi_exp": sum([(x/traded_value_exp)**2 
                                for x in degree_weighted_exp]),
                "hhi_imp": sum([(x/traded_value_imp)**2 
                                for x in degree_weighted_imp]),
                # Compute Shannon index for exports and imports
                'shannon_exp': -sum(
                    [((x / traded_value_exp) * np.log(x / traded_value_exp)) 
                     for x in degree_weighted_exp if x > 0]),
                'shannon_imp': -sum(
                    [((x / traded_value_imp) * np.log(x / traded_value_imp)) 
                     for x in degree_weighted_imp if x > 0])
            }
    )

    return unit_market_concentration

def market_concentration(
        edge_list_dict: dict, 
        weight: list = 'primary_value') -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of total traded value and market 
    concentration indices for exports and imports based on a dictionary of edge 
    lists describing a trade network for each year and product of trade 
    considered. The weight to be taken into account when calculating the total 
    traded value and market concentration indices is given by the weight 
    parameter. The network is directed and weighted.

    Parameters
    ----------
    edge_list_dict : dictionnary
        A dictionnary that associates, for all years and products of trade 
        covered, (i) a tuple (product, year) of the product code and the year of 
        trade and (ii) the associated edge list describing the network and on 
        which network total traded value and market concentration indices are 
        calculated.
    weight : string
        The weight on which network total traded value and market concentration 
        indices are calculated. The default value is "primary_value".

    Returns
    -------
    market_concentration : polars data frame
        A polars data frame of the total traded value and market concentration
        indices for exports and imports for each year and product considered.
        
    '''

    # Divide global dictionary into list of unit edge list dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Apply unit_market_concentration to every edge list dictionnary
    market_concentration = pl.concat(
        [
            unit_market_concentration(
                unit_edge_list_dict = unit_edge_list_dict, 
                weight = weight
            )
            for unit_edge_list_dict in edge_lists
        ],
        how = 'vertical_relaxed'
    )

    # Sort by cmd and year
    market_concentration = market_concentration.sort(['cmd', 'period'])
    
    return market_concentration

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute market concentration stats based on dictionnary of edge lists
market_concentration = market_concentration(
    edge_list_dict= edge_list_dict,
    weight= snakemake.params['weight']
)

# Save market concentration stats
market_concentration.write_csv(
    snakemake.output[0],
    separator=';'
    )