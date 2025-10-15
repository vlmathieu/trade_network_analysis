from snakemake.script import snakemake
import logging
import polars as pl
import pickle
import numpy as np
import networkx as nx

def unit_network_contribution(
        unit_edge_list_dict: dict,
        to_omit: str,
        weight: str = 'primary_value') -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of a country's contribution to 
    total traded value for exports and imports and total number of edges based 
    on an edge list describing a trade network for a given year and traded 
    product. The weight to be taken into account when calculating the total 
    traded value is given by the weight parameter. The network refers to the 
    trade of one year and one product and is directed and weighted.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (cmd, period) of the commodity
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network total traded value and market 
        concentration indices are calculated.
    to_omit : string
        The name of the country on which network contribution is calculated,
        i.e., the country to omit in order to measure its contribution in the
        netork total traded value and number of edges.
    weight : string
        The weight on which network total traded value and market concentration 
        indices are calculated. The default value is "primary_value".

    Returns
    -------
    unit_market_concentration : polars data frame
        A polars data frame of a country's contribution to total traded value 
        for exports and imports and total number of edges for a given year and 
        traded product.
        
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

    # Build dictionnary of tot traded value and market concentration indices
    unit_network_contribution = pl.from_dict(
            {
                "period": period,
                "cmd": cmd,
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

    return unit_network_contribution

def unit_network_contribution_all_nodes(
        unit_edge_list_dict: dict,
        weight: str = 'primary_value') -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of every country's contribution to 
    total traded value for exports and imports and total number of edges based 
    on an edge list describing a trade network for a given year and traded 
    product. The weight to be taken into account when calculating the total 
    traded value is given by the weight parameter. The network refers to the 
    trade of one year and one product and is directed and weighted.

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
    unit_market_concentration_all_nodes : polars data frame
        A polars data frame of every country's contribution to total traded 
        value for exports and imports and total number of edges for a given year
        and traded product.
        
    '''

    # Extract keys and edge list from dict
    [[_, edge_list]] = unit_edge_list_dict.items()

    # Build directed network based on edge_list
    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)

    # List trading countries
    countries = list(net.nodes())
    countries.sort()

    # Compute country's contribution for every trading country
    unit_market_concentration_all_nodes = pl.concat(
        [
            unit_network_contribution(
                unit_edge_list_dict = unit_edge_list_dict,
                to_omit = country,
                weight = weight
            )
            for country in countries
        ]
    )

    return unit_market_concentration_all_nodes

def network_contribution(
        edge_list_dict: dict, 
        weight: list = 'primary_value') -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of every country's contribution to 
    total traded value for exports and imports and total number of edges for 
    every year and commodity of trade based on a dictionary of edge lists 
    describing a trade network for each year and product of trade considered. 
    The weight to be taken into account when calculating the total traded value 
    is given by the weight parameter. The network is directed and weighted.

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
    network_contribution : polars data frame
        A polars data frame of every country's contribution to total traded 
        value for exports and imports and total number of edges for each year 
        and product considered.
        
    '''

    # Divide global dictionary into list of unit edge list dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Apply unit_network_contribution_all_nodes to every edge list dictionnary
    network_contribution = pl.concat(
        [
            unit_network_contribution_all_nodes(
                unit_edge_list_dict = unit_edge_list_dict, 
                weight = weight
            )
            for unit_edge_list_dict in edge_lists
        ],
        how = 'vertical_relaxed'
    )

    # Sort by cmd and year
    network_contribution = network_contribution.sort(['cmd', 'period'])
    
    return network_contribution

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute market concentration stats based on dictionnary of edge lists
network_contribution = network_contribution(
    edge_list_dict= edge_list_dict,
    weight= snakemake.params['weight']
)
logging.info(f"\nNetwork contribution:\n {network_contribution}\n")

# Save market concentration stats
network_contribution.write_csv(
    snakemake.output[0],
    separator=';'
)
