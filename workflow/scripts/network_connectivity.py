from snakemake.script import snakemake
import logging
import polars as pl
import pickle
import numpy as np
import networkx as nx
from scipy.stats import kurtosis
from scipy.stats import skew


def unit_network_connectivity(
        unit_edge_list_dict: dict) -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of the network connectivity metrics 
    (total number of edges, mean, variance, skewness and kurtosis of the number
    of edge) for exports and imports based on an edge list describing a trade 
    network for a given year and traded product. The network refers to the trade
    of one year and one product and is directed and unweighted.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (cmd, period) of the commodity
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network connectivity metrics are calculated.

    Returns
    -------
    unit_network_connectivity : pl.dataframe.frame.DataFrame
        A polars data frame of the network connectivity metrics for exports and
        imports for a given year of trade and traded product.
        
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

    # Build unweighted degree lists for exp=out_degree | imp=in_degree
    degree_unweighted_exp = [
        net.out_degree(x) for x in net.nodes()
        # Consider only non-None edges = > 0 degree
        if net.out_degree(x) >= 1]
    degree_unweighted_imp = [
        net.in_degree(x) for x in net.nodes()
        # Consider only non-None edges = > 0 degree
        if net.in_degree(x) >= 1
    ]

    # Compute network statistics based on degree_list
    unit_network_connectivity = pl.from_dict(
        {
            "period": period,
            "cmd": cmd,
            # Compute total number of edges for export and import
            'nb_edge_exp': sum(degree_unweighted_exp),
            'nb_edge_imp': sum(degree_unweighted_imp),
            # Compute average number of edges for export and import
            'mean_exp': np.mean(degree_unweighted_exp),
            'mean_imp': np.mean(degree_unweighted_imp),
            # Compute variance of the number of edges for export and import
            'var_exp': np.var(degree_unweighted_exp),
            'var_imp': np.var(degree_unweighted_imp),
            # Compute skewness of the number of edges for export and import
            'skew_exp': skew(degree_unweighted_exp),
            'skew_imp': skew(degree_unweighted_imp),
            # Compute kurtosis of the number of edges for export and import
            'kurt_exp': kurtosis(degree_unweighted_exp),
            'kurt_imp': kurtosis(degree_unweighted_imp)
        }
    )

    return unit_network_connectivity

def network_connectivity(
        edge_list_dict: dict) -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of the network connectivity metrics 
    (total number of edges, mean, variance, skewness and kurtosis of the number
    of edge) for exports and imports based on a dictionary of edge lists 
    describing a trade network for each year and product of trade considered. 
    The network is directed and unweighted.

    Parameters
    ----------
    edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (cmd, period) of the commodity
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network connectivity metrics are calculated.

    Returns
    -------
    network_connectivity : pl.dataframe.frame.DataFrame
        A polars data frame of the network connectivity metrics for exports and
        imports for each year of trade and traded product considered.
        
    '''

    # Divide global dictionary into list of unit dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Apply unit_network_connectivity to every edge list dictionnary
    network_connectivity = pl.concat(
        [
            unit_network_connectivity(
                unit_edge_list_dict = unit_edge_list_dict
            )
            for unit_edge_list_dict in edge_lists
        ],
        how = 'vertical_relaxed'
    )

    # Sort by cmd and year
    network_connectivity = network_connectivity.sort(['cmd', 'period'])
    
    return network_connectivity

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute network connectivity metrics based on dictionnary of edge lists
network_connectivity = network_connectivity(
    edge_list_dict= edge_list_dict
)
logging.info(f"\nNetwork connectivity:\n {network_connectivity}\n")

# Save network_connectivity metrics
network_connectivity.write_csv(
    snakemake.output[0],
    separator=';'
    )
