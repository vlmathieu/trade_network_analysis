from snakemake.script import snakemake
import logging
import polars as pl
import pickle
import networkx as nx

def unit_network_composition(
        unit_edge_list_dict: dict) -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of the network composition 
    descriptive statistics (number of trading countries, number of pure 
    exporters, number of pure importers, number of countries that are both 
    exporters and importers) based on an edge list describing a trade network 
    for a given year and traded product. The network refers to the trade of one
    year and one product and is directed and unweighted.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (cmd, period) of the commodity
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network composition descriptive statistics are
        calculated.

    Returns
    -------
    unit_network_composition : pl.dataframe.frame.DataFrame
        A polars data frame of the network composition descriptive statistics 
        for a given year of trade and traded product.
        
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
    
    # List sources (origin of edge = exporters)
    sources = [x for x in net.nodes() if net.out_degree(x) >= 1]

    # List targets (destination of edge = importers)
    targets = [x for x in net.nodes() if net.in_degree(x) >= 1]

    # Compute list of total number of trading countries
    tot_nb_nodes = list(set(sources + targets))

    # Compute lists of pure sources = exporters & of pure targets = importers
    pure_sources = list(set(sources) - set(targets))
    pure_targets = list(set(targets) - set(sources))

    # Compute list of countries that are both sources and targets = exp & imp
    mixed_src_tgt = list(set(sources).intersection(targets))

    # Compute network composition descriptive statistics
    unit_network_composition = pl.from_dict(
        {
            "period": period,
            "cmd": cmd,
            'tot_nb_nodes': len(tot_nb_nodes),  # Number of trading countries
            'nb_pure_exp': len(pure_sources),   # Number of pure exporters
            'nb_pure_imp': len(pure_targets),   # Number of pure importers
            'nb_mixed': len(mixed_src_tgt),     # Number of mixed countries
        }
    )

    return unit_network_composition

def network_composition(edge_list_dict: dict) -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of the network composition 
    descriptive statistics (number of trading countries, number of pure 
    exporters, number of pure importers, number of countries that are both 
    exporters and importers) based on a dictionary of edge lists describing a 
    trade network for each year and product of trade considered. The network is 
    directed and unweighted.

    Parameters
    ----------
    edge_list_dict : dictionnary
        A dictionnary that associates, for all years and products of trade 
        covered, (i) a tuple (product, year) of the product code and the year of 
        trade and (ii) the associated edge list describing the network and on 
        which network composition descriptive statistics are calculated.

    Returns
    -------
    network_composition : polars data frame
        A polars data frame of the network composition descriptive statistics 
        for each year of trade and traded product considered.
        
    '''

    # Divide global dictionary into list of unit edge list dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Apply unit_network_composition to every edge list dictionnary
    network_composition = pl.concat(
        [
            unit_network_composition(
                unit_edge_list_dict = unit_edge_list_dict
            )
            for unit_edge_list_dict in edge_lists
        ],
        how = 'vertical_relaxed'
    )

    # Sort by cmd and period
    network_composition = network_composition.sort(['cmd', 'period'])

    return network_composition

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute network composition desc stats based on dictionnary of edge lists
network_composition = network_composition(
    edge_list_dict= edge_list_dict
)
logging.info(f"\nNetwork composition:\n {network_composition}\n")

# Save network composition
network_composition.write_csv(
    snakemake.output[0],
    separator=';'
    )