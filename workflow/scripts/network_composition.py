from snakemake.script import snakemake
import logging
import polars as pl
import pickle
import networkx as nx # pyright: ignore[reportMissingModuleSource]

def unit_network_composition(
        unit_edge_list_dict: dict,
        threshold: float = 0.8) -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of the network composition 
    descriptive statistics based on an edge list describing a trade network for 
    a given year and traded product. The network refers to the trade of one year 
    and one product and is directed and weighted.
 
    Countries are classified by the share of their total traded value (exports 
    + imports) that is attributable to exports. A country is a "main exporter" 
    if its export share meets or exceeds the threshold, a "main importer" if its 
    import share meets or exceeds the threshold, and "balanced" otherwise.
 
    For each bilateral flow, all volume measures (net weight, primary value) are 
    estimated as the average of the exporter-reported and importer-reported 
    values when both are available, and as the single available value otherwise.
    This reconciles mirror trade data from both reporters.
 
    Flow volumes are attributed to country categories in two ways:
      - Source attribution: each flow is attributed to the category of its 
        source country (exporter). Captures the supply-side concentration of 
        trade, i.e. what share of world exports originates from each category.
      - Target attribution: each flow is attributed to the category of its 
        target country (importer). Captures the demand-side concentration of 
        trade, i.e. what share of world imports is absorbed by each category.
 
    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (cmd, period) of the commodity
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network composition descriptive statistics are
        calculated.
    threshold : float, optional
        Minimum share of total traded value (exports + imports) attributable to 
        exports (or imports) for a country to be classified as a main exporter 
        (or main importer). Must be between 0 and 1. Default is 0.8.
 
    Returns
    -------
    unit_network_composition : pl.dataframe.frame.DataFrame
        A polars data frame with the following metrics for a given year of trade 
        and traded product:
        
        Country counts:
          tot_nb_nodes  : Total number of trading countries
          nb_main_exp   : Number of main exporters
          nb_main_imp   : Number of main importers
          nb_balanced   : Number of balanced countries
 
        Network totals (denominators for share computation):
          tot_net_wgt        : Total net weight traded in the network
          tot_primary_value  : Total primary value traded in the network
 
        Source-attributed flow totals (supply-side perspective):
          src_main_exp_net_wgt       : Net weight of flows originating from main exporters
          src_main_imp_net_wgt       : Net weight of flows originating from main importers
          src_balanced_net_wgt       : Net weight of flows originating from balanced countries
          src_main_exp_primary_value : Primary value of flows originating from main exporters
          src_main_imp_primary_value : Primary value of flows originating from main importers
          src_balanced_primary_value : Primary value of flows originating from balanced countries
 
        Target-attributed flow totals (demand-side perspective):
          tgt_main_exp_net_wgt       : Net weight of flows absorbed by main exporters
          tgt_main_imp_net_wgt       : Net weight of flows absorbed by main importers
          tgt_balanced_net_wgt       : Net weight of flows absorbed by balanced countries
          tgt_main_exp_primary_value : Primary value of flows absorbed by main exporters
          tgt_main_imp_primary_value : Primary value of flows absorbed by main importers
          tgt_balanced_primary_value : Primary value of flows absorbed by balanced countries
        
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
 
    # --- Step 1: compute reconciled flow volumes for each edge ---
    # For each bilateral flow, reconcile mirror trade data by averaging the
    # exporter- and importer-reported values when both are non-zero, and using
    # the single available value otherwise. Done for both net weight and primary
    # value. Results are stored as new edge attributes for reuse in steps 2 & 3.
    for src, tgt, d in net.edges(data=True):
        for measure in ['net_wgt', 'primary_value']:
            v_exp = d[f'{measure}_exp']
            v_imp = d[f'{measure}_imp']
            if v_exp > 0 and v_imp > 0:
                d[f'{measure}_reconciled'] = (v_exp + v_imp) / 2
            else:
                d[f'{measure}_reconciled'] = v_exp + v_imp  # One of the two is 0
 
    # --- Step 2: classify countries ---
    # Accumulate per-country export and import primary values to determine each
    # country's trade orientation. Primary value drives classification; net
    # weight is used only for the flow metrics in step 3.
    country_exp_value = {}
    country_imp_value = {}
 
    for src, tgt, d in net.edges(data=True):
        flow_value = d['primary_value_reconciled']
        country_exp_value[src] = country_exp_value.get(src, 0) + flow_value
        country_imp_value[tgt] = country_imp_value.get(tgt, 0) + flow_value
 
    # Compute list of total number of trading countries
    tot_nb_nodes = list(set(country_exp_value) | set(country_imp_value))
 
    # Classify each country based on its export share of total traded value
    main_exporters = set()
    main_importers = set()
    balanced = set()
 
    for country in tot_nb_nodes:
        exp_val = country_exp_value.get(country, 0)
        imp_val = country_imp_value.get(country, 0)
        total_val = exp_val + imp_val
 
        if total_val == 0:
            # Country appears in the network but has no recorded trade value;
            # classify as balanced to avoid zero-division
            balanced.add(country)
        else:
            exp_share = exp_val / total_val
            if exp_share >= threshold:
                main_exporters.add(country)
            elif exp_share <= (1 - threshold):  # Equivalent: imp_share >= threshold
                main_importers.add(country)
            else:
                balanced.add(country)
 
    # --- Step 3: accumulate flow volumes by category and attribution ---
    # Initialise accumulators for source- and target-attributed flow totals,
    # for each measure and each country category.
    categories = ['main_exp', 'main_imp', 'balanced']
    measures   = ['net_wgt', 'primary_value']
 
    src_flows = {f'src_{cat}_{m}': 0.0 for cat in categories for m in measures}
    tgt_flows = {f'tgt_{cat}_{m}': 0.0 for cat in categories for m in measures}
    tot_flows = {f'tot_{m}': 0.0 for m in measures}
 
    # Helper: map a country to its category label
    def get_category(country: str) -> str:
        if country in main_exporters:
            return 'main_exp'
        elif country in main_importers:
            return 'main_imp'
        else:
            return 'balanced'
 
    for src, tgt, d in net.edges(data=True):
        src_cat = get_category(src)
        tgt_cat = get_category(tgt)
 
        for m in measures:
            flow = d[f'{m}_reconciled']
            src_flows[f'src_{src_cat}_{m}'] += flow
            tgt_flows[f'tgt_{tgt_cat}_{m}'] += flow
            tot_flows[f'tot_{m}']            += flow
 
    # tot_flows is accumulated once per edge but independently for src and tgt,
    # so it is identical in both — use it as the single network total.
 
    # --- Step 4: assemble output dataframe ---
    unit_network_composition = pl.from_dict(
        {
            "period": period,
            "cmd": cmd,
            # Country counts
            'tot_nb_nodes': len(tot_nb_nodes),
            'nb_main_exp' : len(main_exporters),
            'nb_main_imp' : len(main_importers),
            'nb_balanced' : len(balanced),
            # Network totals
            **tot_flows,
            # Source-attributed flow totals (supply-side perspective)
            **src_flows,
            # Target-attributed flow totals (demand-side perspective)
            **tgt_flows,
        }
    )
 
    return unit_network_composition

# def unit_network_composition(
#         unit_edge_list_dict: dict) -> pl.dataframe.frame.DataFrame:
#     '''
#     Function that returns a polar data frame of the network composition 
#     descriptive statistics (number of trading countries, number of pure 
#     exporters, number of pure importers, number of countries that are both 
#     exporters and importers) based on an edge list describing a trade network 
#     for a given year and traded product. The network refers to the trade of one
#     year and one product and is directed and unweighted.

#     Parameters
#     ----------
#     unit_edge_list_dict : dictionnary
#         A dictionnary that associates (i) a tuple (cmd, period) of the commodity
#         code and the year of trade and (ii) the associated edge list describing 
#         the network and on which network composition descriptive statistics are
#         calculated.

#     Returns
#     -------
#     unit_network_composition : pl.dataframe.frame.DataFrame
#         A polars data frame of the network composition descriptive statistics 
#         for a given year of trade and traded product.
        
#     '''

#     # Extract keys and edge list from dict
#     [[keys, edge_list]] = unit_edge_list_dict.items()

#     # Extract commodity code (=cmd) and year (=period)
#     cmd, period = keys

#     # Build directed network based on edge_list
#     net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)

#     # Replace None weights by 0
#     for _,_,d in net.edges(data=True):
#         for key in d:
#             if d[key] is None:
#                 d[key] = 0
    
#     # List sources (origin of edge = exporters)
#     sources = [x for x in net.nodes() if net.out_degree(x) >= 1]

#     # List targets (destination of edge = importers)
#     targets = [x for x in net.nodes() if net.in_degree(x) >= 1]

#     # Compute list of total number of trading countries
#     tot_nb_nodes = list(set(sources + targets))

#     # Compute lists of pure sources = exporters & of pure targets = importers
#     pure_sources = list(set(sources) - set(targets))
#     pure_targets = list(set(targets) - set(sources))

#     # Compute list of countries that are both sources and targets = exp & imp
#     mixed_src_tgt = list(set(sources).intersection(targets))

#     # Compute network composition descriptive statistics
#     unit_network_composition = pl.from_dict(
#         {
#             "period": period,
#             "cmd": cmd,
#             'tot_nb_nodes': len(tot_nb_nodes),  # Number of trading countries
#             'nb_pure_exp': len(pure_sources),   # Number of pure exporters
#             'nb_pure_imp': len(pure_targets),   # Number of pure importers
#             'nb_mixed': len(mixed_src_tgt),     # Number of mixed countries
#         }
#     )

#     return unit_network_composition

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
edge_list_dicts = []
for p in snakemake.input:
    with open(p, 'rb') as f:
        edge_list_dicts.append(pickle.load(f))

# Compute network composition desc stats based on dictionnary of edge lists
network_composition = [
    network_composition(
        edge_list_dict= edge_list_dict
    ) 
    for edge_list_dict in edge_list_dicts
]
logging.info(f"\nNetwork composition:\n {network_composition}\n")

# Save network composition
for data, path in zip(network_composition, snakemake.output):
    data.write_csv(path, separator=';')
