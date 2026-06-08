from snakemake.script import snakemake
import logging
import polars as pl
import pickle
import networkx as nx  # pyright: ignore[reportMissingModuleSource]


def _replace_none_weights(graph: nx.DiGraph) -> None:
    """Replace None edge attributes with 0 in-place."""
    for _, _, d in graph.edges(data=True):
        for key in d:
            if d[key] is None:
                d[key] = 0


def network_contribution_single(
        edge_list: list,
        cmd,
        period,
        weight: str = 'primary_value') -> pl.DataFrame:
    """
    Compute every country's contribution to total traded value and edge count
    for one (cmd, period) network.

    The contribution of a country is the reduction in total trade value that
    would result from removing it from the network. Analytically this equals:

      contrib_exp(c) = out_degree(c, weight_exp)          # c's own exports
                     + sum of weight_exp on c's in-edges  # partners lose exports to c
      contrib_imp(c) = in_degree(c, weight_imp)           # c's own imports
                     + sum of weight_imp on c's out-edges # partners lose imports from c
      contrib_edges(c) = out_degree(c) + in_degree(c)     # all incident edges

    Network-level totals (tot_exp, tot_imp, nb_edges) are computed once per
    network rather than once per country, giving an O(N + E) algorithm instead
    of O(N * (N + E)).

    Parameters
    ----------
    edge_list : list
        Edge list for this (cmd, period) network.
    cmd :
        Commodity code.
    period :
        Year.
    weight : str
        Edge attribute prefix ('primary_value' or 'net_wgt').

    Returns
    -------
    pl.DataFrame
        One row per country with columns: period, cmd, weight, country,
        nb_edges, nb_edges_country, tot_exp, tot_imp, tot_exp_country,
        tot_imp_country, contrib_nb_edges, contrib_tot_exp, contrib_tot_imp.
    """
    # Build directed graph
    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)
    _replace_none_weights(net)

    w_exp = f'{weight}_exp'
    w_imp = f'{weight}_imp'

    # Compute network-level totals once
    nb_edges = net.number_of_edges()
    tot_exp  = float(sum(net.out_degree(x, weight=w_exp) for x in net.nodes()))
    tot_imp  = float(sum(net.in_degree(x,  weight=w_imp) for x in net.nodes()))

    # Compute per-country contributions analytically
    countries = sorted(net.nodes())
    rows = {k: [] for k in [
        'period', 'cmd', 'weight', 'country',
        'nb_edges', 'nb_edges_country',
        'tot_exp', 'tot_imp',
        'tot_exp_country', 'tot_imp_country',
        'contrib_nb_edges', 'contrib_tot_exp', 'contrib_tot_imp',
    ]}

    for country in countries:
        # Contribution to edge count: all edges incident to this country
        nb_edges_country = net.out_degree(country) + net.in_degree(country)

        # Contribution to export value:
        #   own outgoing exports + in-neighbours lose their outgoing export to country
        tot_exp_country = float(
            net.out_degree(country, weight=w_exp)
            + sum(d.get(w_exp, 0) or 0
                  for _, _, d in net.in_edges(country, data=True))
        )

        # Contribution to import value:
        #   own incoming imports + out-neighbours lose their incoming import from country
        tot_imp_country = float(
            net.in_degree(country, weight=w_imp)
            + sum(d.get(w_imp, 0) or 0
                  for _, _, d in net.out_edges(country, data=True))
        )

        rows['period'].append(period)
        rows['cmd'].append(str(cmd))
        rows['weight'].append(weight)
        rows['country'].append(country)
        rows['nb_edges'].append(nb_edges)
        rows['nb_edges_country'].append(nb_edges_country)
        rows['tot_exp'].append(tot_exp)
        rows['tot_imp'].append(tot_imp)
        rows['tot_exp_country'].append(tot_exp_country)
        rows['tot_imp_country'].append(tot_imp_country)
        rows['contrib_nb_edges'].append(
            nb_edges_country / nb_edges if nb_edges else 0.0)
        rows['contrib_tot_exp'].append(
            tot_exp_country / tot_exp if tot_exp else 0.0)
        rows['contrib_tot_imp'].append(
            tot_imp_country / tot_imp if tot_imp else 0.0)

    return pl.from_dict(rows)


def network_contribution(
        edge_list_dict: dict,
        weight: str = 'primary_value') -> pl.DataFrame:
    """
    Compute every country's contribution for all (cmd, period) networks in
    edge_list_dict, for a single weight type.

    Parameters
    ----------
    edge_list_dict : dict
        Maps (cmd, period) tuples to edge lists.
    weight : str
        Edge attribute prefix ('primary_value' or 'net_wgt').

    Returns
    -------
    pl.DataFrame
        Concatenation of per-network results, sorted by cmd and period.
    """
    frames = [
        network_contribution_single(el, cmd, period, weight)
        for (cmd, period), el in edge_list_dict.items()
    ]
    return pl.concat(frames, how='vertical_relaxed').sort(['cmd', 'period'])


# ---------------------------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------------------------

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)

# Load all edge list dictionaries (one per aggregation level)
edge_list_dicts = []
for p in snakemake.input:
    with open(p, 'rb') as f:
        edge_list_dicts.append(pickle.load(f))

# Compute contributions for each aggregation level across all weights
results = []
for edge_list_dict in edge_list_dicts:
    frames = [
        network_contribution(edge_list_dict=edge_list_dict, weight=wgt)
        for wgt in snakemake.params['weight']
    ]
    results.append(pl.concat(frames, how='vertical_relaxed'))

logging.info(f"\nNetwork contribution country level:\n{results[0]}\n")
if len(results) > 1:
    logging.info(f"\nNetwork contribution aggregated EU:\n{results[1]}\n")

# Save results
for data, path in zip(results, snakemake.output):
    data.write_csv(path, separator=';')