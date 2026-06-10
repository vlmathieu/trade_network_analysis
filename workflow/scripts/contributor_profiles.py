from snakemake.script import snakemake
import logging
import polars as pl
import pickle
import networkx as nx # pyright: ignore[reportMissingModuleSource]

def unit_node_metrics(unit_edge_list_dict: dict) -> pl.DataFrame:
    '''
    Returns a Polars DataFrame of node-level metrics for all nodes in a
    directed weighted trade network for a given (cmd, period) tuple.

    Parameters
    ----------
    unit_edge_list_dict : dict
        A dict with one entry: {(cmd, period): edge_list}.

    Returns
    -------
    pl.DataFrame
        One row per node with degree, traded value, net weight,
        orientation ratios, and per-partner intensity columns.
    '''

    [[keys, edge_list]] = unit_edge_list_dict.items()
    cmd, period = keys

    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)

    for _, _, d in net.edges(data=True):
        for key in d:
            if d[key] is None:
                d[key] = 0

    rows = []
    for x in net.nodes():
        pv_exp  = net.out_degree(x, weight='primary_value_exp')
        pv_imp  = net.in_degree(x,  weight='primary_value_imp')
        wgt_exp = net.out_degree(x, weight='net_wgt_exp')
        wgt_imp = net.in_degree(x,  weight='net_wgt_imp')
        nb_exp  = net.out_degree(x)
        nb_imp  = net.in_degree(x)
        pv_tot  = pv_exp + pv_imp

        rows.append({
            "period": period,
            "cmd": cmd,
            "country": x,
            "nb_edge_exp": nb_exp,
            "nb_edge_imp": nb_imp,
            "primary_value_exp": pv_exp,
            "primary_value_imp": pv_imp,
            "net_wgt_exp": wgt_exp,
            "net_wgt_imp": wgt_imp,
            "exp_share": pv_exp / pv_tot if pv_tot > 0 else None,
            "imp_share": pv_imp / pv_tot if pv_tot > 0 else None,
            "primary_value_per_partner_exp": pv_exp  / nb_exp  if nb_exp  > 0 else None,
            "primary_value_per_partner_imp": pv_imp  / nb_imp  if nb_imp  > 0 else None,
            "net_wgt_per_partner_exp":       wgt_exp / nb_exp  if nb_exp  > 0 else None,
            "net_wgt_per_partner_imp":       wgt_imp / nb_imp  if nb_imp  > 0 else None,
        })

    return pl.from_dicts(rows)


def contributor_profiles(edge_list_dict: dict) -> pl.DataFrame:
    '''
    Returns a Polars DataFrame of node-level metrics for all nodes across
    all (cmd, period) combinations in the supplied edge list dictionary.
    No threshold filtering is applied here — country selection is delegated
    to the downstream plotting rule so the threshold can be tuned without
    re-running this rule.

    Parameters
    ----------
    edge_list_dict : dict
        {(cmd, period): edge_list, ...} for all years and products.

    Returns
    -------
    pl.DataFrame
        All node metrics sorted by cmd, period, country.
    '''

    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    result = pl.concat(
        [unit_node_metrics(d) for d in edge_lists],
        how='vertical_relaxed'
    )

    return result.sort(['cmd', 'period', 'country'])


logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

edge_list_dicts = []
for p in snakemake.input:
    with open(p, 'rb') as f:
        edge_list_dicts.append(pickle.load(f))

results = [contributor_profiles(edge_list_dict=d) for d in edge_list_dicts]

logging.info(f"\nContributor profiles:\n{results[0]}\n")
if len(results) > 1:
    logging.info(f"\nContributor profiles (agg EU):\n{results[1]}\n")

for data, path in zip(results, snakemake.output):
    data.write_csv(path, separator=';')
