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
    Classification requires a single value per country, so each flow is
    measured by the country's OWN report (its exports as it reports them, its
    imports as it reports them), falling back on the partner's mirror report
    only when the own report is absent (0). Reports are never averaged.

    Flow volumes are never reconciled: the exporter-reported and importer-
    reported values of each flow are accumulated separately, because the gap
    between the two mirror reports is a data-quality signal to preserve, not
    noise to average away. Every flow statistic therefore comes in two variants
    along two ORTHOGONAL dimensions:
      - Attribution (src/tgt): which country's category the flow counts toward.
          - Source attribution: each flow is attributed to the category of its
            source country (exporter). Captures the supply-side concentration
            of trade, i.e. what share of world exports originates from each
            category.
          - Target attribution: each flow is attributed to the category of its
            target country (importer). Captures the demand-side concentration
            of trade, i.e. what share of world imports is absorbed by each
            category.
      - Report (_exp/_imp suffix): whose declaration measures the flow — the
        exporter's (FOB valuation for primary value) or the importer's (CIF).
        Note the suffix refers to the REPORTING side, not to the main_exp /
        main_imp categories.
    Crossing them yields four cases per category and measure; the gap between
    the _exp and _imp variants of a same attribution is the mirror-report
    discrepancy for that category. A zero report contributes zero to that
    report's series — missingness is part of the signal (no fallback here,
    unlike classification, which needs one number per country).

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

        Country lists (sorted alphabetically, serialised as strings in CSV):
          list_main_exp : Sorted list of main exporter country codes
          list_main_imp : Sorted list of main importer country codes
          list_balanced : Sorted list of balanced country codes

        Network totals (denominators for share computation), one per report:
          tot_net_wgt_exp        : Total net weight, as reported by exporters
          tot_net_wgt_imp        : Total net weight, as reported by importers
          tot_primary_value_exp  : Total primary value, as reported by exporters (FOB)
          tot_primary_value_imp  : Total primary value, as reported by importers (CIF)

        Source-attributed flow totals (supply-side perspective), for each
        category {main_exp|main_imp|balanced}, measure {net_wgt|primary_value}
        and report {exp|imp} (12 columns):
          src_{category}_{measure}_{report} : Flow volume originating from
            countries of that category, as measured by that report side

        Target-attributed flow totals (demand-side perspective), same crossing
        (12 columns):
          tgt_{category}_{measure}_{report} : Flow volume absorbed by countries
            of that category, as measured by that report side

        For each (measure, report) pair, the three src_* columns and the three
        tgt_* columns each sum to the corresponding tot_* column.

        Price-consistent value totals (for implied prices), for each attribution
        {src|tgt}, category and report (14 columns incl. tot):
          {src|tgt}_{category}_primary_value_wp_{report} / tot_primary_value_wp_{report}
            : primary value summed only over flows whose SAME-side net weight is
            reported (wp = "weight present"). Dividing these by the matching
            {..}_net_wgt_{report} totals gives a value-to-weight ratio taken over
            one flow set for numerator and denominator; the plain primary_value
            totals would instead inflate the ratio wherever net weight is missing.

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
 
    # --- Step 1: classify countries ---
    # Accumulate per-country export and import primary values to determine each
    # country's trade orientation. Primary value drives classification; net
    # weight is used only for the flow metrics in step 2.
    # Classification needs a single value per country, so each flow is measured
    # by the country's OWN report, falling back on the partner's mirror report
    # only when the own report is absent (0) — e.g. for non-reporting countries
    # that appear in the network solely through their partners' declarations.
    # Reports are never averaged.
    country_exp_value = {}
    country_imp_value = {}

    for src, tgt, d in net.edges(data=True):
        v_exp = d['primary_value_exp']
        v_imp = d['primary_value_imp']
        country_exp_value[src] = (country_exp_value.get(src, 0)
                                  + (v_exp if v_exp > 0 else v_imp))
        country_imp_value[tgt] = (country_imp_value.get(tgt, 0)
                                  + (v_imp if v_imp > 0 else v_exp))
 
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
 
    # --- Step 2: accumulate flow volumes by category, attribution and report ---
    # Initialise accumulators for source- and target-attributed flow totals,
    # for each measure, each report side and each country category. The
    # exporter- and importer-reported values are accumulated separately (no
    # reconciliation): a zero report contributes zero to that report's series,
    # so mirror-report gaps stay visible in the output.
    categories = ['main_exp', 'main_imp', 'balanced']
    measures   = ['net_wgt', 'primary_value']
    reports    = ['exp', 'imp']

    src_flows = {f'src_{cat}_{m}_{r}': 0.0
                 for cat in categories for m in measures for r in reports}
    tgt_flows = {f'tgt_{cat}_{m}_{r}': 0.0
                 for cat in categories for m in measures for r in reports}
    tot_flows = {f'tot_{m}_{r}': 0.0 for m in measures for r in reports}

    # Price-consistent value accumulators: primary value summed ONLY over flows
    # whose matching-side net weight is reported (> 0). Dividing these by the
    # net-weight totals gives each group's implied price (US$/tonne) over the
    # SAME set of flows for numerator and denominator. The plain primary_value
    # totals instead count flows whose net weight is missing (null, replaced by 0
    # above), which would inflate the ratio over the 2000-2006 reporting gap;
    # these masked columns are what plot_prices_figures.R (LOOP 1) divides by.
    # See the "Net-weight reporting gaps" supplementary section.
    src_price_val = {f'src_{cat}_primary_value_wp_{r}': 0.0
                     for cat in categories for r in reports}
    tgt_price_val = {f'tgt_{cat}_primary_value_wp_{r}': 0.0
                     for cat in categories for r in reports}
    tot_price_val = {f'tot_primary_value_wp_{r}': 0.0 for r in reports}

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
            for r in reports:
                flow = d[f'{m}_{r}']
                src_flows[f'src_{src_cat}_{m}_{r}'] += flow
                tgt_flows[f'tgt_{tgt_cat}_{m}_{r}'] += flow
                tot_flows[f'tot_{m}_{r}']           += flow

        # Price-consistent value: count value only where the same-side net
        # weight is reported (> 0), keeping numerator and denominator of the
        # implied price on one flow set (there are no genuine zero-weight flows,
        # so "> 0" isolates reported-weight flows from the null-set ones).
        for r in reports:
            if d[f'net_wgt_{r}'] > 0:
                v = d[f'primary_value_{r}']
                src_price_val[f'src_{src_cat}_primary_value_wp_{r}'] += v
                tgt_price_val[f'tgt_{tgt_cat}_primary_value_wp_{r}'] += v
                tot_price_val[f'tot_primary_value_wp_{r}']           += v

    # tot_flows is accumulated once per edge but independently for src and tgt,
    # so it is identical in both — use it as the single network total per
    # (measure, report) pair.

    # --- Step 3: assemble output dataframe ---
    unit_network_composition = pl.from_dict(
        {
            "period": period,
            "cmd": cmd,
            # Country counts
            'tot_nb_nodes': len(tot_nb_nodes),
            'nb_main_exp' : len(main_exporters),
            'nb_main_imp' : len(main_importers),
            'nb_balanced' : len(balanced),
            # Country lists (sorted alphabetically for reproducibility)
            'list_main_exp': [sorted(main_exporters)],
            'list_main_imp': [sorted(main_importers)],
            'list_balanced': [sorted(balanced)],
            # Network totals
            **tot_flows,
            **tot_price_val,
            # Source-attributed flow totals (supply-side perspective)
            **src_flows,
            **src_price_val,
            # Target-attributed flow totals (demand-side perspective)
            **tgt_flows,
            **tgt_price_val,
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
    data.with_columns([
        pl.col('list_main_exp').list.join('|'),
        pl.col('list_main_imp').list.join('|'),
        pl.col('list_balanced').list.join('|'),
    ]).write_csv(path, separator=';')