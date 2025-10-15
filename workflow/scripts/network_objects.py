from snakemake.script import snakemake
import logging
import polars as pl
import pickle

def one_side_flows(data: pl.dataframe.frame.DataFrame, 
                   flow_code: str,
                   flow_type: list = ['wgt', 'value'],
                   trader_role: list = ['reporter', 'partner']):
    '''
    Function that extracts and formats one part of the mirror flows (either
    exporter reports or importer reports) and return one-side trade flows.

    Parameters
    ----------
    data : polars dataframe
        The UNComtrade data.
    flow_code : string
        Code that specifies if a trade flow report is an import or an export. 
        Must be 'M' for import or 'X' for export.
    flow_type : list of string
        Type of flows considered for building mirror flows. The items of 
        flow_type must be contained in column names. Default value is 
        ['wgt', 'value'] for commodities weight and trade values.
    trader_role : list of string
        The role of reporter or partner of traders in trade flow reporting.
        Default value is ['reporter', 'partner'].
    
    Raises
    ------
    ValueError
        Raises ValueError is flowDesc is not 'M' or 'X'.

    Returns
    -------
    extraction_formated : polars dataframe
        Either import- or export-side trade flows formated to build mirro flows.

    '''
    
    # Values of flow_code must be 'M' for import or 'X' for export, error is not
    valid = {'M','X'}
    if flow_code not in valid:
        raise ValueError('one_side_flow: flow_code must be one of %r' % valid)
    
    # List column names of flow types
    flows = [_ for _ in data.columns if any(s in _ for s in flow_type)]
    traders = [_ for _ in data.columns if any(s in _ for s in trader_role)]

    # Format one_side_flow columns names depending on flow_code
    if flow_code == 'X':
        mapping = (
            # If export, reporter is the exporter and partner the importer
            dict(zip(traders,
                     [_.replace('reporter', 'exporter')
                       .replace('partner', 'importer') for _ in traders])) |
            # Add suffix _exp to flow type columns
            dict(zip(flows, 
                     [f'{_}_exp' for _ in flows]))
        )
    elif flow_code == 'M':
        mapping = (
            # If import, reporter is the importer and partner the exporter
            dict(zip(traders,
                     [_.replace('reporter', 'importer')
                       .replace('partner', 'exporter') for _ in traders])) |
            # Add suffix _imp to flow type columns
            dict(zip(flows, 
                     [f'{_}_imp' for _ in flows]))
        )
    # Extract data and format columns
    extraction_formated = (
        data
        .filter(pl.col('flow_code') == flow_code)
        .rename(mapping)
        .drop('flow_code')
    )
    
    return extraction_formated

def compute_mirror_flows(data: pl.dataframe.frame.DataFrame,
                         report_suffix: list = ['_imp', '_exp']):
    '''
    Function that assembles one-side trade flows into mirror flows based on 
    cleaned UN Comtrade data.

    Parameters
    ----------
    data : polars dataframe
        The UNComtrade data.
    report_suffix : list of string
        List of suffixes contained in column names referring to country reports. 
        Default is ['_imp', '_exp'] for importer or exporter reports, 
        respectively.

    Returns
    -------
    mirror_flows : polars dataframe
        The trade mirror flows.

    '''
    
    exp_flow = one_side_flows(data,'X')
    imp_flow = one_side_flows(data,'M')
    join_col = [_ for _ in exp_flow.columns 
                if not any(s in _ for s in report_suffix)]

    mirror_flows = exp_flow.join(imp_flow,
                                 on=join_col,
                                 how='full',
                                 coalesce=True)
    
    return mirror_flows

def unit_edge_list(unit_mirror_flows: pl.dataframe.frame.DataFrame,
                   year_col: str = 'period',
                   cmd_col: str = 'fao_code',
                   country_desc: str = '_desc'):
    '''
    Functions that computes nodes, edge list, adjacency list, and adjacency 
    matrices (unweighted, and weighted by export or import reports) based on
    mirror flows of trade for one year and one product.

    Parameters
    ----------
    unit_mirror_flows : polars dataframe
        The mirror flows for one product and one year of trade.
    year_col : string
        Column name of years. Default value is 'period'.
    cmd_col : string
        Column name of commodity codes. Default value is 'fao_code'.
    country_desc : string
        Suffix of the exporter and importer columns (= nodes of the network) to 
        be considered . Can take values such as '_desc' (i.e., exporters and 
        importers are described by their name) or '_iso' (i.e., exporters and 
        importers are described by their ISO code). Default value is '_desc'.

    Raises
    ------
    ValueError
        Raises ValueError if mirror flows cover more than one product or year.

    Returns
    -------
    result : dictionary
        A dictionary of the edge list (weighted by export or import reports) 
        describing the trade network for a given product and year.
    '''
    
    # Check if input unit_mirror_flows cover only one product and one year
    if (len(unit_mirror_flows.select(year_col).unique()) > 1 
        or len(unit_mirror_flows.select(cmd_col).unique()) > 1):
        raise ValueError(
            'edge_list: data must refer to one year and one product only'
            )
    
    # Get product code
    product = int(unit_mirror_flows.select(cmd_col).unique().item())

    # Get year of trade
    year = int(unit_mirror_flows.select(year_col).unique().item())
    
    # Build weights dictionnary
    wgt_col = [_ for _ in unit_mirror_flows.columns 
               if any(s in _ for s in ['_exp', '_imp'])]
    wgt_dicts = unit_mirror_flows.select(wgt_col).to_dicts()

    # Build edge_list
    desc_col = [_ for _ in unit_mirror_flows.columns if country_desc in _]
    edge_list = [(exp, imp) 
                 for exp, imp in unit_mirror_flows.select(desc_col).rows()]

    # Add weights to edge_list
    edge_list_wgt = [
        edge + (wgt,)
        for edge, wgt in zip(edge_list, wgt_dicts)
    ]

    # Edit result to return
    result = {(product, year): edge_list_wgt}

    return result

def edge_lists(data: pl.dataframe.frame.DataFrame,
               year_col: str = 'period',
               cmd_col: str = 'fao_code',
               country_desc: str = '_desc'):
    '''
    Functions that computes edge list (weighted by export or import reports)
    describing the trade network for the different products and years covered 
    in the trade data and stores the results in a dictionary.

    Parameters
    ----------
    data : polars dataframe
        The trade data.
    year_col : string
        Column name of years. Default value is 'period'.
    cmd_col : string
        Column name of commodity codes. Default value is 'fao_code'.
    country_desc : string
        Suffix of the exporter and importer columns (= nodes of the network) to 
        be considered . Can take values such as '_desc' (i.e., exporters and 
        importers are described by their name) or '_iso' (i.e., exporters and 
        importers are described by their ISO code). Default value is '_desc'.

    Returns
    -------
    result : dictionary
        Dictionary of edge list describing the trade network for the different 
        products and years considered.

    '''
    
    # Build mirror flows based on trade data
    mirror = compute_mirror_flows(data)

    # Split mirror flows into a list of unit_mirror_flows (single year and cmd)
    split_mirror = [rows for _, rows in mirror.group_by([year_col, cmd_col])]
    
    # Build a list of edge_list by year and cmd, based on unit_mirror_flows
    result_list = (
        list(
            map(
                lambda df: unit_edge_list(df, 
                                          year_col=year_col, 
                                          cmd_col=cmd_col, 
                                          country_desc=country_desc), 
                                          split_mirror))
    )

    # Structure result as a dictionnary
    result = {k:v for elem in result_list for (k,v) in elem.items()}
    
    return result

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
# Load input trade data
input_data = pl.read_parquet(snakemake.input[0])

# Create mirror flows and edge lists
mirror_flows = compute_mirror_flows(input_data)
edge_lists = edge_lists(input_data)
logging.info(f"\nMirror flows:\n {mirror_flows}\n")

# Save mirror flows
mirror_flows.write_csv(
    snakemake.output[0],
    separator=';'
    )

# Save edge list dictionnary
with open(snakemake.output[1], 'wb') as f:
    pickle.dump(edge_lists, f)
