from snakemake.script import snakemake
import polars as pl
import polars.selectors as cs
import comtradeapicall
from tqdm import tqdm 

def get_uncomtrade(apikey: str, year: str, cmd: str, flow: str):
    '''
    Function that downloads UN Comtrade data for a given year, a given 
    commodity, and a given trade flow. Need an API key.

    Parameters
    ----------
    apikey : sting
        The API subscription key to download data.
    year : string
        The year of trade.
    cmd : string
        The commodity code.
    flow : string
        The trade flow to download (import, export, re-import, re-export...).

    Returns
    -------
    data : polars dataframe
        The UN Comtrade data for a given year, commodity, and trade flow.

    '''
    
    uncomtrade_data = comtradeapicall.getFinalData(
        apikey,
        typeCode        = 'C',          # typeCode(str) : Product type. Goods (C) or Services (S)
        freqCode        = 'A',          # freqCode(str) : The time interval at which observations occur. Annual (A) or Monthly (M)
        clCode          = 'HS',         # clCode(str) : Indicates the product classification used and which version (HS, SITC)
        period          = year,         # period(str) : Combination of year and month (for monthly), year for (annual)
        reporterCode    = None,         # reporterCode(str) : The country or geographic area to which the measured statistical phenomenon relates
        cmdCode         = cmd,          # cmdCode(str) : Product code in conjunction with classification code
        flowCode        = flow,         # flowCode(str) : Trade flow or sub-flow (exports, re-exports, imports, re-imports, etc.)
        partnerCode     = None,         # partnerCode(str) : The primary partner country or geographic area for the respective trade flow
        partner2Code    = None,         # partner2Code(str) : A secondary partner country or geographic area for the respective trade flow
        customsCode     = None,         # customsCode(str) : Customs or statistical procedure
        motCode         = None,         # motCode(str) : The mode of transport used when goods enter or leave the economic territory of a country
        format_output   = 'JSON',       # format_output(str) : The output format. CSV or JSON
        breakdownMode   = 'classic',    # breakdownMode(str) : Option to select the classic (trade by partner/product) or plus (extended breakdown) mode
        includeDesc     = True          # includeDesc(bool) : Option to include the description or not
        )
        
    return uncomtrade_data

def chunks(lst: list, n: int):
    '''
    Yield successive n-sized chunks from a list.

    Parameters
    ----------
    lst : list
        List to divide into chunks.
    n : integer
        Size of the chunks.
    '''
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_uncomtrade_bulk(apikey: str, years: list, cmdCode: list, flowCode: list):
    '''
    Function that downloads UN Comtrade data for a several years, commodities, 
    and trade flows. Need an API key.

    Parameters
    ----------
    apikey : string
        The API subscription key to download data.
    years : list of strings
        The years of trade.
    cmdCode : list of strings
        The commodity codes.
    flowCode : list of strings
        The trade flow to download (import, export, re-import, re-export...).

    Returns
    -------
    uncomtrade_data : polars dataframe
        The UN Comtrade data for a several years, commodity, and trade flows.

    '''

    # Download data by batch of 5 years and 3 commodities
    uncomtrade_years_batch = (
        [pl.from_pandas(
            get_uncomtrade(
                str(apikey),
                ','.join([str(_) for _ in years_batch]),
                ','.join([str(_) for _ in cmd_batch]),
                ','.join([str(_) for _ in flowCode])
            ))
        for years_batch in chunks(years, 5)
        for cmd_batch in chunks(cmdCode, 3)]
    )

    # Concatenate all batch
    uncomtrade_data = pl.concat(
        [df for df in uncomtrade_years_batch if df.shape != (0,0)],
        how='vertical_relaxed'
    )

    return uncomtrade_data

def process_param(start: int, 
                  stop: int, 
                  HS_versions: list, 
                  FAO_HS_json: str,
                  fao_divisions: list = ['012']):
    '''
    Function that process parameters for downloading uncomtrade data for wood 
    products, namely : years of download and HS codes to download, while taking 
    care of the evolution of Harmonized System classification through the years.

    Parameters
    ----------
    start : integer
        The first year of the trade period coverage.
    stop : integer
        The last year of the trade period coverage. This year is not considered
        for downloading (data will stop at stop-1).
    HS_versions : list of strings
        The year of the differents version of the Harmonized System: 1996, 2002,
        2007, 2012, 2017, 2022...
    FAO_HS_json : string
        Path to the json file of the correspondance between FAO and HS codes.
    fao_divisions : list of strings
        List of FAO division codes from which HS codes are derived. It 
        corresponds to the product of interest. Default value is ['012'] for 
        wood in the rough (roundwood).

    Returns
    -------
    zip(years, codes) : zip iterator of tuples
        The zip iterator of tuples of years under which the HS version remained
        unchanged and the corresponding HS codes for wood products.

    '''

    # HS classification versions used during the trading period covered
    HS_versions_keep = HS_versions[
        next(x[0] for x in enumerate(HS_versions) if x[1] > start)-1:
    ]

    # List of break years based on HS version and trade period coverage
    break_years = sorted(set(
        [_ for _ in HS_versions + [start, stop] if _ >= start]
    ))

    # Years over which HS version is unchanged
    years = []
    for i in range(0, len(break_years)-1):
        years.append(list(range(break_years[i], break_years[i+1])))

    # Codes for wood products under a common HS version for each year chunk
    FAO_HS_codes = (
        pl.read_json(FAO_HS_json)
          .filter(pl.col('FAO Code').is_in(fao_divisions)))
    codes = [cmdCode.to_list()
        for HS in HS_versions_keep
        for cmdCode 
        in FAO_HS_codes.select(cs.contains(str(HS))).unique()]
    
    return zip(years, codes)

def get_uncomtrade_all(apikey: str, start: int, stop: int, HS_versions: list, 
                       FAO_HS_json: str, flowCode: list):
    '''
    Function that downloads UNComtrade data for all wood products over a desired 
    trade period and for specific trade flows, taking into account the evolution 
    of the Harmonized System classification.

    Parameters
    ----------
    apikey : string
        The API subscription key to download data.
    start : integer
        The first year of the trade period coverage.
    stop : integer
        The last year of the trade period coverage. This year is not considered
        for downloading (data will stop at stop-1).
    HS_versions : list of strings
        The year of the differents version of the Harmonized System: 1996, 2002,
        2007, 2012, 2017, 2022...
    FAO_HS_json : string
        Path to a json file of the correspondance between FAO and HS codes.
    flowCode: list of strings
        The trade flow to download (import, export, re-import, re-export...).

    Returns
    -------
    uncomtrade_data : polars dataframe
        The UN Comtrade data for all wood products over a desired trade period 
        and for specific trade flows.
    '''

    # List years and codes for checks and tqdm
    years, codes = (
    map(list, 
        zip(*[list(_) 
              for _ in process_param(start, stop, HS_versions, FAO_HS_json)]))
    )

    # Dowload data by batch of HS version
    uncomtrade_batch = (
        [get_uncomtrade_bulk(apikey, years, cmdCode, flowCode)
         for years, cmdCode 
         in tqdm(process_param(start, 
                               stop, 
                               HS_versions, 
                               FAO_HS_json),
                               total=len(years))]
    )

    # Concatenate all batch
    uncomtrade_data = pl.concat(
        [df for df in uncomtrade_batch if df.shape != (0,0)],
        how='vertical_relaxed'
    )

    # Check of years, commodities, and different flows considered
    check_list = (
        sorted(set(uncomtrade_data['period'].unique())) == sorted(set([str(_) for years_ in years for _ in years_])),
        sorted(set(uncomtrade_data['cmdCode'].unique())) == sorted(set([str(_) for codes_ in codes for _ in codes_])),
        sorted(set(uncomtrade_data['flowCode'].unique())) == sorted(set([str(_) for _ in flowCode]))
    )

    # Checks for development (unnecessary)
    # print(sorted(set(uncomtrade_data['period'].unique())), '\n')
    # print(sorted(set([str(_) for years_ in years for _ in years_])), '\n')

    # print(sorted(set(uncomtrade_data['cmdCode'].unique())), '\n')
    # print(sorted(set([str(_) for codes_ in codes for _ in codes_])), '\n')

    # print(sorted(set(uncomtrade_data['flowCode'].unique())), '\n')
    # print(sorted(set([str(_) for _ in flowCode])), '\n')

    return uncomtrade_data, check_list

UN_Comtrade_data, check_list = get_uncomtrade_all(
    snakemake.params['apikey'],
    snakemake.params['year_start'],
    snakemake.params['year_stop'],
    snakemake.params['HS_version'],
    snakemake.input[0],
    snakemake.params['flowCode']
)

print("\nDataframe head: \n\n", UN_Comtrade_data.head(5), "\n")
print("\nDataframe size (rows, columns): ", UN_Comtrade_data.shape, "\n")

# Save data if check list passed
if all(check_list):
    print('Data have been checked.\n')   
    UN_Comtrade_data.write_parquet(
        snakemake.output[0],
        compression='gzip'
        )
else:
    print('Issues found in data download.\n')
