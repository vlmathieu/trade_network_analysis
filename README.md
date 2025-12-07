# Rise (and Fall?) of China's Demand for Imported Timber: A Network Analysis of the International Roundwood Trade

## Paper in brief

<!-- **Authors**: Valentin Mathieu, David Shanafelt -->

**Abstract**:

The international trade of wood products is an increasingly complex and interdependent component of the global economy, deeply linked to environmental issues, geopolitical shifts, changing competitive landscape, and macroeconomic trends.
While the literature offers comprehensive studies on the international timber trade, few to date consider the physical structure of the trade network in their analyses, leaving its structural evolution poorly understood.
To address this gap, we employ a network-theoretic approach to model global roudnwood trade flows between countries, providing a unique diagnosis of the trade's key structural properties and their temporal evolution from 1996 to 2023.
The analysis reveals that the network has developed short-term resilience while becoming consistently more interconnected and structurally concentrated over the long term.
Crucially, the structure is highly polarized around China's market dominance, which rapidly grew to concentrate over 56% of total trade value by 2021, due to its import demand.
This centralization, coupled with the structural withdrawal of major exporter Russia, created a single point of vulnerability in the global supply chain.
The recent decline in China's contribution post-2021—a result of geopolitical shocks, domestic self-sufficiency policies, and structural economic shifts—signals a fundamental shift in the international roundwood trade.
These findings underscore the need for exporters to implement market diversification strategies to mitigate policy risks in the global roundwood value chain.

**Keywords**: Globalization; Trade; Wood products; Network theory; Structural resilience; Market concentration; Forest policy

## Information

<!-- ### Funding statement

This study was supported by two grants from the French Grand Est Region (N° 19_GE8_019 20P05044) and the Lab of Excellence ARBRE (N° WP4/20PN17), by the French National Research Institute for Agriculture, Food and Environment (INRAE), and by AgroParisTech.

### CRediT authorship contribution statement

Valentin Mathieu: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data Curation, Writing - Original Draft, Writing - Review & Editing, Visualization, Supervision, Project administration. 

David W. Shanafelt: Conceptualization, Methodology, Formal analysis, Investigation, Writing - Original Draft, Writing - Review & Editing, Supervision. -->

### Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

### Declaration of generative AI and AI-assisted technologies in the writing process

During the preparation of this work the authors used (1) Asta in order to review the litterature aside traditional litterature search tools such as Scopus and Google Scholar and (2) Gemini and DeepL tools to improve English writing. After using this tool/service, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication.

### Data availability

Data and code are available on [GitHub](https://github.com/vlmathieu/trade_network_analysis).

<!-- ### Acknowledgements

We thank Felix Bastit and Clément Nedoncelle for their helpful and relevant comments and suggestions, as well as the following audience: BETA internal environment, Labex ARBRE, Chaire RENEL, the 9th annual conference of the French Association of Environmental and Resource Economists (FAERE), DEEPSURF 2022 Energy and Ecological transition international conference, and the 16th social science research days (JRSS), Paris-Saclay Conference on Trade and Environment, 26th World Congress of the International Union of Forest Research Organizations (IUFRO) - Forests and Society towards 2050. All remaining errors are our own. -->

## Folder structure
The folder is organized as follows:

```bash
├── config
│   ├── config.yaml
│   └── requirements.txt
├── dag.pdf
├── preprint.pdf
├── README.md
├── resources
│   ├── inhouse
│   └── public
├── results
│   ├── global
│   └── network_analysis
│       ├── input
│       ├── intermediary
│       ├── output
│       └── plot
└── workflow
    ├── envs
    ├── rules
    ├── scripts
    └── Snakefile
```

README.md provides information on the repository structuration and explains the data analysis workflow.
dag.pdf presents the directed graph of jobs in the workflow.
preprint.pdf is the preprint of the paper, also available on HAL.
 
The workflow code goes into a subfolder `workflow`, while the configuration is stored in a subfolder `config`. 
Inside of the `workflow` subfolder, the central `Snakefile` marks the entrypoint of the workflow (it will be automatically discovered when running snakemake from the root of above structure). 
In addition to the central `Snakefile`, rules are stored in a modular way, using the subfolder `workflow/rules`. 
Such modules should end with `.smk`, the recommended file extension of Snakemake. 
Further, scripts are stored in a subfolder `workflow/scripts`. 
Conda environments are stored in the subfolder `workflow/envs` (they are kept as finegrained as possible to improve transparency and maintainability).
The `config` subfolder contains a `requirement.txt` text file that lists the required packages, and a `config.yaml` file that contains all parameters needed for the workflow.

All output files generated in the workflow are stored under `results`, unless they are rather retrieved resources, in which case they should be stored under `resources`. 
The latter subfolder also contains small resources that shall be delivered along with the workflow via git, sorted into inhouse data, public data, and- if it is the case- private data.
The `results` subfolder separate `global` output that consists into resources that have been sorted, cleaned, and possibly arranged. The `results/network_analysis` sort workflow output in input data, intermediary data (typically intermediary object that yield final outputs), output, and plots. 

## How to run the workflow

This work relies on [Snakemake workflow management system](https://snakemake.readthedocs.io/en/stable/) and requires an API key to download in bulk the UN Comtrade data.

However, since raw UN Comtrade data are available in `resources` folder, you can render a random API key as an environment parameter and all subsequent jobs can be run one by one using snakemake commands.

To run a particular rule, you need to have [installed Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), then run in your terminal the following command:

```zsh
snakemake rule_name --sdm conda
```

<!-- Further information can be asked to the authors: [valentin.mathieu@agroparistech.fr](mailto:valentin.mathieu@agroparistech.fr) -->