# The Rise (and Fall?) of China's Demand for Imported Timber

**A network analysis of the international timber trade.**
This repository holds the code, input data, and reproducible workflow behind the paper
by Valentin Mathieu and David W. Shanafelt (Forest Policy and Economics, under review).

## About

We model the global timber trade as a network of directed flows between countries and
track how its structure evolved from **1996 to 2023**. The analysis covers three
wood-product groups — **roundwood, sawnwood, and wood-based panels** (FAO divisions 01,
05, 07) — built from UN Comtrade bilateral trade records.

The headline finding: the network grew steadily more interconnected and concentrated,
polarising around China's import demand (over half of total trade value by 2021) while a
major exporter, Russia, structurally withdrew — creating, then partly unwinding, a single
point of vulnerability in the global supply chain.

The full study, abstract, and results are in the paper (DOI added on publication).

## Repository layout

```
trade_network_analysis/
├── workflow/            # the Snakemake pipeline (Snakefile, rules, scripts, conda envs)
├── config/             # config.yaml (all parameters) + pinned requirements
├── resources/          # pre-seeded input data (UN Comtrade, World Bank, FAO↔HS table)
├── Dockerfile          # builds a container with the whole software stack baked in
├── docker-tutorial.md  # step-by-step reproduction with Docker (start here to reproduce)
├── dag.pdf             # the workflow's job dependency graph
└── LICENSE             # MIT
```

Every tracked folder carries its own `README.md`. Start with
[`workflow/README.md`](workflow/README.md) for the pipeline,
[`config/README.md`](config/README.md) for parameters, and
[`resources/README.md`](resources/README.md) for the input data.

> The pipeline writes its outputs to `results/` (figures, metric tables). That folder is
> **not** stored in git — it is regenerated when you run the workflow.

## Reproduce the analysis

The pre-seeded data in `resources/` lets the entire pipeline run **without a UN Comtrade
API key**. There are two ways to run it.

### Option A — Docker (recommended, no setup)

One command regenerates every figure and table in the paper inside a container that
carries the exact operating system, R, Python, and package versions used. No R, Python,
or conda needed on your machine.

```zsh
docker build -t trade_network_analysis:1.0 .
mkdir -p ~/repro_out
docker run --rm -v ~/repro_out:/workspace/results trade_network_analysis:1.0
```

Full walkthrough, a pre-built image option, and troubleshooting are in
[`docker-tutorial.md`](docker-tutorial.md).

### Option B — Snakemake + conda

If you already use [Snakemake](https://snakemake.readthedocs.io/), run from the project root:

```zsh
snakemake --sdm conda --cores 4                 # full pipeline (builds conda envs on first run)
snakemake <rule> --sdm conda --cores 4          # a single rule
snakemake --sdm conda --cores 4 -n              # dry-run: list jobs, run nothing
```

There is no test suite: verify a rule by running it and inspecting `results/` and
`workflow/logs/<rule>.log`. Only the `get_uncomtrade` rule needs the `comtrade_apikey`
environment variable, and only if you want to re-download the raw trade data — every other
rule runs from the pre-seeded data. See [`workflow/README.md`](workflow/README.md) for the
pipeline structure.

## Data and software availability

- **Code** — this repository (MIT-licensed, see [`LICENSE`](LICENSE)); reuse and adapt freely.
- **Input data** — third-party records from [UN Comtrade](https://comtrade.un.org/) and the
  [World Bank](https://data.worldbank.org/), plus an in-house FAO↔HS product correspondence
  table, all pre-seeded under [`resources/`](resources/).
- **Archived deposit** — on publication, the input data and a citeable Docker image are
  deposited on [recherche.data.gouv.fr](https://recherche.data.gouv.fr/) with a DOI (added
  here at that point). The Docker image tarball is deliberately not stored in git.

## How to cite

Please cite the paper (citation and DOI added on publication) if you use this code or data.

## Contact

Valentin Mathieu — open an [issue](https://github.com/vlmathieu/trade_network_analysis/issues)
or reach the authors for questions about the workflow.
