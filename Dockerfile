# check=skip=SecretsUsedInArgOrEnv
# ^ comtrade_apikey ENV below is a public dummy to satisfy Snakemake's envvars: parse check,
#   not a secret. Pre-seeded data means the API is never called. Override at run time to re-download.

# Base pinned by manifest-list digest (multi-arch: arm64 native build, amd64 artifact).
# Rebuilds start from an identical OS layer regardless of what :latest points to later.
FROM condaforge/miniforge3:latest@sha256:70ce395025aeec0d011eb64b6d6dc870fcc1a774a163965bc6dd8e8f6bf5da0d

# --- 0. Provenance metadata (VCS_REF passed at build: --build-arg VCS_REF=$(git rev-parse --short HEAD)) ---
ARG VERSION=1.0
ARG VCS_REF=unknown
LABEL org.opencontainers.image.title="trade_network_analysis" \
      org.opencontainers.image.description="Reproducible Snakemake workflow for the timber trade network analysis. Regenerates all figures from pre-seeded data, no API key required." \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

# --- 1. Snakemake itself, pinned to the version that produced the results.
#        Installed in a DEDICATED env, not base: this pinned base image ships python=3.13
#        in base, which snakemake 8.25.5 cannot solve against. A separate env (python 3.12)
#        decouples the driver from the base python. PATH prepend makes it the default;
#        /opt/conda/bin stays on PATH so `conda`/`mamba` remain available for --sdm conda. ---
RUN conda create -n smk -y -c conda-forge -c bioconda snakemake=8.25.5 \
 && conda clean --all -y
ENV PATH=/opt/conda/envs/smk/bin:$PATH

WORKDIR /workspace

# --- 2. Everything Snakemake needs to READ the workflow ---
COPY workflow/Snakefile  /workspace/workflow/Snakefile
COPY workflow/rules/     /workspace/workflow/rules/
COPY workflow/envs/      /workspace/workflow/envs/
COPY config/             /workspace/config/
COPY resources/          /workspace/resources/
# utils.R is a declared input: of plot_network_connectivity + plot_network_composition,
# so the DAG build in Block 4 needs it present. Copied here (not with the rest of the
# scripts in Block 5) because it changes rarely; editing other scripts still skips env re-solve.
COPY workflow/scripts/utils.R /workspace/workflow/scripts/utils.R

# --- 3. Satisfy the Snakefile's `envvars:` declaration ---
ENV comtrade_apikey=unused-by-default

# --- 4. Pre-build all 7 conda environments INTO the image ---
RUN snakemake --sdm conda --conda-prefix /conda-envs \
      --conda-create-envs-only --cores 1 \
 && conda clean --all -y

# --- 5. The analysis scripts (copied last, on purpose) ---
COPY workflow/scripts/   /workspace/workflow/scripts/

CMD ["snakemake", "--sdm", "conda", "--conda-prefix", "/conda-envs", "--cores", "4"]
