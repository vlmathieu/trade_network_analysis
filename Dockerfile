FROM condaforge/miniforge3:latest

# --- 1. Snakemake itself, pinned to the version that produced your results ---
RUN conda install -n base -y -c conda-forge -c bioconda snakemake=8.25.5 \
 && conda clean --all -y

WORKDIR /workspace

# --- 2. Everything Snakemake needs to READ the workflow ---
COPY workflow/Snakefile  /workspace/workflow/Snakefile
COPY workflow/rules/     /workspace/workflow/rules/
COPY workflow/envs/      /workspace/workflow/envs/
COPY config/             /workspace/config/
COPY resources/          /workspace/resources/

# --- 3. Satisfy the Snakefile's `envvars:` declaration ---
ENV comtrade_apikey=unused-by-default

# --- 4. Pre-build all 7 conda environments INTO the image ---
RUN snakemake --sdm conda --conda-prefix /conda-envs \
      --conda-create-envs-only --cores 1 \
 && conda clean --all -y

# --- 5. The analysis scripts (copied last, on purpose) ---
COPY workflow/scripts/   /workspace/workflow/scripts/

CMD ["snakemake", "--sdm", "conda", "--conda-prefix", "/conda-envs", "--cores", "4"]
