FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="c49c4d20d6b206410ecbb95a1168ceabbec146d2bc538a47cf1c35c10b1f11a7"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/comtradeapicall.yaml
#   prefix: /conda-envs/b6f2f106f77298b30003ec2aa3c923e5
#   channels:
#     - conda-forge
#     - defaults
#     - bioconda
#   dependencies:
#     - polars =1.33.1
#     - pandas ==2.2.1
#     - urllib3 ==2.2.1
#     - tqdm =4.67.1
#     - pip:
#       - comtradeapicall ==1.2.1
RUN mkdir -p /conda-envs/b6f2f106f77298b30003ec2aa3c923e5
COPY workflow/envs/comtradeapicall.yaml /conda-envs/b6f2f106f77298b30003ec2aa3c923e5/environment.yaml

# Conda environment:
#   source: workflow/envs/network_connectivity.yaml
#   prefix: /conda-envs/a0c85d73c54d2364b071250667958461
#   channel:
#     - conda-forge
#   dependencies:
#     - polars =1.33.1
#     - numpy =2.3.3
#     - networkx =3.5
#     - scipy =1.16.2
RUN mkdir -p /conda-envs/a0c85d73c54d2364b071250667958461
COPY workflow/envs/network_connectivity.yaml /conda-envs/a0c85d73c54d2364b071250667958461/environment.yaml

# Conda environment:
#   source: workflow/envs/network_metrics.yaml
#   prefix: /conda-envs/4d9eef157aa80593cafacea682da88da
#   channel:
#     - conda-forge
#   dependencies:
#     - polars =1.33.1
#     - numpy =2.3.3
#     - networkx =3.5
RUN mkdir -p /conda-envs/4d9eef157aa80593cafacea682da88da
COPY workflow/envs/network_metrics.yaml /conda-envs/4d9eef157aa80593cafacea682da88da/environment.yaml

# Conda environment:
#   source: workflow/envs/polars.yaml
#   prefix: /conda-envs/2b2e512950e2cbade72b71bfc52c0a53
#   channels:
#     - conda-forge
#   dependencies:
#     - polars =1.33.1
RUN mkdir -p /conda-envs/2b2e512950e2cbade72b71bfc52c0a53
COPY workflow/envs/polars.yaml /conda-envs/2b2e512950e2cbade72b71bfc52c0a53/environment.yaml

# Conda environment:
#   source: workflow/envs/r_plots.yaml
#   prefix: /conda-envs/a27bef1e6a0b7ef398b437474553e66b
#   channels:
#    - conda-forge
#   dependencies:
#    - r-base=4.3.3
#    - r-tidyverse=2.0.0
#    - r-dplyr=1.1.4
#    - r-ggplot2=3.5.1
#    - r-svglite=2.1.3
#    - r-ggrepel=0.9.6
#    - r-patchwork=1.3.0
#    - r-ggh4x=0.3.0
#    - r-ggrepel=0.9.6
#    - r-ggpubr=0.6.0
#    - r-hrbrthemes=0.8.7
#    - r-patchwork=1.3.*
#    - r-maps=3.4.2.1
#    - r-geosphere=1.5
#   #  - r-CoordinateCleaner=3.0.1
#    - r-reshape2=1.4.4
RUN mkdir -p /conda-envs/a27bef1e6a0b7ef398b437474553e66b
COPY workflow/envs/r_plots.yaml /conda-envs/a27bef1e6a0b7ef398b437474553e66b/environment.yaml

# Conda environment:
#   source: workflow/envs/wbgapi.yaml
#   prefix: /conda-envs/9660e0bdab442ec46009dabfb59b68ea
#   channel:
#     - conda-forge
#   dependencies:
#     - polars =1.33.1
#     - pandas =2.2.1
#     - wbgapi =1.0.12
RUN mkdir -p /conda-envs/9660e0bdab442ec46009dabfb59b68ea
COPY workflow/envs/wbgapi.yaml /conda-envs/9660e0bdab442ec46009dabfb59b68ea/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/b6f2f106f77298b30003ec2aa3c923e5 --file /conda-envs/b6f2f106f77298b30003ec2aa3c923e5/environment.yaml && \
    conda env create --prefix /conda-envs/a0c85d73c54d2364b071250667958461 --file /conda-envs/a0c85d73c54d2364b071250667958461/environment.yaml && \
    conda env create --prefix /conda-envs/4d9eef157aa80593cafacea682da88da --file /conda-envs/4d9eef157aa80593cafacea682da88da/environment.yaml && \
    conda env create --prefix /conda-envs/2b2e512950e2cbade72b71bfc52c0a53 --file /conda-envs/2b2e512950e2cbade72b71bfc52c0a53/environment.yaml && \
    conda env create --prefix /conda-envs/a27bef1e6a0b7ef398b437474553e66b --file /conda-envs/a27bef1e6a0b7ef398b437474553e66b/environment.yaml && \
    conda env create --prefix /conda-envs/9660e0bdab442ec46009dabfb59b68ea --file /conda-envs/9660e0bdab442ec46009dabfb59b68ea/environment.yaml && \
    conda clean --all -y
