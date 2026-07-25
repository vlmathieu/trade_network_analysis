# Containerising this workflow with Docker

A step-by-step guide to freezing the `trade_network_analysis` Snakemake workflow into a
Docker image, so that the exact software stack that produced the published figures can be
re-run by anyone, on any machine, indefinitely.

Written for someone new to Docker. Each step explains both what it does *technically* and
what it means *in plain terms*.

The `Dockerfile` and `.dockerignore` in this repo are the finished result of this guide, and
were **built and run end-to-end on 2026-07-25** against the current code (24/24 jobs green,
66 PNG + 66 SVG + 13 CSV produced; the five `agg_eu` metric CSVs reproduce the tracked results
to a worst relative difference of ~1e-15). You do not have to edit anything to reproduce — the
steps below explain what is already in place.

---

## Why we are doing this

A conda environment pins your **packages**. It does not pin the operating system underneath
them — glibc, fontconfig, BLAS, the installed fonts. That matters here more than usual:

- This machine is **arm64** (Apple Silicon); conda-forge ships different binaries for
  `osx-arm64` and `linux-64`.
- The plotting code is sensitive to **text metrics** — `fit_bar_labels()` filters segments
  by a ~15 px/char width heuristic, and `compute_label_positions()` does collision
  avoidance on label extents. Different fonts on a reviewer's Linux box can change which
  labels appear.

A container pins the layer conda cannot reach.

```
┌─ CPU architecture              arm64 vs x86_64          ┐
├─ OS + system libraries         glibc, fontconfig, BLAS  ┤ ← the CONTAINER pins these
├─ packages                      r-base 4.3.3, ggplot2 3.5.2  ← conda pins these
└─ your code                     workflow/scripts/*           ← git pins these
```

---

## Docker in 60 seconds

| Word | Technically | In plain terms |
|---|---|---|
| **Image** | An immutable, layered filesystem identified by a hash | A frozen snapshot of a whole Linux computer — OS, R, Python, your code. Like a `.iso`, but composable. |
| **Container** | A running process tree isolated onto an image | A disposable machine booted from that snapshot. Delete it; the image is untouched. |
| **Build context** | The directory tree sent to the Docker engine before building | The pile of files you hand to Docker. Keeping it small matters — see Step 1. |

The model: **`Dockerfile` is a recipe → `docker build` cooks it into an image →
`docker run` serves a container from it.**

---

## Step 0 — Prerequisites (already satisfied in the code)

Two things must be true for the workflow to build a DAG from a **clean checkout** — i.e.
inside a fresh container, with nothing but `workflow/` + `config/` + `resources/`. Both are
already handled in the current code; this section is here so you understand *why* they
matter, not because you need to change anything.

**1. `rule all` points at paths a rule can actually produce.** The target list in
`workflow/Snakefile` uses the current `results/network_analysis/{agg_lvl}/plot/{fao_div}/…`
layout. (An earlier version pointed at a stale pre-`{agg_lvl}` path with an extra `{wgt}`
directory level that no rule produced — from a clean checkout that made the DAG unbuildable.
Fixed.)

**2. The API download is not re-triggered.** In a fresh container `correspondence_FAO_HS`
would regenerate its JSON, which — if that JSON were a live input to `get_uncomtrade` — would
make Snakemake re-download from the UN Comtrade API (needing a real key, and letting
back-data drift break reproducibility). `get_uncomtrade` wraps that input in `ancient(...)`,
so its timestamp is ignored and the pre-seeded `resources/public/uncomtrade_data.parquet.gzip`
is used as-is. A container dry-run confirms `get_uncomtrade`, `correspondence_FAO_HS` and
`wb_data` are **not** scheduled — the run starts at `deflate_uncomtrade`.

---

## Step 1 — `.dockerignore`

The file at the project root:

```
# Ignore everything by default
*

# Allow back only what the image actually needs
!workflow/
!config/
!resources/

# ...and never these, even inside the allowed folders
workflow/_sandbox/
**/__pycache__/
**/.DS_Store
```

**Technically:** `.dockerignore` filters the build context before it is transmitted to the
Docker daemon. The deny-all-then-allow-list pattern (`*` followed by `!`) is safer than
listing exclusions, because anything added to the project later is excluded by default
rather than silently shipped.

**In plain terms:** you are telling Docker *"only look at these three folders."* Without it,
Docker sends the *entire* project — `.snakemake/` alone is tens of GB — before building
anything. With it, the context is ~41 MB, from unusable to instant.

It also protects you: `_writing/` is confidential and gitignored, and this guarantees it can
never be baked into an image you might publish.

---

## Step 2 — The Dockerfile

```dockerfile
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

# --- 1. Snakemake in a DEDICATED env (not base — see note below) ---
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
```

### What each block does

**Base — pinned by digest.**
*Technically:* `FROM …@sha256:70ce3950…` pins the exact base image content, not the floating
`:latest` tag. The digest is a **multi-arch manifest list**, so the same line resolves to the
arm64 layer on Apple Silicon and the amd64 layer when you build the portable artifact.
*In plain terms:* the OS layer conda cannot reach is now frozen too. Rebuilding a year from
now starts from the identical Linux, no matter where `:latest` has moved. **This pin is not
cosmetic** — see "Two failures the clean build caught" below.

**Block 0 — provenance labels.**
Standard OCI labels so the image is self-describing (`docker inspect`). Pass the source commit
at build time with `--build-arg VCS_REF=$(git rev-parse --short HEAD)`; it defaults to
`unknown` if omitted.

**Block 1 — Snakemake in a dedicated env.**
*Technically:* `conda create -n smk … snakemake=8.25.5`, then prepend `/opt/conda/envs/smk/bin`
to `PATH`. Snakemake is **not** installed into base.
*In plain terms:* the pinned base ships **python 3.13** in its base env, and snakemake 8.25.5
cannot solve against it (see below). A separate env (python 3.12) decouples the driver from
the base python. `/opt/conda/bin` stays later on `PATH`, so `conda`/`mamba` are still found
for `--sdm conda`. Pinning to **8.25.5** — the version that produced the published results —
makes the image immune to future Snakemake releases.

**Block 2 — copy what is needed to parse the workflow.**
*Technically:* Snakemake reads the Snakefile, `.smk` rule files, `config.yaml` and the env
YAMLs to build the DAG. `workflow/scripts/utils.R` is copied here too — it is a declared
`input:` of `plot_network_connectivity` and `plot_network_composition`, so the DAG build in
Block 4 fails without it (see below). It is copied separately from the rest of the scripts
because it changes rarely; editing any other script still skips the env re-solve.

**Block 3 — `ENV comtrade_apikey`.**
*Technically:* `workflow/Snakefile` has `envvars: 'comtrade_apikey'`. Snakemake raises a
`WorkflowError` at *parse* time if that variable is undefined — before any rule runs.
*In plain terms:* **without this line the build fails immediately**, even though the pre-seeded
data means the API is never called. A dummy value satisfies the check. Override it at run time
if you ever do want to re-download.

**Block 4 — pre-build the environments. This is the heart of it.**
*Technically:* `--conda-create-envs-only` makes Snakemake construct every rule's environment
without executing any job, writing them to `/conda-envs`.
*In plain terms:* all seven environments — R 4.3.3, ggplot2 3.5.2, polars, circlize,
everything — are solved *once, now*, and frozen into the image. Every future run uses these
exact binaries. No solving, no drift.

> **Why `--conda-prefix /conda-envs` and not the default?**
> The default is `.snakemake/conda` *inside the working directory*. In Step 4 a host folder is
> mounted over `/workspace/results`, and if the envs lived under `/workspace` a mount could
> hide them and silently force a rebuild. Putting them at `/conda-envs`, outside the working
> directory, keeps them untouchable. **This is the most important flag in this file.**

**Block 5 — scripts copied last.**
*Technically:* each Dockerfile instruction is a cached layer; changing one invalidates every
layer below it.
*In plain terms:* scripts change often, environments almost never. Copying scripts *after* env
creation means editing an R figure script rebuilds in seconds instead of re-solving all seven
environments. (The one exception is `utils.R`, promoted to Block 2 because the DAG build needs
it — editing `utils.R` therefore *does* re-solve, but it changes rarely.)

### Two failures the clean build caught

Building this from a clean checkout on 2026-07-25 surfaced two breaks that an incremental
build on the original author's machine never hit. Both are now fixed above, and both are
exactly the kind of latent problem a reproducibility container exists to expose:

1. **Base drift → python 3.13.** Installing snakemake 8.25.5 into base failed with a solver
   conflict, because today's `condaforge/miniforge3` base pins python 3.13. The original
   Dockerfile used floating `:latest` and would fail identically today. Fix: dedicated `smk`
   env (Block 1). Pinning the base by digest is what makes this reproducible going forward.
2. **`utils.R` input timing.** `--conda-create-envs-only` still builds the full DAG, which
   needs every declared input file present. The `utils.R` refactor made two plot rules declare
   `workflow/scripts/utils.R` as an `input:`, but the scripts were copied *last* → the DAG
   build hit `MissingInputException`. Fix: copy `utils.R` before Block 4 (Block 2).

> **Why not `snakemake --containerize`?**
> That tool generated the original Dockerfile. We hand-write instead because `--containerize`
> bakes in a dependency hash that must match what the host Snakemake computes later — a fragile
> coupling that fails *silently* (it rebuilds environments against today's dependencies instead
> of erroring). Letting Snakemake build its own environments, at its own paths, with its own
> hashing, removes that failure mode and picks up all 7 env files automatically.

---

## Step 3 — Build the image

```zsh
docker build --build-arg VCS_REF=$(git rev-parse --short HEAD) -t trade_network_analysis:1.0 .
```

**Technically:** reads `Dockerfile`, sends the (~41 MB) context, executes each instruction
into a cached layer, tags the result, and stamps the source commit into the image labels.

**In plain terms:** cooks the recipe into the frozen snapshot. **The first build takes
15–30 minutes**, mostly the R environments in Block 4. Later builds reuse cached layers and
are much faster.

Check the result:

```zsh
docker images trade_network_analysis          # ~6 GB is normal for two R environments
docker inspect trade_network_analysis:1.0 --format '{{json .Config.Labels}}'
```

---

## Step 4 — Run the workflow

**Dry run first**, to confirm the workflow parses, the environments are found, and no API
download is scheduled:

```zsh
docker run --rm trade_network_analysis:1.0 \
  snakemake --sdm conda --conda-prefix /conda-envs --cores 4 -n
```

A correct plan lists **24 jobs starting at `deflate_uncomtrade`**, with
`get_uncomtrade` / `correspondence_FAO_HS` / `wb_data` appearing only under "Rules with
missing metadata" (i.e. pre-seeded, not re-run).

**Then the real run:**

```zsh
mkdir -p ~/archival_run
docker run --rm \
  -v ~/archival_run:/workspace/results \
  trade_network_analysis:1.0
```

**Technically:** `-v host:container` bind-mounts a host directory over the container's
`results/`, so writes land on the host filesystem and survive. `--rm` deletes the container
afterwards.

**In plain terms:** *"run the analysis, and put the outputs in a folder I can actually see."*
Without `-v`, results vanish when the container exits. Only `results/` is mounted — not the
whole project — precisely so the mount cannot shadow anything else (including `/conda-envs`).

A green run reports `24 of 24 steps (100%) done` and writes 66 PNG + 66 SVG + 13 CSV under
`~/archival_run/network_analysis/{agg_eu,country_lvl}/`.

Other invocations:

```zsh
# just one rule
docker run --rm -v ~/archival_run:/workspace/results \
  trade_network_analysis:1.0 \
  snakemake plot_trade_flows --sdm conda --conda-prefix /conda-envs --cores 4

# poke around inside
docker run --rm -it trade_network_analysis:1.0 bash
```

That last one drops you into a shell inside the image. Look at `/conda-envs`, run
`R --version`, then `exit`. Nothing you do persists.

> **On exact reproducibility.** Metric CSVs reproduce the tracked results to ~15 significant
> figures (worst relative difference ~1e-15); the residual is ordinary floating-point
> nondeterminism across BLAS/arch builds, not a logic difference. Raster figures (`.png`) can
> differ at the pixel level across architectures because of font rendering — which is *why*
> the archived artifact is a specific image, verified on the architecture it ships for
> (Step 5), rather than a promise of bit-identical PNGs on every host.

---

## Step 5 — Produce the archival artifact

### Architecture

Steps 3–4 build a native **arm64** image, which x86_64 Linux users (most reviewers) cannot
run natively. Once the native build is confirmed working, rebuild for the portable
architecture:

```zsh
docker buildx build --platform linux/amd64 \
  --build-arg VCS_REF=$(git rev-parse --short HEAD) \
  -t trade_network_analysis:1.0-amd64 --load .
```

Do this *second*, not first — it runs under emulation on Apple Silicon and is slower. Iterate
on the fast native build; produce the artifact once. Confirm it before shipping:

```zsh
docker inspect trade_network_analysis:1.0-amd64 --format 'arch={{.Architecture}}'   # amd64
docker run --rm --platform linux/amd64 trade_network_analysis:1.0-amd64 \
  bash -lc 'ls /conda-envs/*.env_setup_done | wc -l'                                # 7
```

### Distribution

⚠️ A `sha256:` *registry* digest only exists once an image is pushed to a registry; local
images do not have one. Two options:

**Option A — push to a registry.** GitHub Container Registry is free for public images. Cite
`ghcr.io/<user>/trade_network_analysis@sha256:…` in the data-availability statement. (Hold
this until publication — see SD-4.)

**Option B — export a file.** Better suited to recherche.data.gouv.fr, which hosts files
rather than images:

```zsh
docker save trade_network_analysis:1.0-amd64 | gzip > trade_network_analysis-1.0-amd64.tar.gz
shasum -a 256 trade_network_analysis-1.0-amd64.tar.gz
```

Deposit the tarball beside the data and cite the checksum. Anyone can then
`docker load < trade_network_analysis-1.0-amd64.tar.gz` and reproduce the figures exactly.

> **Cite the image ID for *identity*, the tarball checksum for *download integrity*.**
> `docker save | gzip` embeds a timestamp, so re-saving the *same* image yields a *different*
> tarball checksum — the checksum verifies one specific deposited file, not the image. The
> content-addressed image ID (`docker images --no-trunc --format '{{.ID}}'`) is the stable
> identity to cite in prose.

**Reference build (2026-07-25, `amd64`, commit `77dec67`):**

| | value |
|---|---|
| Image ID | `sha256:bc2ce83046c760ada4bac6ae966cc1023549c74621c2c5de54859f1f9a2492ca` |
| Tarball sha256 | `e7e79d2cdfa0634a54a464f95b0ce5e8fe130916e93703ca6d489fff717b7f7d` |

Regenerate these from the final committed state at publication; the values above pin *this*
build and will change with any code, base-image, or dependency change.

### Data-availability statement (template)

> The complete analysis workflow is packaged as a Docker image
> (`trade_network_analysis:1.0`, linux/amd64, image ID `sha256:bc2ce830…2492ca`), archived at
> [recherche.data.gouv.fr DOI] together with the pre-seeded input data. Loading the image
> (`docker load < trade_network_analysis-1.0-amd64.tar.gz`; sha256 `e7e79d2c…7b7f7d`) and
> running `docker run --rm -v <out>:/workspace/results trade_network_analysis:1.0-amd64`
> regenerates every figure and metric table without a UN Comtrade API key.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Build hangs at "sending build context" | `.dockerignore` missing or wrong | Step 1 — verify it is at the project root |
| `conda install … snakemake` solver conflict on base | base image ships a python snakemake can't solve against | Install snakemake in a dedicated env (Block 1) |
| `MissingInputException … utils.R` at env creation | a declared script `input:` is copied after Block 4 | Copy it before Block 4 (Block 2) |
| `environment variable comtrade_apikey is requested` | the `envvars:` directive | the `ENV` line in Block 3 |
| Environments rebuild on every run (slow) | something shadowed `/conda-envs` | check `--conda-prefix` is on *both* build and run commands |
| `results/` empty after a run | no `-v` mount | Step 4 |
| Editing a script triggers a 20-min rebuild | a `COPY scripts` line too high in the Dockerfile | keep Block 5 last (only `utils.R` is promoted) |
| `no space left on device` | image layers accumulating | `docker system prune` |

---

## One honest caveat

The base is pinned by **digest**, so the OS layer no longer floats — but the digest is still
attached to the `:latest` *tag*, which will eventually move. That is fine: the artifact you
archive is the **built image**, identified by its image ID or the tarball checksum, and it is
immutable regardless of what `:latest` points to later.

The Dockerfile is the recipe; the image is the specimen. **Archive the specimen.**

---

## Appendix — Apptainer, and why Docker here

Snakemake has no `--sdm docker`; its per-rule container support is Apptainer/Singularity only.
What this guide does instead is run the **whole workflow inside one Docker container**, which
needs no Snakemake container support at all.

This is not a lock-in decision: Docker and Apptainer consume the same images
(`apptainer pull docker://…`). Apptainer becomes the right choice if this workflow ever moves
to an HPC cluster, where the root-owned Docker daemon is usually forbidden.

Note also that Docker Desktop on macOS runs its own hidden Linux VM. Using Docker does not
eliminate the VM — it just means you no longer have to manage one yourself.
