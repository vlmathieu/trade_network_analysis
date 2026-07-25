# Reproducing the analysis with Docker

This repository ships a Docker recipe that regenerates **every figure and table in the
paper** from the pre-seeded data — with no UN Comtrade API key and no manual setup. One
command boots a container that carries the entire software stack (operating system, R,
Python, all packages, and the analysis code) and runs the whole workflow.

Written for someone new to Docker. If you have never used Docker before, read "Docker in
60 seconds" first, then follow "Quick start".

---

## Why a container (30 seconds)

A conda environment pins your **packages**, but not the operating system underneath them —
`glibc`, `fontconfig`, `BLAS`, the installed fonts. That matters here: the plotting code is
sensitive to **text metrics** (it filters and places chart labels by their measured pixel
width), so different system fonts can change which labels appear. A container freezes the
layer conda cannot reach.

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
| **Image** | An immutable, layered filesystem identified by a hash | A frozen snapshot of a whole Linux computer — OS, R, Python, the code. Like a `.iso`, but composable. |
| **Container** | A running process tree isolated onto an image | A disposable machine booted from that snapshot. Delete it; the image is untouched. |
| **Volume (`-v`)** | A host directory bind-mounted into the container | A shared folder, so results the container writes land on your real disk. |

The model: **an image is a frozen computer → `docker run` boots a throwaway machine from it →
it does the work → it disappears, leaving only the files you asked it to save.**

---

## Before you start

Install Docker:

- **macOS / Windows:** [Docker Desktop](https://www.docker.com/products/docker-desktop/).
- **Linux:** [Docker Engine](https://docs.docker.com/engine/install/).

Confirm it works:

```zsh
docker --version
```

You do **not** need R, Python, conda, or an API key on your machine — they all live inside
the image.

---

## Quick start — reproduce the figures

Two ways to get the image. **Pick one.**

### Option A — load the published image (fastest)

The data deposit includes a file named like `trade_network_analysis-1.0-amd64.tar.gz`.

```zsh
# 1. (recommended) check the download is intact — compare to the published checksum
shasum -a 256 trade_network_analysis-1.0-amd64.tar.gz

# 2. load it into Docker
docker load < trade_network_analysis-1.0-amd64.tar.gz
```

This image is built for **Intel/AMD (`amd64`)** chips. It runs natively on Intel/AMD Linux,
on Windows, and on Intel Macs. On an **Apple Silicon Mac** it runs only under slow emulation
and one library may crash — use **Option B** there instead.

### Option B — build the image from source (any machine)

From the root of this repository:

```zsh
docker build -t trade_network_analysis:1.0 .
```

The first build takes **15–30 minutes** (it solves the R environments once, then freezes
them in). It works on any chip, because it builds the version native to your machine — this
is the right choice on an Apple Silicon Mac.

### Run it

```zsh
mkdir -p ~/repro_out
docker run --rm -v ~/repro_out:/workspace/results trade_network_analysis:1.0
```

- `--rm` throws the container away when it finishes (the image stays).
- `-v ~/repro_out:/workspace/results` is the shared folder: everything the workflow writes
  appears in `~/repro_out` on your disk. **Without it, the results vanish** when the
  container exits.

When it's done you'll see `24 of 24 steps (100%) done`, and your outputs are under
`~/repro_out/network_analysis/` — **66 PNG + 66 SVG figures and 13 CSV tables**, split into
`agg_eu/` (the EU-as-one-node analysis in the paper) and `country_lvl/` (the supplementary
country-level analysis).

**Preview without running** (lists the jobs, changes nothing):

```zsh
docker run --rm trade_network_analysis:1.0 \
  snakemake --sdm conda --conda-prefix /conda-envs --cores 4 -n
```

**Run a single figure** instead of everything, e.g. the trade-flow diagrams:

```zsh
docker run --rm -v ~/repro_out:/workspace/results trade_network_analysis:1.0 \
  snakemake plot_trade_flows --sdm conda --conda-prefix /conda-envs --cores 4
```

---

## What "reproduce" means here

- **Tables (CSV):** numbers match the published results to about **15 significant figures**.
  The tiny differences beyond that are ordinary floating-point rounding that varies between
  CPUs — not an error.
- **Figures (PNG):** may differ by a pixel or two across different chips, because of font
  rendering. That is why the archived, citeable object is a **specific image verified on the
  architecture it ships for**, not a promise of pixel-identical PNGs on every computer.

---

## Verifying what you ran

Every image is self-describing:

```zsh
docker inspect trade_network_analysis:1.0 --format '{{json .Config.Labels}}'
```

The published reference build (native `amd64`) has these fixed identifiers:

| | value |
|---|---|
| Image ID | `sha256:392e8841e6a1202e9daee08728c9815cd5185c60130ebf105d8760fd3e99b76a` |
| Tarball sha256 | `cfac8a1a7b9eb1682898d2b86a5c4626920d785d1daa415c729a2a64f26716d4` |

The **image ID** is the stable identity of the image itself. The **tarball checksum** verifies
that the file you downloaded is intact (re-saving the same image produces a byte-different
tarball, so the checksum pins the deposited file, not the image).

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `~/repro_out` is empty after a run | no `-v` shared folder | add `-v ~/repro_out:/workspace/results` (see "Run it") |
| `exec format error`, or extremely slow on Apple Silicon | you loaded the Intel/AMD image on an ARM Mac | build from source instead (Option B) |
| `no space left on device` | old images/layers piling up | `docker system prune` |
| build/run seems to hang at the start | first-time image download or env solve | wait — the first build solves the R environments once (15–30 min) |

---

## Citing this artifact (data-availability statement)

> The complete analysis workflow is packaged as a Docker image (`trade_network_analysis:1.0`,
> `linux/amd64`, image ID `sha256:392e8841…9b76a`), archived at [DOI] together with the
> pre-seeded input data. Loading the image
> (`docker load < trade_network_analysis-1.0-amd64.tar.gz`; sha256 `cfac8a1a…16716d4`) and
> running `docker run --rm -v <out>:/workspace/results trade_network_analysis:1.0` regenerates
> every figure and metric table without a UN Comtrade API key.

---

## Advanced notes

**What's inside the image.** The `Dockerfile` at the repository root is the full recipe:
it starts from a digest-pinned `condaforge/miniforge3` base, installs Snakemake 8.25.5, and
bakes in the seven conda environments the workflow needs (R 4.3.3, ggplot2 3.5.2, polars,
circlize, and the rest). The pre-seeded data in `resources/` lets every downstream rule run
offline. Read the Dockerfile's inline comments if you want to adapt or rebuild it.

**Running on an HPC cluster (Apptainer/Singularity).** Docker and Apptainer consume the same
images, so you can `apptainer pull docker://…` (or convert the tarball) and run there without
a Docker daemon — the usual situation on clusters, where Docker is not permitted.

**Longevity.** The base image is pinned by digest, so rebuilds start from an identical OS
layer. But the durable, citeable object is the **built image** (by image ID or tarball
checksum), which is immutable regardless of what upstream tags do later. The Dockerfile is the
recipe; the image is the specimen — archive the specimen.
