# Shared plotting constants and helpers.
# Sourced by the plot_*.R scripts via:
#   source(file.path(snakemake@scriptdir, "utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(hrbrthemes)
})

# --- Palette (see CLAUDE.md "Shared visualization conventions") ---
PAL_EXP      <- "#C32E5A"   # main exporters / exports
PAL_IMP      <- "#3D85F7"   # main importers / imports
PAL_BALANCED <- "#A49393"
PAL_TOTAL    <- "black"

PAL_DIRECTION <- c(imp = PAL_IMP, exp = PAL_EXP)

# Trader-type palette (keyed on the classification labels). Named, so plots
# match colours by key regardless of the order they list the levels in.
PAL_TRADER <- c(main_exp = PAL_EXP, balanced = PAL_BALANCED, main_imp = PAL_IMP)

# --- Benchmark gridlines ---
BENCHMARK_YEARS <- c(2000, 2005, 2010, 2015, 2020)

benchmark_lines <- function(years = BENCHMARK_YEARS, ...) {
  geom_vline(xintercept = years, linetype = "dotted",
             colour = "grey70", linewidth = 0.3, ...)
}

# --- Shared theme ---
theme_trade <- function(base_size = 9, axis_title_size = 10, ...) {
  theme_ipsum(base_size = base_size, axis_title_size = axis_title_size, ...)
}

# --- Aggregation level from an input path ---
# 'results/network_analysis/agg_eu/output/x.csv' -> 'agg_eu'
agg_lvl_from_path <- function(path) basename(dirname(dirname(path)))

# --- Save wrapper enforcing project conventions ---
save_figure <- function(plot, path, width, height, dpi = 300) {
  ggsave(filename = path, plot = plot, width = width, height = height,
         dpi = dpi, bg = "white", create.dir = TRUE)
}
