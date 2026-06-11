library("ggplot2")
library("hrbrthemes")
suppressPackageStartupMessages(library("dplyr"))
library("ggrepel")
library("patchwork")
library("scales")
library("reshape2")
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("igraph"))

pal              <- c("#3D85F7", "#C32E5A")
trader_type_pal  <- c("main_exp" = "#C32E5A", "balanced" = "#A49393", "main_imp" = "#3D85F7")

# Same threshold logic as trajectory_plot() in plot_contributor_profiles.R
select_contributors <- function(profiles, prod, year, threshold) {
  d <- profiles %>% filter(cmd == prod, period == year) # nolint
  tot_exp <- sum(d$primary_value_exp, na.rm = TRUE)
  tot_imp <- sum(d$primary_value_imp, na.rm = TRUE)
  d %>%
    filter(primary_value_exp / tot_exp >= threshold | # nolint
           primary_value_imp / tot_imp >= threshold) %>%
    pull(country) # nolint
}

# Classify from exp_share column (already in contributor_profiles.csv)
get_trader_type <- function(profiles, prod, year, countries) {
  d <- profiles %>%
    filter(cmd == prod, period == year, country %in% countries) # nolint
  types <- ifelse(d$exp_share >= 0.8, "main_exp",
                  ifelse(d$exp_share <= 0.2, "main_imp", "balanced"))
  setNames(types, d$country)
}

# Build square flow matrix from long-format bilateral flows
build_flow_matrix <- function(d_flows, value_col, countries) {
  mat <- matrix(0, nrow = length(countries), ncol = length(countries),
                dimnames = list(countries, countries))
  d_agg <- d_flows %>%
    group_by(exporter_desc, importer_desc) %>% # nolint
    summarise(val = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  for (k in seq_len(nrow(d_agg))) {
    e <- d_agg$exporter_desc[k]
    m <- d_agg$importer_desc[k]
    if (e %in% countries && m %in% countries) mat[e, m] <- d_agg$val[k]
  }
  mat
}

# Two side-by-side chord diagrams: (a) exporter reports / (b) importer reports
chord_plot <- function(flows, profiles, prod, year, threshold, output_path, ext) {
  contributors <- select_contributors(profiles, prod, year, threshold)
  if (length(contributors) < 2) return(invisible(NULL))

  d <- flows %>%
    filter(fao_code_agg == prod, period == year, # nolint
           exporter_desc %in% contributors,
           importer_desc %in% contributors) %>%
    mutate(
      pv_exp = ifelse(is.na(primary_value_exp), 0, primary_value_exp) / 1e6,
      pv_imp = ifelse(is.na(primary_value_imp), 0, primary_value_imp) / 1e6
    )

  if (nrow(d) == 0) return(invisible(NULL))

  # Country order: descending total exp flow (applied to both panels for comparability)
  flow_totals <- bind_rows(
    d %>% rename(country = exporter_desc) %>% group_by(country) %>% # nolint
      summarise(val = sum(pv_exp, na.rm = TRUE), .groups = "drop"),
    d %>% rename(country = importer_desc) %>% group_by(country) %>% # nolint
      summarise(val = sum(pv_exp, na.rm = TRUE), .groups = "drop")
  ) %>%
    group_by(country) %>% summarise(total = sum(val), .groups = "drop") # nolint
  country_order <- flow_totals$country[order(flow_totals$total, decreasing = TRUE)]
  # Add any contributor with no flows in this year at the end
  missing       <- setdiff(contributors, country_order)
  country_order <- c(country_order, missing)

  mat_exp <- build_flow_matrix(d, "pv_exp", country_order)
  mat_imp <- build_flow_matrix(d, "pv_imp", country_order)

  types    <- get_trader_type(profiles, prod, year, country_order)
  arc_cols <- trader_type_pal[types[country_order]]
  arc_cols[is.na(arc_cols)] <- "#A49393"
  names(arc_cols) <- country_order

  country_labels <- function() {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      circos.text(mean(xlim), 0.5, get.cell.meta.data("sector.index"),
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0.5), cex = 0.7)
    }, bg.border = NA)
  }

  if (ext == "png") {
    png(output_path, width = 2 * 2480, height = 2480, res = 300, bg = "white")
  } else {
    svg(output_path, width = 2 * 2480 / 300, height = 2480 / 300, bg = "white")
  }

  par(mfrow = c(1, 2))

  circos.clear()
  chordDiagram(mat_exp, grid.col = arc_cols, transparency = 0.4,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.1))
  country_labels()
  title("(a) Exporter reports (million USD)", cex.main = 1.2, font.main = 2)
  circos.clear()

  chordDiagram(mat_imp, grid.col = arc_cols, transparency = 0.4,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.1))
  country_labels()
  title("(b) Importer reports (million USD)", cex.main = 1.2, font.main = 2)
  circos.clear()

  dev.off()
}

# One ggplot2 network panel: nested bubbles (same encoding as contributor_profiles)
# with straight-line edges underneath
build_network_panel <- function(flows, profiles, prod, year, contributors,
                                layout_df, size_max, panel_letter) {
  d_flows <- flows %>%
    filter(fao_code_agg == prod, period == year, # nolint
           exporter_desc %in% contributors,
           importer_desc %in% contributors) %>%
    mutate(pv_exp = ifelse(is.na(primary_value_exp), 0, primary_value_exp) / 1e6) %>%
    filter(pv_exp > 0) %>%
    left_join(layout_df %>% rename(exporter_desc = country, x0 = layout_x, y0 = layout_y),
              by = "exporter_desc") %>%
    left_join(layout_df %>% rename(importer_desc = country, x1 = layout_x, y1 = layout_y),
              by = "importer_desc") %>%
    filter(!is.na(x0), !is.na(x1))

  d_nodes <- profiles %>%
    filter(cmd == prod, period == year, country %in% contributors) %>% # nolint
    mutate(size_exp = primary_value_exp / 1e6,
           size_imp = primary_value_imp / 1e6)

  # Long format: two rows per country, larger circle rendered first
  d_long <- bind_rows(
    d_nodes %>% mutate(trade_type = "Total exported value", size_val = size_exp),
    d_nodes %>% mutate(trade_type = "Total imported value", size_val = size_imp)
  ) %>%
    arrange(country, desc(size_val)) %>% # nolint
    left_join(layout_df, by = "country")

  d_labels <- d_nodes %>% left_join(layout_df, by = "country")

  ggplot() +
    geom_segment(data = d_flows,
                 aes(x = x0, y = y0, xend = x1, yend = y1, linewidth = pv_exp), # nolint
                 alpha = 0.15, color = "grey50",
                 show.legend = FALSE) +
    geom_point(data = d_long,
               aes(x = layout_x, y = layout_y, size = size_val, color = trade_type), # nolint
               alpha = 0.6) +
    geom_text_repel(data = d_labels,
                    aes(x = layout_x, y = layout_y, label = country), # nolint
                    size = 3, color = "black",
                    bg.color = "white", bg.r = 0.2,
                    segment.colour = "grey55", segment.size = 0.25,
                    min.segment.length = 0,
                    box.padding = 0.6, point.padding = 0.8,
                    force = 4, seed = 42,
                    direction = "both", max.overlaps = Inf,
                    show.legend = FALSE) +
    scale_color_manual(
      values = c("Total exported value" = pal[2], "Total imported value" = pal[1]),
      name   = NULL
    ) +
    scale_size_continuous(range  = c(1, 20),
                          limits = c(0, size_max),
                          labels = scales::label_comma(),
                          name   = "Traded value (million US$)") +
    scale_linewidth_continuous(range = c(0.1, 2), guide = "none") +
    theme_void() +
    theme(legend.position  = "bottom",
          plot.title       = element_text(size = 10, face = "bold"),
          plot.background  = element_rect(fill = "white", color = NA),
          plot.margin      = margin(10, 20, 10, 20)) +
    labs(title = paste0(panel_letter, " ", year))
}

# Two-panel network figure (year_start left, year_end right)
network_plot <- function(flows, profiles, prod, year_start, year_end, threshold) {
  contributors <- select_contributors(profiles, prod, year_end, threshold)
  if (length(contributors) < 2) return(NULL)

  # Layout from year_end graph; reused for year_start (stable positions)
  d_edges_end <- flows %>%
    filter(fao_code_agg == prod, period == year_end, # nolint
           exporter_desc %in% contributors,
           importer_desc %in% contributors) %>%
    mutate(pv_exp = ifelse(is.na(primary_value_exp), 0, primary_value_exp) / 1e6) %>%
    filter(pv_exp > 0) %>%
    select(from = exporter_desc, to = importer_desc, weight = pv_exp)

  if (nrow(d_edges_end) == 0) return(NULL)

  g <- igraph::graph_from_data_frame(
    d_edges_end, directed = TRUE,
    vertices = data.frame(name = contributors)
  )
  set.seed(42)
  coords     <- igraph::layout_with_fr(g)
  layout_df  <- data.frame(country  = igraph::V(g)$name,
                            layout_x = coords[, 1],
                            layout_y = coords[, 2])

  # Shared size maximum across both years
  d_both   <- profiles %>%
    filter(cmd == prod, period %in% c(year_start, year_end), # nolint
           country %in% contributors)
  size_max <- ceiling(max(d_both$primary_value_exp, d_both$primary_value_imp,
                          na.rm = TRUE) / 1e6)

  p_start <- build_network_panel(flows, profiles, prod, year_start, contributors,
                                  layout_df, size_max, "(a)")
  p_end   <- build_network_panel(flows, profiles, prod, year_end,   contributors,
                                  layout_df, size_max, "(b)")

  (p_start + p_end) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
}

# ── Main loop ───────────────────────────────────────────────────────────────

threshold  <- snakemake@params$threshold
year_start <- snakemake@params$year_start
year_end   <- snakemake@params$year_end

n_levels      <- length(snakemake@input) / 2
flows_files   <- snakemake@input[seq_len(n_levels)]
profile_files <- snakemake@input[seq_len(n_levels) + n_levels]

for (i in seq_len(n_levels)) {
  flows_file   <- flows_files[[i]]
  profile_file <- profile_files[[i]]

  agg_lvl     <- basename(dirname(dirname(flows_file)))
  output_root <- dirname(dirname(dirname(flows_file)))

  flows    <- read.csv(flows_file,   sep = ";", na.strings = c("", "NA"),
                       stringsAsFactors = FALSE)
  profiles <- read.csv(profile_file, sep = ";", stringsAsFactors = FALSE)

  for (fao_division in snakemake@params$fao_divisions) {
    prod <- as.numeric(fao_division)

    # ── Chord diagram ──────────────────────────────────────────────────────
    for (ext in snakemake@params$ext) {
      out_dir    <- file.path(output_root, agg_lvl, "plot", fao_division)
      chord_path <- file.path(out_dir, paste0("chord_diagram.", ext))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      chord_plot(flows, profiles, prod, year_end, threshold, chord_path, ext)
    }

    # ── Trade network ──────────────────────────────────────────────────────
    net_plot <- network_plot(flows, profiles, prod, year_start, year_end, threshold)
    if (!is.null(net_plot)) {
      for (ext in snakemake@params$ext) {
        ggsave(
          filename   = paste0("trade_network.", ext),
          plot       = net_plot,
          device     = ext,
          path       = file.path(output_root, agg_lvl, "plot", fao_division),
          create.dir = TRUE,
          scale      = 1,
          width      = 2480 * 2.2,
          height     = 1240 * 1.5,
          units      = "px",
          dpi        = 300,
          bg         = "white"
        )
      }
    }
  }
}
