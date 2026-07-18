library("ggplot2")
library("hrbrthemes")
suppressPackageStartupMessages(library("dplyr"))
library("ggrepel")
library("patchwork")
library("scales")
library("reshape2")
library("legendry")
suppressPackageStartupMessages(library("circlize")) # nolint
suppressPackageStartupMessages(library("igraph"))

pal              <- c("#3D85F7", "#C32E5A")
trader_type_pal  <- c("main_exp" = "#C32E5A",
                      "balanced" = "#A49393",
                      "main_imp" = "#3D85F7")
row_col          <- "#E8D5B0"
tol_muted        <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", # nolint
                      "#DDCC77", "#CC6677", "#882255", "#AA4499", "#EE7733")

# Top-n traders by their share of global export or import value (ties broken by
# total traded value); every other country is pooled into "Rest of the World"
select_contributors <- function(profiles, prod, year, n_top) {
  d <- profiles %>% filter(cmd == prod, period == year) # nolint
  tot_exp <- sum(d$primary_value_exp, na.rm = TRUE)
  tot_imp <- sum(d$primary_value_imp, na.rm = TRUE)
  d %>% # nolint
    mutate(pv_exp = ifelse(is.na(primary_value_exp), 0, primary_value_exp), # nolint
           pv_imp = ifelse(is.na(primary_value_imp), 0, primary_value_imp), # nolint
           rank_key = pmax(pv_exp / tot_exp, pv_imp / tot_imp), # nolint
           tot_key  = pv_exp + pv_imp) %>% # nolint
    arrange(desc(rank_key), desc(tot_key)) %>% # nolint
    head(n_top) %>% # nolint
    pull(country) # nolint
}

# Classify from exp_share column (already in contributor_profiles.csv)
get_trader_type <- function(profiles, prod, year, countries) {
  d <- profiles %>% # nolint
    filter(cmd == prod, period == year, country %in% countries) # nolint
  types <- ifelse(d$exp_share >= 0.8, "main_exp",
                  ifelse(d$exp_share <= 0.2, "main_imp", "balanced"))
  setNames(types, d$country)
}

# Build square flow matrix from long-format bilateral flows
build_flow_matrix <- function(d_flows, value_col, countries) {
  mat <- matrix(0, nrow = length(countries), ncol = length(countries),
                dimnames = list(countries, countries))
  d_agg <- d_flows %>% # nolint
    group_by(exporter_desc, importer_desc) %>% # nolint
    summarise(val = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  for (k in seq_len(nrow(d_agg))) {
    e <- d_agg$exporter_desc[k]
    m <- d_agg$importer_desc[k]
    if (e %in% countries && m %in% countries) mat[e, m] <- d_agg$val[k]
  }
  mat
}

# Chord diagram. panels = "all" → 2×2 (primary value / net weight × exporter / importer reports); # nolint
#               panels = "fob" → 1×2 exporter reports only (primary value, net weight) # nolint
chord_plot <- function(flows, profiles, prod, year, n_top, output_path, ext,
                       panels = "all") { # nolint
  # Display alias applied to flows and profiles alike so selection stays consistent
  flows <- flows %>% # nolint
    mutate(exporter_desc = dplyr::recode(exporter_desc, "Russian Federation" = "Russia"), # nolint
           importer_desc = dplyr::recode(importer_desc, "Russian Federation" = "Russia")) # nolint
  profiles <- profiles %>% # nolint
    mutate(country = dplyr::recode(country, "Russian Federation" = "Russia")) # nolint

  contributors <- select_contributors(profiles, prod, year, n_top)
  if (length(contributors) < 2) return(invisible(NULL))

  d <- flows %>% # nolint
    filter(fao_code_agg == prod, period == year) %>% # nolint
    mutate(
      exporter_desc = ifelse(exporter_desc %in% contributors, # nolint
                             exporter_desc, "Rest of the World"),
      importer_desc = ifelse(importer_desc %in% contributors, # nolint
                             importer_desc, "Rest of the World"),
      pv_exp = ifelse(is.na(primary_value_exp), 0, primary_value_exp) / 1e6, # nolint
      pv_imp = ifelse(is.na(primary_value_imp), 0, primary_value_imp) / 1e6, # nolint
      nw_exp = ifelse(is.na(net_wgt_exp),       0, net_wgt_exp)       / 1e3, # nolint
      nw_imp = ifelse(is.na(net_wgt_imp),       0, net_wgt_imp)       / 1e3  # nolint
    )

  if (nrow(d) == 0) return(invisible(NULL))

  # arc_cols keyed by country name; reindexed per panel after metric-specific ordering # nolint
  n_contrib <- length(contributors)
  arc_cols  <- setNames(
    tol_muted[((seq_len(n_contrib) - 1L) %% length(tol_muted)) + 1L],
    contributors
  )
  arc_cols["Rest of the World"] <- row_col

  # Per-panel ordering: rank countries by total involvement in that metric (source + target) # nolint
  # so ranking shifts between panels highlight data discrepancies across metrics/reports # nolint
  all_ctry <- c(contributors, "Rest of the World")
  make_ordered_matrix <- function(d_flows, value_col) {
    d_agg <- d_flows %>% # nolint
      group_by(exporter_desc, importer_desc) %>% # nolint
      summarise(val = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>% # nolint
      filter(val > 0) # nolint
    t_vec <- rep(0, length(all_ctry))
    names(t_vec) <- all_ctry
    for (k in seq_len(nrow(d_agg))) {
      e <- d_agg$exporter_desc[k]
      i <- d_agg$importer_desc[k]
      if (e %in% all_ctry) t_vec[e] <- t_vec[e] + d_agg$val[k]
      if (i %in% all_ctry) t_vec[i] <- t_vec[i] + d_agg$val[k]
    }
    non_row <- setdiff(all_ctry, "Rest of the World")
    ord <- c(non_row[order(t_vec[non_row], decreasing = TRUE)], "Rest of the World") # nolint
    build_flow_matrix(d_flows, value_col, ord)
  }

  mat_pv_exp <- make_ordered_matrix(d, "pv_exp")
  mat_pv_imp <- make_ordered_matrix(d, "pv_imp")
  mat_nw_exp <- make_ordered_matrix(d, "nw_exp")
  mat_nw_imp <- make_ordered_matrix(d, "nw_imp")

  # Long country names split to two lines for legibility around the arc
  label_map <- c(
    "European Union"    = "European\nUnion",
    "New Zealand"       = "New\nZealand",
    "Papua New Guinea"  = "Papua N.G.",
    "Rest of the World" = "Rest of\nthe World"
  )
  wrap_label <- function(x) ifelse(x %in% names(label_map), label_map[x], x)

  country_labels <- function() {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) { # nolint
      xlim <- get.cell.meta.data("xlim") # nolint
      lbl  <- wrap_label(get.cell.meta.data("sector.index")) # nolint
      circos.text(mean(xlim), 0.5, lbl, # nolint
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0.5), cex = 0.7)
    }, bg.border = NA)
  }

  draw_panel <- function(mat, panel_label, subtitle,
                         mar = c(0.3, 0.3, 3.0, 0.3), canvas_lim = 1.15,
                         cex_main = 0.85, title_line = NA) {
    panel_arc_cols <- arc_cols[rownames(mat)]  # reorder colours to match panel ordering # nolint
    circos.clear() # nolint
    par(mar = mar)
    circos.par(canvas.xlim = c(-canvas_lim, canvas_lim),
               canvas.ylim = c(-canvas_lim, canvas_lim)) # nolint
    chordDiagram(mat, grid.col = panel_arc_cols, transparency = 0.4, # nolint
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.1))
    country_labels()
    title(paste(panel_label, subtitle), adj = 0, cex.main = cex_main, font.main = 2,
          line = title_line)
    circos.clear() # nolint
  }

  if (panels == "fob") {
    # Main-text figure: 1×2, exporter reports only (FOB), 190 mm × 100 mm at 600 DPI
    px_w <- round(190 / 25.4 * 600); px_h <- round(100 / 25.4 * 600)
    if (ext == "png") {
      png(output_path, width = px_w, height = px_h, res = 600, bg = "white")
    } else {
      svg(output_path, width = 190 / 25.4, height = 100 / 25.4, bg = "white")
    }
    par(mfrow = c(1, 2), oma = c(0, 0, 0, 0))
    draw_panel(mat_pv_exp, "(a)",
               "Primary value — exporter reports (FOB, million USD)",
               mar = c(0.3, 0.6, 2.5, 0.7), canvas_lim = 1.15, cex_main = 0.65,
               title_line = 1.0)
    draw_panel(mat_nw_exp, "(b)",
               "Net weight — exporter reports (tonnes)",
               mar = c(0.7, 0.7, 2.5, 0.6), canvas_lim = 1.15, cex_main = 0.65,
               title_line = 1.0)
  } else {
    # Full four-panel figure (supplementary): 190 mm × 190 mm at 600 DPI, 2×2 square grid
    px_side <- round(190 / 25.4 * 600)  # 4488 px
    if (ext == "png") {
      png(output_path, width = px_side, height = px_side, res = 600, bg = "white")
    } else {
      in_side <- 190 / 25.4
      svg(output_path, width = in_side, height = in_side, bg = "white")
    }
    par(mfrow = c(2, 2), oma = c(0, 0, 0, 0))
    draw_panel(mat_pv_exp, "(a)",
               "Primary value — exporter reports (FOB, million USD)")
    draw_panel(mat_pv_imp, "(b)",
               "Primary value — importer reports (CIF, million USD)")
    draw_panel(mat_nw_exp, "(c)",
               "Net weight — exporter reports (tonnes)")
    draw_panel(mat_nw_imp, "(d)",
               "Net weight — importer reports (tonnes)")
  }

  dev.off()
}

# Size-aware layout repulsion.
# coords:     n×2 matrix of raw igraph coordinates.
# node_sizes: per-node max bubble value (same units as size_max).
# size_max:   global max used for the ggplot size scale (range c(1,20)).
# panel_mm:   panel width in mm (converts ggplot size units → coord units).
# margin:     clearance factor on top of the sum-of-radii minimum distance.
repel_layout <- function(coords, node_sizes, size_max,
                         panel_mm = 80, margin = 1.3, max_iter = 1000) {
  n <- nrow(coords)
  if (n < 2) return(coords)
  # Normalise preserving aspect ratio (same scale for x and y)
  sc <- max(diff(range(coords[, 1])), diff(range(coords[, 2])))
  if (sc > 0) {
    coords[, 1] <- (coords[, 1] - min(coords[, 1])) / sc
    coords[, 2] <- (coords[, 2] - min(coords[, 2])) / sc
  }
  # Map node values to ggplot size units (linear scale c(1,20))
  ggplot_sizes <- 1 + (pmin(node_sizes, size_max) / size_max) * 19
  # ggplot size ≈ diameter in mm  →  radius in [0,1] coord units
  radii <- (ggplot_sizes / 2) / panel_mm
  for (iter in seq_len(max_iter)) {
    moved <- FALSE
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        dx <- coords[j, 1] - coords[i, 1]
        dy <- coords[j, 2] - coords[i, 2]
        d  <- sqrt(dx^2 + dy^2)
        needed <- (radii[i] + radii[j]) * margin
        if (d < needed) {
          if (d < 1e-9) { dx <- 1e-4; dy <- 1e-4; d <- sqrt(2) * 1e-4 }
          push <- (needed - d) / 2
          coords[i, 1] <- coords[i, 1] - push * dx / d
          coords[i, 2] <- coords[i, 2] - push * dy / d
          coords[j, 1] <- coords[j, 1] + push * dx / d
          coords[j, 2] <- coords[j, 2] + push * dy / d
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  coords
}

# One ggplot2 network panel: nested bubbles (same encoding as contributor_profiles) # nolint
# with straight-line arrows underneath. Node sizes computed from flows_raw.
# show_legend = TRUE only on panel (b) to avoid duplicate guide_circles crash in patchwork. # nolint
build_network_panel <- function(flows, flows_raw, prod, year, contributors,
                                layout_df, size_max, panel_letter,
                                show_legend = FALSE, node_threshold = 0.01,
                                panel_title = paste0(panel_letter, " ", year)) {
  # Node sizes from flows_raw (covers all countries, not just profiled contributors) # nolint
  d_exp_agg <- flows_raw %>% # nolint
    group_by(country = exporter_desc) %>% # nolint
    summarise(size_exp = sum(ifelse(is.na(primary_value_exp), 0, primary_value_exp)) / 1e6, # nolint
              .groups = "drop")
  d_imp_agg <- flows_raw %>% # nolint
    group_by(country = importer_desc) %>% # nolint
    summarise(size_imp = sum(ifelse(is.na(primary_value_imp), 0, primary_value_imp)) / 1e6, # nolint
              .groups = "drop")
  d_nodes <- full_join(d_exp_agg, d_imp_agg, by = "country") %>% # nolint
    mutate(size_exp = ifelse(is.na(size_exp), 0, size_exp), # nolint
           size_imp = ifelse(is.na(size_imp), 0, size_imp)) # nolint

  # Countries above threshold: nodes, edges, and labels all filtered to this set
  tot_exp_m <- sum(flows_raw$primary_value_exp, na.rm = TRUE) / 1e6
  tot_imp_m <- sum(flows_raw$primary_value_imp, na.rm = TRUE) / 1e6
  above_thr <- d_nodes %>% # nolint
    filter((tot_exp_m > 0 & size_exp / tot_exp_m >= node_threshold) | # nolint
           (tot_imp_m > 0 & size_imp / tot_imp_m >= node_threshold)) %>% # nolint
    pull(country) # nolint

  # Edges: binary (link exists if either report > 0), one row per directed pair
  d_flows <- flows %>% # nolint
    filter(fao_code_agg == prod, period == year, # nolint
           exporter_desc %in% above_thr, importer_desc %in% above_thr) %>% # nolint
    filter((!is.na(primary_value_exp) & primary_value_exp > 0) | # nolint
           (!is.na(primary_value_imp) & primary_value_imp > 0)) %>% # nolint
    distinct(exporter_desc, importer_desc) %>% # nolint
    left_join(layout_df %>% rename(exporter_desc = country, x0 = layout_x, y0 = layout_y), # nolint
              by = "exporter_desc") %>% # nolint
    left_join(layout_df %>% rename(importer_desc = country, x1 = layout_x, y1 = layout_y), # nolint
              by = "importer_desc") %>% # nolint
    filter(!is.na(x0), !is.na(x1)) %>% # nolint
    mutate(xm = (x0 + x1) / 2, ym = (y0 + y1) / 2) # nolint

  # Long format: two rows per country; 
  # split so smallest renders on top at alpha = 1
  d_long <- bind_rows(
    d_nodes %>% filter(country %in% above_thr) %>% # nolint
      mutate(trade_type = "Total exported value", size_val = size_exp), # nolint
    d_nodes %>% filter(country %in% above_thr) %>% # nolint
      mutate(trade_type = "Total imported value", size_val = size_imp) # nolint
  ) %>% # nolint
    arrange(country, desc(size_val)) %>% # nolint
    group_by(country) %>% # nolint
    mutate(is_top = row_number() == 2L) %>% # nolint
    ungroup() %>% # nolint
    left_join(layout_df, by = "country")
  d_long_bottom <- d_long %>% filter(!is_top) # nolint
  d_long_top    <- d_long %>% filter(is_top)  # nolint

  d_labels <- d_nodes %>% # nolint
    filter(country %in% above_thr) %>% # nolint
    left_join(layout_df, by = "country") %>% # nolint
    mutate(country_label = dplyr::recode(country,
      "European Union"   = "European\nUnion",
      "Papua New Guinea"  = "Papua New\nGuinea",
      "Rep. of Korea"    = "Rep. of\nKorea"
    )) # nolint

  # White background disc sized to the larger of the two bubbles per country
  d_white <- d_nodes %>% # nolint
    filter(country %in% above_thr) %>% # nolint
    mutate(size_white = pmax(size_exp, size_imp)) %>% # nolint
    left_join(layout_df, by = "country")

  ggplot() +
    geom_segment(data = d_flows,
                 aes(x = x0, y = y0, xend = x1, yend = y1), # nolint
                 linewidth = 0.3, alpha = 0.25, color = "grey50",
                 show.legend = FALSE) +
    geom_point(data = d_white,
               aes(x = layout_x, y = layout_y, size = size_white), # nolint
               color = "white", alpha = 1.0,
               show.legend = FALSE) +
    geom_point(data = d_long_bottom,
               aes(x = layout_x, y = layout_y, size = size_val, color = trade_type), # nolint
               alpha = 0.6,
               show.legend = FALSE) +
    geom_point(data = d_long_top,
               aes(x = layout_x, y = layout_y, size = size_val, color = trade_type), # nolint
               alpha = 0.6,
               show.legend = show_legend) +
    geom_text_repel(data = d_labels,
                    aes(x = layout_x, y = layout_y, label = country_label), # nolint
                    size = 2.5, color = "black",
                    bg.color = "white", bg.r = 0.2,
                    segment.colour = "grey80", segment.size = 0.25,
                    min.segment.length = 0,
                    box.padding = 0.6, point.padding = 0.4,
                    force = 4, seed = 42,
                    direction = "both", max.overlaps = Inf,
                    show.legend = FALSE) +
    scale_color_manual(
      values = c("Total exported value" = pal[2], "Total imported value" = pal[1]), # nolint
      name   = NULL
    ) +
    { # nolint
      if (show_legend) {
        size_breaks <- unique(round(c(size_max * 0.2, size_max * 0.5, size_max) / 1000) * 1000) # nolint
        size_breaks <- size_breaks[size_breaks > 0]
        scale_size_continuous(range  = c(1, 20),
                              breaks = size_breaks,
                              limits = c(0, size_max),
                              labels = scales::label_comma(),
                              name   = "Traded value (million US$)")
      } else {
        scale_size_continuous(range  = c(1, 20),
                              limits = c(0, size_max),
                              guide  = "none")
      }
    } +
    { # nolint
      if (show_legend) {
        guides(
          colour = guide_legend(title        = NULL,
                               override.aes = list(size = 8, alpha = 0.7), # nolint
                               theme        = theme( # nolint
                                 legend.margin = margin(t = 20, unit = "pt") # nolint
                               )),
          size   = guide_circles(
            theme         = theme(
              legendry.legend.key.margin = margin(5, 5, 5, 5) # nolint
            ),
            text_position = "right",
            override.aes  = aes(colour = "grey50", alpha = 0.5)
          )
        )
      } else {
        guides(colour = "none", size = "none")
      }
    } +
    theme_void() +
    theme(legend.position    = "bottom",
          legend.box         = "horizontal",
          legend.box.just    = "center",
          legend.margin      = margin(5, 5, 5, 5),
          plot.subtitle      = element_text(size = 7, face = "bold",
                                              margin = margin(b = 8)), # nolint
          plot.background    = element_rect(fill = "white", color = NA),
          plot.margin        = margin(10, 20, 10, 20),
          legend.text        = element_text(size = 7),
          legend.title       = element_text(size = 7),
          legend.key.size    = unit(0.40, "cm"),
          legend.spacing.x   = unit(0.5, "cm")) +
    labs(subtitle = panel_title)
}

# Two-panel network figure (year_start left, year_end right); all countries, no threshold # nolint
network_plot <- function(flows, prod, year_start, year_end, node_threshold = 0.01) { # nolint
  flows <- flows %>% # nolint
    mutate(
      exporter_desc = dplyr::recode(exporter_desc, # nolint
                                    "Russian Federation" = "Russia",
                                    "China, Hong Kong SAR" = "Hong Kong"),
      importer_desc = dplyr::recode(importer_desc, # nolint
                                    "Russian Federation" = "Russia",
                                    "China, Hong Kong SAR" = "Hong Kong")
    )
  flows_raw_start <- flows %>% filter(fao_code_agg == prod, period == year_start) # nolint
  flows_raw_end   <- flows %>% filter(fao_code_agg == prod, period == year_end)   # nolint
  if (nrow(flows_raw_end) == 0) return(NULL)

  # Above-threshold countries for one year (mirrors build_network_panel logic)
  year_above_thr <- function(d_raw) {
    d_e <- d_raw %>% group_by(country = exporter_desc) %>% # nolint
      summarise(se = sum(ifelse(is.na(primary_value_exp), 0, primary_value_exp)) / 1e6, # nolint
                .groups = "drop")
    d_i <- d_raw %>% group_by(country = importer_desc) %>% # nolint
      summarise(si = sum(ifelse(is.na(primary_value_imp), 0, primary_value_imp)) / 1e6, # nolint
                .groups = "drop")
    dn  <- full_join(d_e, d_i, by = "country") %>%
      mutate(se = ifelse(is.na(se), 0, se), si = ifelse(is.na(si), 0, si)) # nolint
    te  <- sum(d_raw$primary_value_exp, na.rm = TRUE) / 1e6
    ti  <- sum(d_raw$primary_value_imp, na.rm = TRUE) / 1e6
    dn %>% filter((te > 0 & se / te >= node_threshold) | # nolint
                  (ti > 0 & si / ti >= node_threshold)) %>% pull(country) # nolint
  }

  # Layout uses union of above-threshold countries from both years for stable positions # nolint
  layout_countries <- union(year_above_thr(flows_raw_start), year_above_thr(flows_raw_end)) # nolint
  if (length(layout_countries) < 2) return(NULL)

  # Layout anchored to year_end; binary edges (either report > 0)
  d_edges_end <- flows_raw_end %>% # nolint
    filter(exporter_desc %in% layout_countries, importer_desc %in% layout_countries) %>% # nolint
    filter((!is.na(primary_value_exp) & primary_value_exp > 0) | # nolint
           (!is.na(primary_value_imp) & primary_value_imp > 0)) %>% # nolint
    distinct(exporter_desc, importer_desc) %>% # nolint
    rename(from = exporter_desc, to = importer_desc)

  if (nrow(d_edges_end) == 0) return(NULL)

  # Size maximum across both years — needed before layout for size-aware repulsion # nolint
  d_both_flows <- bind_rows(flows_raw_start, flows_raw_end) %>% # nolint
    filter(exporter_desc %in% layout_countries | importer_desc %in% layout_countries) # nolint
  size_max <- ceiling(max(
    tapply(ifelse(is.na(d_both_flows$primary_value_exp), 0, d_both_flows$primary_value_exp), # nolint
           d_both_flows$exporter_desc, sum),
    tapply(ifelse(is.na(d_both_flows$primary_value_imp), 0, d_both_flows$primary_value_imp), # nolint
           d_both_flows$importer_desc, sum),
    na.rm = TRUE
  ) / 1e6)

  # Per-node max bubble size (max of exp and imp across both years)
  d_nsize_e <- d_both_flows %>% # nolint
    group_by(country = exporter_desc) %>% # nolint
    summarise(v = sum(ifelse(is.na(primary_value_exp), 0, primary_value_exp)) / 1e6, # nolint
              .groups = "drop")
  d_nsize_i <- d_both_flows %>% # nolint
    group_by(country = importer_desc) %>% # nolint
    summarise(v = sum(ifelse(is.na(primary_value_imp), 0, primary_value_imp)) / 1e6, # nolint
              .groups = "drop")
  d_nsizes <- bind_rows(d_nsize_e, d_nsize_i) %>% # nolint
    group_by(country) %>% summarise(size_val = max(v), .groups = "drop") # nolint

  g <- igraph::graph_from_data_frame(
    d_edges_end, directed = TRUE,
    vertices = data.frame(name = layout_countries)
  )
  node_names    <- igraph::V(g)$name
  node_sizes_v  <- d_nsizes$size_val[match(node_names, d_nsizes$country)]
  node_sizes_v[is.na(node_sizes_v)] <- 0

  set.seed(42)
  coords    <- igraph::layout_with_kk(g)
  coords    <- repel_layout(coords, node_sizes = node_sizes_v,
                            size_max = size_max, panel_mm = 80)
  layout_df <- data.frame(country  = node_names,
                           layout_x = coords[, 1], # nolint
                           layout_y = coords[, 2]) # nolint

  prod_names <- c(
    "1" = "Roundwood", "5" = "Sawnwood", "7" = "Wood-based panels"
  )
  prod_label <- prod_names[as.character(prod)]
  if (is.na(prod_label)) prod_label <- paste("Product", prod)

  p_start <- build_network_panel(
    flows, flows_raw_start, prod, year_start, layout_countries,
    layout_df, size_max, "(a)",
    show_legend    = FALSE,
    node_threshold = node_threshold,
    panel_title    = paste0("(a) ", prod_label, " trade network in ",
                            year_start)
  )
  p_end <- build_network_panel(
    flows, flows_raw_end, prod, year_end, layout_countries,
    layout_df, size_max, "(b)",
    show_legend    = TRUE,
    node_threshold = node_threshold,
    panel_title    = paste0("(b) ", prod_label, " trade network in ",
                            year_end)
  )

  (p_start + p_end) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "top") # nolint
}

# ── Main loop ───────────────────────────────────────────────────────────────

year_start <- snakemake@params$year_start
year_end   <- snakemake@params$year_end
chord_year <- snakemake@params$chord_year

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
      out_dir        <- file.path(output_root, agg_lvl, "plot", fao_division)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      chord_path     <- file.path(out_dir, paste0("chord_diagram.", ext))
      chord_fob_path <- file.path(out_dir, paste0("chord_diagram_fob.", ext))
      chord_plot(flows, profiles, prod, chord_year, snakemake@params$chord_n, # nolint
                 chord_path, ext, panels = "all")
      chord_plot(flows, profiles, prod, chord_year, snakemake@params$chord_n, # nolint
                 chord_fob_path, ext, panels = "fob")
    }

    # ── Trade network ──────────────────────────────────────────────────────
    net_plot <- network_plot(flows, prod, year_start, year_end,
                             node_threshold = 0.01)
    if (!is.null(net_plot)) {
      for (ext in snakemake@params$ext) {
        ggsave(
          filename   = paste0("trade_network.", ext),
          plot       = net_plot,
          device     = ext,
          path       = file.path(output_root, agg_lvl, "plot", fao_division),
          create.dir = TRUE,
          scale      = 1,
          width      = 190,
          height     = 130,
          units      = "mm",
          dpi        = 600,
          bg         = "white"
        )
      }
    }
  }
}
