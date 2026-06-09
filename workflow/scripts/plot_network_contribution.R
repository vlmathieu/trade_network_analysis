library("ggplot2")
library("hrbrthemes")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
library("patchwork")
library("tibble")
library("purrr")
# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

benchmark_years <- c(2000, 2005, 2010, 2015, 2020)

trader_titles <- c(
  "main_exp" = "Main exporters",
  "balanced" = "Balanced",
  "main_imp" = "Main importers"
)

wgt_titles <- c(
  "primary_value" = "Primary value",
  "net_wgt"       = "Net weight"
)

y_axis_labels <- c(
  "primary_value" = "Share of global trade value (%)",
  "net_wgt"       = "Share of global trade volume (%)"
)

# Base colour per trader type — consistent with plot_network_composition.R
pal_trader <- c(
  "main_exp" = "#C32E5A",
  "balanced" = "#A49393",
  "main_imp" = "#3D85F7"
)

# Sequential ranges (light → base → dark) used to generate per-country shades
pal_ranges <- list(
  "main_exp" = c("#F5B8C8", "#C32E5A", "#7A1130"),
  "balanced" = c("#D9D0D0", "#A49393", "#5C4D4D"),
  "main_imp" = c("#B8D0FC", "#3D85F7", "#1040A0")
)

# ---------------------------------------------------------------------------
# compute_label_positions() # nolint
# ---------------------------------------------------------------------------
# Resolves vertical collisions for N end-of-line labels placed in a right-hand
# strip. Applies a greedy upward push then shifts the whole stack down if it
# overflows y_max, matching the strategy used in plot_market_concentration.R.

compute_label_positions <- function(y_vals, y_max, min_gap_frac = 0.06) {
  if (length(y_vals) == 0) return(y_vals)
  ord <- order(y_vals)
  y_s <- y_vals[ord]
  nms <- names(y_vals)[ord]
  gap <- y_max * min_gap_frac
  for (i in seq_along(y_s)[-1]) {
    if (y_s[i] - y_s[i - 1] < gap) y_s[i] <- y_s[i - 1] + gap
  }
  if (y_s[length(y_s)] > y_max) y_s <- y_s - (y_s[length(y_s)] - y_max)
  if (y_s[1] < 0)               y_s <- y_s - y_s[1]
  setNames(y_s, nms)
}

# ---------------------------------------------------------------------------
# select_countries() # nolint
# ---------------------------------------------------------------------------
# Returns a character vector of countries that were major contributors over
# the last `time_span` years, based on primary value only.
# Selection criterion: max(contrib_tot_exp, contrib_tot_imp) exceeds the
# (1 - min_contrib) quantile across all countries in the time window.

select_countries <- function(data,
                             prod,
                             time_span = 5,
                             top_frac  = 0.03) {

  max_year <- max(data$period)

  data %>% # nolint
    dplyr::filter(
      cmd    == prod,                        # nolint
      weight == "primary_value",             # nolint
      period >= max_year - time_span         # nolint
    ) %>% # nolint
    rowwise() %>% # nolint
    mutate(max_contrib = max(contrib_tot_exp, contrib_tot_imp)) %>% # nolint
    ungroup() %>% # nolint
    dplyr::filter(
      max_contrib > quantile(.$max_contrib, probs = 1 - top_frac) # nolint
    ) %>% # nolint
    dplyr::select(country) %>% # nolint
    dplyr::distinct() %>% # nolint
    dplyr::pull(country)
}

# ---------------------------------------------------------------------------
# classify_countries() # nolint
# ---------------------------------------------------------------------------
# Returns a named character vector (country -> group) where group is one of
# "main_exp", "balanced", "main_imp".
#
# Classification is read directly from network_composition.csv: for each
# (cmd, period) row, countries are listed in one of three columns:
#   list_main_exp, list_main_imp, list_balanced (comma-separated strings).
#
# The per-year classification is aggregated to a single label per country
# by taking the modal category (most frequent across years). Ties are broken
# in favour of the more extreme category: main_exp > main_imp > balanced.

classify_countries <- function(composition,
                               prod,
                               selected_countries) {

  # Helper: modal category with tiebreak main_exp > main_imp > balanced
  modal_group <- function(groups) {
    tab <- table(groups)
    priority <- c("main_exp" = 1, "main_imp" = 2, "balanced" = 3)
    max_count <- max(tab)
    tied <- names(tab)[tab == max_count]
    tied[which.min(priority[tied])]
  }

  # Filter composition to this product
  comp_prod <- composition %>% # nolint
    dplyr::filter(cmd == prod) # nolint

  # Build long-format table: one row per (period, country, group)
  # by splitting the comma-separated list columns
  rows <- list()
  for (i in seq_len(nrow(comp_prod))) {
    period_i <- comp_prod$period[i] # nolint
    for (col_grp in list(
      list(col = "list_main_exp", grp = "main_exp"),
      list(col = "list_main_imp", grp = "main_imp"),
      list(col = "list_balanced", grp = "balanced")
    )) {
      raw <- comp_prod[[col_grp$col]][i]
      if (!is.na(raw) && nchar(trimws(raw)) > 0) {
        countries_in_list <- trimws(strsplit(raw, ",")[[1]])
        matched <- intersect(countries_in_list, selected_countries)
        if (length(matched) > 0) {
          rows[[length(rows) + 1]] <- data.frame(
            country = matched,
            group   = col_grp$grp,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(rows) == 0) return(setNames(character(0), character(0)))

  long <- dplyr::bind_rows(rows)

  # Modal group per country
  long %>% # nolint
    dplyr::group_by(country) %>%                              # nolint
    dplyr::summarise(group = modal_group(group), # nolint
                     .groups = "drop") %>% # nolint
    tibble::deframe()  # named vector: country -> group
}

# ---------------------------------------------------------------------------
# shape_data() # nolint
# ---------------------------------------------------------------------------
# Filters to the selected countries and weight type, computes contrib/
# contrib_min/contrib_max (x100), and fits three LOESS smooths per country.
# Countries with fewer than 5 non-NA observations are silently dropped before
# fitting to avoid LOESS failures.

shape_data <- function(data,
                       prod,
                       wgt,
                       selected_countries,
                       span = 0.5) {

  shaped <- data %>% # nolint
    dplyr::filter(
      cmd    == prod,                      # nolint
      weight == wgt,                       # nolint
      country %in% selected_countries      # nolint
    ) %>% # nolint
    tibble::as_tibble() %>% # nolint
    dplyr::arrange(country, period) %>%    # nolint
    mutate(
      contrib     = (contrib_tot_exp + contrib_tot_imp) / 2 * 100,              # nolint
      contrib_min = pmin(contrib_tot_exp, contrib_tot_imp, na.rm = TRUE) * 100, # nolint
      contrib_max = pmax(contrib_tot_exp, contrib_tot_imp, na.rm = TRUE) * 100  # nolint
    )

  # Drop countries with too few observations for LOESS
  shaped <- shaped %>% # nolint
    dplyr::group_by(country) %>%           # nolint
    dplyr::filter(sum(!is.na(contrib)) >= 5) %>% # nolint
    dplyr::ungroup()

  # Fit three LOESS smooths per country via nest/map.
  # contrib_min/max: fit on non-NA rows only, then predict on the full grid so
  # that years where only one trade direction is available don't create gaps.
  shaped <- shaped %>% # nolint
    dplyr::group_by(country) %>%           # nolint
    tidyr::nest() %>% # nolint
    mutate(
      contrib_pred     = purrr::map(data, function(d) {
        predict(loess(contrib ~ period, data = d, span = span), newdata = d)
      }),
      contrib_min_pred = purrr::map(data, function(d) {
        d_fit <- d[!is.na(d$contrib_min), ]
        if (nrow(d_fit) < 4) return(rep(NA_real_, nrow(d)))
        predict(loess(contrib_min ~ period, data = d_fit, span = span), newdata = d) # nolint
      }),
      contrib_max_pred = purrr::map(data, function(d) {
        d_fit <- d[!is.na(d$contrib_max), ]
        if (nrow(d_fit) < 4) return(rep(NA_real_, nrow(d)))
        predict(loess(contrib_max ~ period, data = d_fit, span = span), newdata = d) # nolint
      })
    ) %>% # nolint
    tidyr::unnest(cols = c(data, contrib_pred,             # nolint
                           contrib_min_pred, contrib_max_pred)) %>% # nolint
    mutate(contrib_min_pred = pmax(0, contrib_min_pred)) # nolint

  return(shaped) # nolint
}

# ---------------------------------------------------------------------------
# build_panel() # nolint
# ---------------------------------------------------------------------------
# Builds one ggplot panel for a given (trader_type, wgt) combination.
# Returns NULL if no countries from group_countries are present in plot_data.

build_panel <- function(plot_data,
                        group_countries,
                        trader_type,
                        wgt,
                        y_max,
                        x_max,
                        panel_letter) {

  panel_data <- plot_data %>% # nolint
    dplyr::filter(country %in% group_countries) # nolint

  if (nrow(panel_data) == 0) return(NULL)

  label_data <- panel_data %>% # nolint
    dplyr::filter(period == x_max, !is.na(contrib_pred)) # nolint

  y_actual  <- setNames(label_data$contrib_pred, label_data$country)
  y_nudged  <- compute_label_positions(y_actual, y_max)
  label_df  <- label_data %>% # nolint
    dplyr::mutate(y_nudged = y_nudged[country]) # nolint

  panel_title <- paste0(
    panel_letter, " ",
    trader_titles[trader_type], " \u2014 ", wgt_titles[wgt]
  )

  # Per-country colour shades and shapes derived from the trader-type base colour # nolint
  countries      <- sort(unique(panel_data$country))
  n_countries    <- length(countries)
  country_cols   <- if (n_countries == 1) {
    setNames(pal_trader[trader_type], countries)
  } else {
    setNames(colorRampPalette(pal_ranges[[trader_type]])(n_countries), countries) # nolint
  }
  country_shapes <- setNames(
    c(16, 17, 15, 18, 8, 25)[seq_len(n_countries)],
    countries
  )

  ggplot(panel_data,
         aes(x = period, group = country)) + # nolint
    # HS-revision discontinuity band (net weight panels only)
    (if (wgt == "net_wgt") geom_rect(
      aes(xmin = 2000, xmax = 2006, ymin = 0, ymax = y_max), # nolint
      fill        = "grey85",
      alpha       = 0.08,
      inherit.aes = FALSE
    ) else NULL) +
    # Benchmark year gridlines
    geom_vline(xintercept = benchmark_years,
               linetype   = "dashed",
               color      = "grey70",
               linewidth  = 0.3) +
    # Smoothed ribbon (source vs target gap)
    geom_ribbon(aes(ymin = contrib_min_pred, ymax = contrib_max_pred, # nolint
                    fill = country),                                   # nolint
                alpha = 0.3) +
    # Smoothed midpoint line
    geom_line(aes(y = contrib_pred, color = country)) +    # nolint
    # Smoothed midpoint points \u2014 shape varies by country
    geom_point(aes(y = contrib_pred, color = country, shape = country), # nolint
               size = 0.8, alpha = 0.6) +
    # Connector segments from ribbon-top anchor to nudged label position
    geom_segment(
      data        = label_df, # nolint
      mapping     = aes(x    = x_max + 0.15, xend = x_max + 0.35, # nolint
                        y    = contrib_pred, yend = y_nudged,       # nolint
                        color = country),                           # nolint
      linewidth   = 0.35, # nolint
      inherit.aes = FALSE) +
    # Country name labels at nudged positions
    geom_text(
      data        = label_df, # nolint
      mapping     = aes(x     = x_max + 0.4, y = y_nudged,         # nolint
                        label = country, color = country),          # nolint
      hjust       = 0, # nolint
      fontface    = "bold",
      size        = 4,
      inherit.aes = FALSE) +
    # Scales
    scale_color_manual(values = country_cols) +
    scale_fill_manual(values  = country_cols) +
    scale_shape_manual(values = country_shapes) +
    scale_x_continuous(
      breaks = c(1996, benchmark_years, x_max),
      labels = as.character(c(1996, benchmark_years, x_max)),
      expand = expansion(mult = c(0.02, 0.20))
    ) +
    scale_y_continuous(
      breaks = seq(0, y_max, by = 10),
      expand = c(0, 0)
    ) +
    labs(
      x     = "Year",
      y     = y_axis_labels[wgt],
      title = panel_title
    ) +
    coord_cartesian(ylim = c(0, y_max), clip = "off") +
    theme_ipsum(base_size = 14, axis_title_size = 14) +
    theme(
      legend.position  = "none",
      plot.margin      = margin(3, 35, 3, 5, unit = "mm"),
      plot.title       = element_text(size = 14, face = "bold"),
      panel.ontop      = TRUE,
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
}

# ---------------------------------------------------------------------------
# Main loops
# ---------------------------------------------------------------------------

# Input files: first half are network_contribution.csv (one per agg_lvl),
# second half are network_composition.csv (same order, same agg_lvl).
# snakemake@input is a flat list: [contrib_1, contrib_2, ..., compo_1, compo_2, ...] # nolint
n_levels     <- length(snakemake@input) / 2
contrib_files <- snakemake@input[seq_len(n_levels)]
compo_files   <- snakemake@input[seq_len(n_levels) + n_levels]

# Panel letter sequence: (a)-(f), left-to-right then top-to-bottom
# Row order: main_exp, balanced, main_imp
# Column order: primary_value, net_wgt
trader_types  <- c("main_exp", "balanced", "main_imp")
weights       <- c("primary_value", "net_wgt")
panel_letters <- letters[1:6]  # a-f

for (idx in seq_len(n_levels)) {

  contrib_file <- contrib_files[[idx]]
  compo_file   <- compo_files[[idx]]

  # Derive agg_lvl and output root from the contribution input path
  agg_lvl     <- basename(dirname(dirname(contrib_file)))
  output_root <- dirname(dirname(dirname(contrib_file)))

  # Load data
  network_contribution <- read.csv(file   = contrib_file,
                                   header = TRUE,
                                   sep    = ";")
  network_composition  <- read.csv(file   = compo_file,
                                   header = TRUE,
                                   sep    = ";",
                                   colClasses = c(cmd = "character"))

  for (fao_division in snakemake@params$fao_divisions) {

    prod <- as.numeric(fao_division)

    # --- Country selection (primary value only) ---
    selected_countries <- select_countries(
      data      = network_contribution,
      prod      = prod,
      time_span = snakemake@params$time_span,
      top_frac  = snakemake@params$top_frac
    )

    # --- Country classification from network_composition ---
    # cmd in network_composition is zero-padded (e.g. "01"), so format prod
    country_groups <- classify_countries(
      composition         = network_composition,
      prod                = sprintf("%02d", prod),
      selected_countries  = selected_countries
    )

    # --- Shape data for both weight types ---
    shaped <- list(
      primary_value = shape_data(network_contribution, prod,
                                 "primary_value", selected_countries),
      net_wgt       = shape_data(network_contribution, prod,
                                 "net_wgt",       selected_countries)
    )

    # --- Shared y-axis maximum across all 6 panels (data-driven, multiple of 10) --- # nolint
    y_max <- ceiling(max(
      shaped$primary_value$contrib_max_pred,
      shaped$primary_value$contrib_max,
      shaped$net_wgt$contrib_max_pred,
      shaped$net_wgt$contrib_max,
      na.rm = TRUE
    ) / 10) * 10

    x_max <- max(shaped$primary_value$period)

    # --- Build 6 panels ---
    # Letter index runs left-to-right, top-to-bottom:
    # (a) main_exp/pv  (b) main_exp/nw
    # (c) balanced/pv  (d) balanced/nw
    # (e) main_imp/pv  (f) main_imp/nw
    panels     <- list()
    letter_idx <- 1

    for (trader_type in trader_types) {
      group_countries <- names(country_groups)[country_groups == trader_type]

      for (wgt in weights) {
        panels[[paste(trader_type, wgt, sep = "_")]] <- build_panel(
          plot_data       = shaped[[wgt]],
          group_countries = group_countries,
          trader_type     = trader_type,
          wgt             = wgt,
          y_max           = y_max,
          x_max           = x_max,
          panel_letter    = paste0("(", panel_letters[letter_idx], ")")
        )
        letter_idx <- letter_idx + 1
      }
    }

    # --- Assemble patchwork, skipping empty rows ---
    # NULL panels replaced by theme_void() placeholders to keep grid aligned.
    # Rows assembled explicitly (not via Reduce) to avoid patchwork nesting
    # issues with the / operator on pre-assembled row objects.

    fill_null <- function(p) {
      if (is.null(p)) ggplot() + theme_void() else p
    }

    row_list <- list()
    for (tt in trader_types) {
      p_left  <- panels[[paste(tt, "primary_value", sep = "_")]]
      p_right <- panels[[paste(tt, "net_wgt",       sep = "_")]]
      if (!is.null(p_left) || !is.null(p_right)) {
        row_list[[tt]] <- fill_null(p_left) + fill_null(p_right)
      }
    }

    scale_rows <- length(row_list)
    if (scale_rows == 0) next  # nothing to plot for this division

    patchwork_plot <- if (scale_rows == 1) {
      row_list[[1]]
    } else if (scale_rows == 2) {
      row_list[[1]] / row_list[[2]]
    } else {
      row_list[[1]] / row_list[[2]] / row_list[[3]]
    }

    # --- Save ---
    for (ext in snakemake@params$ext) {
      ggsave(
        filename   = paste0("network_contribution.", ext),
        plot       = patchwork_plot,
        device     = ext,
        path       = file.path(output_root, agg_lvl, "plot", fao_division),
        create.dir = TRUE,
        width      = 2480 * 2.2,
        height     = 1240 * 1.5 * scale_rows,
        units      = "px",
        dpi        = 300,
        bg         = "white"
      )
    }
  }
}