library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")
library("patchwork")
library("scales")
library("tibble")

# ── Shared palette and labels ────────────────────────────────────────────────

# Palette keyed on the COUNT columns (used by composition loop)
pal <- c(
  "nb_main_imp" = "#3D85F7",
  "nb_balanced" = "#A49393",
  "nb_main_exp" = "#C32E5A"
)

type_labels <- c(
  "nb_main_imp" = "Main importers",
  "nb_balanced" = "Balanced",
  "nb_main_exp" = "Main exporters"
)

# Stacking order: main_exp adjacent to spine, then balanced, main_imp at tip
# (symmetric on both sides of the mirror)
stack_order <- c("nb_main_imp", "nb_balanced", "nb_main_exp")

# Palette keyed on the VOLUME columns (used by mirrored desc stat loop)
# Same colours, different key names
pal_vol <- c(
  "main_imp" = "#3D85F7",
  "balanced" = "#A49393",
  "main_exp" = "#C32E5A"
)

vol_labels <- c(
  "main_imp" = "Main importers",
  "balanced" = "Balanced",
  "main_exp" = "Main exporters"
)

# Stacking order for volume (spine → tip): main_exp, balanced, main_imp
# Applied symmetrically on both sides
vol_stack_order <- c("main_exp", "balanced", "main_imp")

benchmark_years <- c(2000, 2005, 2010, 2015, 2020)

# ── Helper: auto-scale a volume vector ──────────────────────────────────────
# Returns a list(divisor, unit_label)
auto_scale <- function(values, unit = "") {
  mx  <- max(abs(values), na.rm = TRUE)
  sfx <- if (nchar(unit) > 0) paste0(" ", unit) else ""
  if (mx >= 1e9) {
    list(divisor = 1e9, label = paste0("billions", sfx))
  } else if (mx >= 1e6) {
    list(divisor = 1e6, label = paste0("millions", sfx))
  } else if (mx >= 1e3) {
    list(divisor = 1e3, label = paste0("thousands", sfx))
  } else {
    list(divisor = 1,   label = unit)
  }
}

# ── Helper: resolve vertical collisions for end-of-line labels ──────────────
# Same strategy as plot_network_contribution.R: greedy upward push (minimum
# gap = min_gap_frac of y_max), then shift the whole stack down if it
# overflows y_max (and up if it drops below 0).
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

# ── Helper: build one side of the mirrored bar chart ────────────────────────
# side      : "src" or "tgt" (attribution: exporter's or importer's category)
# direction : "left" (values negated) or "right"
# metric    : "net_wgt" or "primary_value"
# report    : "exp" or "imp" (whose declaration measures the flow — mirror
#             reports are kept separate, never averaged)
# divisor   : scaling factor
# min_seg   : minimum scaled value to show in-bar label
build_mirror_side <- function(dat, side, direction, metric, report,
                              divisor, min_seg) {

  cols <- paste0(side, "_", vol_stack_order, "_", metric, "_", report)

  # For the left side, reverse stacking order so main_exp is adjacent to spine.
  # position_stack stacks from axis outward in factor level order, so reversing
  # for the left (negative) side achieves the symmetric layout.
  stack_levels <- if (direction == "left") rev(vol_stack_order) else vol_stack_order # nolint

  dat_side <- dat %>% # nolint
    select(period, all_of(cols)) %>% # nolint
    setNames(c("period", vol_stack_order)) %>% # nolint
    melt(id.vars      = "period",
         variable.name = "trader_type",
         value.name    = "raw_vol") %>% # nolint
    mutate(
      trader_type = factor(trader_type, levels = stack_levels), # nolint
      vol         = raw_vol / divisor # nolint
    ) %>% # nolint
    group_by(period) %>% # nolint
    mutate(
      total_vol = sum(vol), # nolint
      share     = vol / total_vol * 100 # nolint
    ) %>% # nolint
    ungroup() %>% # nolint
    mutate(
      plot_vol  = if (direction == "left") -vol else vol,
      # Single-line label (value + share) — two-line labels overflow the bar
      # slot vertically and collide with neighbouring bars
      bar_label = ifelse(
        vol >= min_seg,
        paste0(round(vol, 1), " (", round(share, 0), "%)"), # nolint
        ""
      )
    )

  dat_side
}

# ── Helper: asymmetric x limits with data-driven tip padding ────────────────
# Each side of the mirror scales to its OWN maximum (the spine stays at 0,
# possibly off-centre), and is padded by just the room its widest total tip
# label needs (char_px per character plus a small gap), converted to data
# units through the final px-per-unit. Closed form: the axis range R solves
#   R = (M_l + M_r) / (1 - (w_l + w_r) / panel_px)
# with M_s the side maxima and w_s the tip-label widths in px. This replaces
# the former symmetric limits (global max × 1.12 × 1.22), which reserved
# ~27% of each half as padding and squeezed the smaller side's bars.
mirror_x_limits <- function(tot_left, tot_right, panel_px = 1620,
                            char_px = 15, gap_px = 12) {
  m_l <- max(abs(tot_left$tip_val))
  m_r <- max(abs(tot_right$tip_val))
  w_l <- max(nchar(tot_left$total_lab)) * char_px + gap_px
  w_r <- max(nchar(tot_right$total_lab)) * char_px + gap_px
  r   <- (m_l + m_r) / (1 - (w_l + w_r) / panel_px)
  list(left  = -(m_l + w_l * r / panel_px),
       right =   m_r + w_r * r / panel_px,
       range = r)
}

# Named input: "composition" = the per-agg_lvl network_composition.csv files
# (both loops)
composition_inputs <- snakemake@input[["composition"]]

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 1 — Network composition figure
# ══════════════════════════════════════════════════════════════════════════════

for (input_file in composition_inputs) {

  agg_lvl <- basename(dirname(dirname(input_file)))

  data <- read.csv(file   = input_file,
                   header = TRUE,
                   sep    = ";")

  for (fao_division in snakemake@params$fao_divisions) {

    prod <- as.numeric(fao_division)

    dat <- data %>% # nolint
      filter(cmd == prod) %>% # nolint
      arrange(period)

    dat_long <- dat %>% # nolint
      select(period, nb_main_exp, nb_main_imp, nb_balanced) %>% # nolint
      melt(id.vars      = "period",
           variable.name = "trader_type",
           value.name    = "nb_country") %>% # nolint
      mutate(trader_type = factor(trader_type, levels = stack_order))

    dat_shares <- dat_long %>% # nolint
      group_by(period) %>% # nolint
      mutate(
        share     = nb_country / sum(nb_country) * 100,
        bar_label = ifelse(share >= 5, paste0(round(share, 0), "%"), "")
      ) %>% # nolint
      ungroup()

    dat_total  <- dat %>% select(period, tot_nb_nodes) # nolint
    last_year  <- max(dat$period)
    last_vals  <- dat_long %>% filter(period == last_year) # nolint
    last_total <- dat_total %>% filter(period == last_year) # nolint

    # Data-driven y ceiling: round up to next multiple of 50 above the max
    y_max_nodes <- ceiling(max(dat$tot_nb_nodes, na.rm = TRUE) / 50) * 50
    y_breaks_nodes <- seq(0, y_max_nodes, by = 50)

    # End-of-line labels — same collision strategy as plot_network_contribution:
    # nudge overlapping labels apart, connect them to the series end with a
    # short segment. All four labels (3 trader types + total) nudged jointly.
    end_vals <- c(
      setNames(last_vals$nb_country, as.character(last_vals$trader_type)),
      total = last_total$tot_nb_nodes
    )
    lab_pal  <- c(pal, total = "black")
    lab_txt  <- c(type_labels, total = "Total")
    y_nudged <- compute_label_positions(end_vals, y_max_nodes,
                                        min_gap_frac = 0.04)
    label_df <- tibble(
      key   = names(end_vals),
      y_end = as.numeric(end_vals),
      y_lab = as.numeric(y_nudged[names(end_vals)]),
      label = unname(lab_txt[names(end_vals)]),
      col   = unname(lab_pal[names(end_vals)])
    )

    # ── (a) Line chart ───────────────────────────────────────────────────────
    p_line <- ggplot() +
      geom_vline(xintercept = benchmark_years,
                 color = "grey80", linewidth = 0.35, linetype = "dashed") +
      geom_line(data = dat_long,
                aes(x = period, y = nb_country,
                    color = trader_type, group = trader_type),
                linewidth = 0.8) +
      geom_line(data = dat_total,
                aes(x = period, y = tot_nb_nodes),
                color = "black", linewidth = 0.9) +
      # Connector segments from series end to nudged label position
      geom_segment(data = label_df,
                   aes(x = last_year + 0.15, xend = last_year + 0.35,
                       y = y_end, yend = y_lab),
                   color = label_df$col, linewidth = 0.35,
                   inherit.aes = FALSE) +
      geom_text(data = label_df,
                aes(x = last_year + 0.4, y = y_lab, label = label),
                color = label_df$col, hjust = 0, size = 2.8,
                fontface = "bold", inherit.aes = FALSE) +
      scale_color_manual(values = pal, labels = type_labels) +
      scale_x_continuous(
        breaks = c(min(dat$period), benchmark_years, last_year),
        expand = expansion(mult = c(0.02, 0.18))
      ) +
      scale_y_continuous(limits = c(0, y_max_nodes),
                         breaks = y_breaks_nodes,
                         expand = c(0, 0)) +
      labs(x = "Year",
           y = "Number of trading countries",
           title = "(a) Country count by trader type") +
      coord_cartesian(clip = "off") +
      theme_ipsum(axis_title_size = 10, base_size = 9) +
      theme(legend.position = "none",
            plot.margin     = margin(10, 60, 10, 10),
            plot.title      = element_text(size = 10, face = "bold"))

    # ── (b) 100% stacked bar ─────────────────────────────────────────────────
    p_bar <- ggplot(dat_shares,
                    aes(x = factor(period), y = share, fill = trader_type)) +
      geom_col(width = 0.85, alpha = 0.85) +
      geom_text(aes(label = bar_label),
                position = position_stack(vjust = 0.5),
                size = 2.2, color = "white", fontface = "bold") +
      scale_fill_manual(values = pal, labels = type_labels, name = "Trader type") + # nolint
      scale_x_discrete(
        breaks = as.character(c(min(dat$period), benchmark_years, last_year))
      ) +
      scale_y_continuous(expand = c(0, 0),
                         labels = function(x) paste0(x, "%")) +
      labs(x = "Year",
           y = "Share of trading countries",
           title = "(b) Composition by trader type (%)") +
      theme_ipsum(axis_title_size = 10, base_size = 9) +
      theme(legend.position = "right",
            legend.title    = element_text(size = 8),
            legend.text     = element_text(size = 8),
            plot.title      = element_text(size = 10, face = "bold"),
            plot.margin     = margin(10, 10, 10, 10))

    composite <- p_line + p_bar +
      plot_layout(ncol = 2, widths = c(1, 1.4))

    for (ext in snakemake@params$ext) {
      out_path <- file.path(
        "results", "network_analysis", agg_lvl, "plot",
        fao_division,
        paste0("network_composition.", ext)
      )
      ggsave(filename = out_path, plot = composite, device = ext,
             create.dir = TRUE, width = 4960, height = 1860,
             units = "px", dpi = 300, bg = "white")
    }
  } # end fao_division loop
} # end agg_lvl loop

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 2 — Mirrored descriptive statistics figure
# ══════════════════════════════════════════════════════════════════════════════

for (input_file in composition_inputs) {

  agg_lvl <- basename(dirname(dirname(input_file)))

  data <- read.csv(file   = input_file,
                   header = TRUE,
                   sep    = ";")

  for (fao_division in snakemake@params$fao_divisions) {

    prod <- as.numeric(fao_division)

    dat <- data %>% # nolint
      filter(cmd == prod) %>% # nolint
      arrange(period)

    # ── Auto-scale units per commodity ───────────────────────────────────────
    # net_wgt columns are in kg — convert to tonnes before scaling
    # Scaling and label thresholds span both mirror reports (exp and imp)
    wgt_cols <- grep("net_wgt", names(dat), value = TRUE)
    dat_t <- dat
    dat_t[wgt_cols] <- dat[wgt_cols] / 1000

    wgt_scale <- auto_scale(c(dat_t$tot_net_wgt_exp, dat_t$tot_net_wgt_imp),
                            unit = "tonnes")
    val_scale <- auto_scale(c(dat_t$tot_primary_value_exp,
                              dat_t$tot_primary_value_imp),
                            unit = "USD")

    # Minimum segment size to show in-bar label (8% of max total — slivers
    # are unlabelled; their share is the complement of the labelled ones)
    wgt_min_seg <- max(c(dat_t$tot_net_wgt_exp, dat_t$tot_net_wgt_imp) / wgt_scale$divisor, na.rm = TRUE) * 0.08 # nolint
    val_min_seg <- max(c(dat_t$tot_primary_value_exp, dat_t$tot_primary_value_imp) / val_scale$divisor, na.rm = TRUE) * 0.08 # nolint

    # ── Build long-format data for all four sides ─────────────────────────────
    # Natural pairing: each side of the mirror is drawn from that side's OWN
    # reports — supply (source attribution) from exporter reports (FOB for
    # value), demand (target attribution) from importer reports (CIF). The
    # left/right total asymmetry is therefore the aggregate mirror-report gap.
    src_val <- build_mirror_side(dat_t, "src", "left",  "primary_value", "exp",
                                 val_scale$divisor, val_min_seg)
    tgt_val <- build_mirror_side(dat_t, "tgt", "right", "primary_value", "imp",
                                 val_scale$divisor, val_min_seg)
    src_wgt <- build_mirror_side(dat_t, "src", "left",  "net_wgt", "exp",
                                 wgt_scale$divisor, wgt_min_seg)
    tgt_wgt <- build_mirror_side(dat_t, "tgt", "right", "net_wgt", "imp",
                                 wgt_scale$divisor, wgt_min_seg)

    # ── Total labels at bar tips ──────────────────────────────────────────────
    # For each period, one label at the outermost end of each bar
    make_totals <- function(dat_side, direction) {
      dat_side %>% # nolint
        group_by(period) %>% # nolint
        summarise(
          tip_val   = if (direction == "left") -sum(vol) else sum(vol),
          total_lab = as.character(round(sum(vol), 1)),
          .groups   = "drop"
        )
    }

    tot_src_val <- make_totals(src_val, "left")
    tot_tgt_val <- make_totals(tgt_val, "right")
    tot_src_wgt <- make_totals(src_wgt, "left")
    tot_tgt_wgt <- make_totals(tgt_wgt, "right")

    # ── Asymmetric x limits per metric (data-driven tip padding) ─────────────
    val_xlim <- mirror_x_limits(tot_src_val, tot_tgt_val)
    wgt_xlim <- mirror_x_limits(tot_src_wgt, tot_tgt_wgt)

    # ── Suppress in-bar labels that would overflow their segment ─────────────
    # Single-line labels are wide: keep one only if its estimated rendered
    # width fits inside the segment. The axis range maps to roughly 1620 px
    # of panel width at the 1880 px page size; one character at label size
    # 2.6 is ~15 px.
    fit_bar_labels <- function(dat_side, range_units) {
      char_unit <- range_units / 1620 * 15
      dat_side %>% # nolint
        mutate(bar_label = ifelse(
          vol >= nchar(bar_label) * char_unit * 1.1, # nolint
          bar_label,
          ""
        ))
    }

    src_val <- fit_bar_labels(src_val, val_xlim$range)
    tgt_val <- fit_bar_labels(tgt_val, val_xlim$range)
    src_wgt <- fit_bar_labels(src_wgt, wgt_xlim$range)
    tgt_wgt <- fit_bar_labels(tgt_wgt, wgt_xlim$range)

    # ── Mirror bar builder ───────────────────────────────────────────────────
    # Each call builds a STANDALONE full-page figure (one per metric), so
    # there is no panel letter and each figure carries its own legend.
    # left_note / right_note: side annotations naming the report side feeding
    # each half of the mirror (FOB/CIF for value); title_note carries the same
    # information on the title's second line. shade: draw the HS-revision
    # discontinuity band (net weight figure only).
    build_mirror_plot <- function(dat_left, dat_right,
                                  tot_left, tot_right,
                                  x_limits, unit_label, metric_name,
                                  title_note,
                                  left_note, right_note,
                                  shade = FALSE) {

      # Asymmetric limits from mirror_x_limits(): each side scales to its
      # own maximum; the spine stays at 0 (possibly off-centre)
      x_left  <- x_limits$left
      x_right <- x_limits$right

      # Year → position on the discrete y-axis (levels sorted ascending)
      yr_pos  <- function(y) match(y, sort(unique(dat$period)))
      n_years <- length(unique(dat$period))

      # Grid positions — same x breaks as the axis so lines and labels align.
      # Layer order (bottom → top): shaded band < grid < bars, so the theme
      # grid is blanked and redrawn as layers between band and bars.
      x_breaks <- scales::breaks_extended()(c(x_left, x_right))
      x_breaks <- x_breaks[x_breaks >= x_left & x_breaks <= x_right]

      ggplot() +

        # HS-revision discontinuity band — horizontal here because years sit
        # on the y-axis. Same windows as plot_market_concentration.R /
        # plot_network_contribution.R:
        # division 07: disruption 1996–1999 (+ transient 2007 line);
        # divisions 01/05: disruption 2000–2006
        (if (shade && fao_division == "07") annotate("rect",
          xmin = x_left, xmax = x_right,
          ymin = yr_pos(1996) - 0.5, ymax = yr_pos(1999) + 0.5,
          fill = "grey85"
        ) else if (shade) annotate("rect",
          xmin = x_left, xmax = x_right,
          ymin = yr_pos(2000) - 0.5, ymax = yr_pos(2006) + 0.5,
          fill = "grey85"
        ) else NULL) +

        # Grid redrawn above the band (theme grid blanked below)
        geom_hline(yintercept = seq_len(n_years),
                   color = "#cccccc", linewidth = 0.2) +
        geom_vline(xintercept = x_breaks,
                   color = "#cccccc", linewidth = 0.2) +

        # 2007 transient spike (division 07 net weight only)
        (if (shade && fao_division == "07") geom_hline(
          yintercept = yr_pos(2007), color = "grey70",
          linetype   = "solid", linewidth = 0.5
        ) else NULL) +

        # Zero spine
        geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +

        # LEFT bars — separate geom_col so factor levels (and stacking order)
        # are taken from dat_left independently of dat_right
        geom_col(data  = dat_left,
                 aes(y = factor(period), x = plot_vol, fill = trader_type),
                 width = 0.75, alpha = 0.85) +
        geom_text(data     = dat_left,
                  aes(y    = factor(period), x = plot_vol,
                      label = bar_label, fill = trader_type),
                  position  = position_stack(vjust = 0.5),
                  size      = 2.6, color = "white", fontface = "bold") +

        # RIGHT bars
        geom_col(data  = dat_right,
                 aes(y = factor(period), x = plot_vol, fill = trader_type),
                 width = 0.75, alpha = 0.85) +
        geom_text(data     = dat_right,
                  aes(y    = factor(period), x = plot_vol,
                      label = bar_label, fill = trader_type),
                  position  = position_stack(vjust = 0.5),
                  size      = 2.6, color = "white", fontface = "bold") +

        # Total labels at bar tips — the axis limits reserve just enough
        # room beyond each side's longest bar (see mirror_x_limits)
        geom_text(data = tot_left,
                  aes(y = factor(period), x = tip_val, label = total_lab),
                  hjust = 1.1, size = 2.7, color = "grey20") +
        geom_text(data = tot_right,
                  aes(y = factor(period), x = tip_val, label = total_lab),
                  hjust = -0.1, size = 2.7, color = "grey20") +

        # Side annotations — centered on each half's own midpoint
        annotate("text", x = x_left * 0.5, y = nlevels(factor(dat$period)) + 0.7, # nolint
                 label = left_note, vjust = -0.1,
                 size = 3.2, fontface = "italic", color = "grey30") +
        annotate("text", x = x_right * 0.5, y = nlevels(factor(dat$period)) + 0.7, # nolint
                 label = right_note, vjust = -0.1,
                 size = 3.2, fontface = "italic", color = "grey30") +

        scale_fill_manual(values = pal_vol,
                          labels = vol_labels,
                          name   = "Trader type") +
        scale_x_continuous(
          limits = c(x_left, x_right),
          breaks = x_breaks,
          labels = function(x) abs(x),
          expand = c(0, 0)
        ) +
        scale_y_discrete(
          breaks = as.character(c(min(dat$period), benchmark_years,
                                  max(dat$period)))
        ) +
        labs(y     = "Year",
             x     = paste0(metric_name, " (", unit_label, ")"),
             title = paste0(metric_name, " by trader type (", unit_label,
                            ")\n", title_note)) +
        coord_cartesian(clip = "off") +
        theme_ipsum(axis_title_size = 11, base_size = 10) +
        theme(legend.position  = "bottom",
              legend.title     = element_text(size = 9),
              legend.text      = element_text(size = 9),
              plot.title       = element_text(size = 11, face = "bold"),
              plot.margin      = margin(20, 10, 10, 10),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }
    # ── Build the two standalone figures ─────────────────────────────────────
    p_value <- build_mirror_plot(
      dat_left     = src_val,
      dat_right    = tgt_val,
      tot_left     = tot_src_val,
      tot_right    = tot_tgt_val,
      x_limits     = val_xlim,
      unit_label   = val_scale$label,
      metric_name  = "Primary value",
      title_note   = "exporter reports (FOB) vs importer reports (CIF)",
      left_note    = "Supply (source) — FOB",
      right_note   = "Demand (target) — CIF"
    )

    p_weight <- build_mirror_plot(
      dat_left     = src_wgt,
      dat_right    = tgt_wgt,
      tot_left     = tot_src_wgt,
      tot_right    = tot_tgt_wgt,
      x_limits     = wgt_xlim,
      unit_label   = wgt_scale$label,
      metric_name  = "Net weight",
      title_note   = "exporter vs importer reports",
      left_note    = "Supply (source)",
      right_note   = "Demand (target)",
      shade        = TRUE
    )

    # ── Save — one full-page portrait figure per metric ──────────────────────
    # Each metric gets its own page-size figure so the 28 bar rows keep
    # enough vertical space for readable in-bar labels. Portrait A4 inside
    # normal 2.54 cm margins at 300 dpi:
    # (21.0 - 2*2.54) x (29.7 - 2*2.54) cm ≈ 1880 x 2900 px.
    mirror_figs <- list(value = p_value, weight = p_weight)
    for (metric_key in names(mirror_figs)) {
      for (ext in snakemake@params$ext) {
        out_path <- file.path(
          "results", "network_analysis", agg_lvl, "plot",
          fao_division,
          paste0("network_mirrored_desc_stat_", metric_key, ".", ext)
        )
        ggsave(filename = out_path, plot = mirror_figs[[metric_key]],
               device = ext, create.dir = TRUE, width = 1880, height = 2900,
               units = "px", dpi = 300, bg = "white")
      }
    }
  } # end fao_division loop
} # end agg_lvl loop
