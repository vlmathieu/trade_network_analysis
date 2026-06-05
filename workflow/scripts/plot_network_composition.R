library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")
library("patchwork")
library("scales")
library("tibble")
library("ggrepel")

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

# ── Helper: build one side of the mirrored bar chart ────────────────────────
# side      : "src" or "tgt"
# direction : "left" (values negated) or "right"
# metric    : "net_wgt" or "primary_value"
# divisor   : scaling factor
# min_seg   : minimum scaled value to show in-bar label
build_mirror_side <- function(dat, side, direction, metric, divisor, min_seg) {

  cols <- paste0(side, "_", vol_stack_order, "_", metric)

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
      bar_label = ifelse(
        vol >= min_seg,
        paste0(round(vol, 1), "\n(", round(share, 0), "%)"), # nolint
        ""
      )
    )

  dat_side
}

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 1 — Network composition figure
# ══════════════════════════════════════════════════════════════════════════════

for (input_file in snakemake@input) {

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
      geom_text(data = last_vals,
                aes(x = last_year + 0.3, y = nb_country,
                    label = type_labels[as.character(trader_type)],
                    color = trader_type),
                hjust = 0, size = 2.8, fontface = "bold") +
      geom_text(data = last_total,
                aes(x = last_year + 0.3, y = tot_nb_nodes, label = "Total"),
                color = "black", hjust = 0, size = 2.8, fontface = "bold") +
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

for (input_file in snakemake@input) {

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
    wgt_cols <- grep("net_wgt", names(dat), value = TRUE)
    dat_t <- dat
    dat_t[wgt_cols] <- dat[wgt_cols] / 1000

    wgt_scale <- auto_scale(dat_t$tot_net_wgt, unit = "tonnes")
    val_scale <- auto_scale(dat_t$tot_primary_value, unit = "USD")

    # Minimum segment size to show in-bar label (5% of max total)
    wgt_min_seg <- max(dat_t$tot_net_wgt / wgt_scale$divisor, na.rm = TRUE) * 0.05 # nolint
    val_min_seg <- max(dat_t$tot_primary_value / val_scale$divisor, na.rm = TRUE) * 0.05 # nolint

    # ── Build long-format data for all four sides ─────────────────────────────
    src_val <- build_mirror_side(dat_t, "src", "left",  "primary_value",
                                 val_scale$divisor, val_min_seg)
    tgt_val <- build_mirror_side(dat_t, "tgt", "right", "primary_value",
                                 val_scale$divisor, val_min_seg)
    src_wgt <- build_mirror_side(dat_t, "src", "left",  "net_wgt",
                                 wgt_scale$divisor, wgt_min_seg)
    tgt_wgt <- build_mirror_side(dat_t, "tgt", "right", "net_wgt",
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

    # ── Shared x limit for each metric (symmetric around 0) ──────────────────
    wgt_lim <- max(abs(c(tot_src_wgt$tip_val, tot_tgt_wgt$tip_val))) * 1.12
    val_lim <- max(abs(c(tot_src_val$tip_val, tot_tgt_val$tip_val))) * 1.12

    # ── Mirror bar builder ───────────────────────────────────────────────────
    build_mirror_plot <- function(dat_left, dat_right,
                                  tot_left, tot_right,
                                  x_lim, unit_label, metric_name,
                                  panel_letter) {

      # Extra room beyond x_lim for total-label text (20% each side)
      x_plot_lim <- x_lim * 1.22

      ggplot() +

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
                  size      = 1.9, color = "white", fontface = "bold",
                  lineheight = 0.85) +

        # RIGHT bars
        geom_col(data  = dat_right,
                 aes(y = factor(period), x = plot_vol, fill = trader_type),
                 width = 0.75, alpha = 0.85) +
        geom_text(data     = dat_right,
                  aes(y    = factor(period), x = plot_vol,
                      label = bar_label, fill = trader_type),
                  position  = position_stack(vjust = 0.5),
                  size      = 1.9, color = "white", fontface = "bold",
                  lineheight = 0.85) +

        # Total labels at bar tips — placed outside x_lim, inside x_plot_lim
        geom_text(data = tot_left,
                  aes(y = factor(period), x = tip_val, label = total_lab),
                  hjust = 1.15, size = 2.2, color = "grey20") +
        geom_text(data = tot_right,
                  aes(y = factor(period), x = tip_val, label = total_lab),
                  hjust = -0.15, size = 2.2, color = "grey20") +

        # Side annotations
        annotate("text", x = -x_lim * 0.5, y = nlevels(factor(dat$period)) + 0.7, # nolint
                 label = "Supply (source)", vjust = -0.1,
                 size = 2.8, fontface = "italic", color = "grey30") +
        annotate("text", x =  x_lim * 0.5, y = nlevels(factor(dat$period)) + 0.7, # nolint
                 label = "Demand (target)", vjust = -0.1,
                 size = 2.8, fontface = "italic", color = "grey30") +

        scale_fill_manual(values = pal_vol,
                          labels = vol_labels,
                          name   = "Trader type") +
        scale_x_continuous(
          limits = c(-x_plot_lim, x_plot_lim),
          labels = function(x) abs(x),
          expand = c(0, 0)
        ) +
        scale_y_discrete(
          breaks = as.character(c(min(dat$period), benchmark_years,
                                  max(dat$period)))
        ) +
        labs(y     = "Year",
             x     = paste0(metric_name, " (", unit_label, ")"),
             title = paste0(panel_letter, " ", metric_name,
                            " by trader type (", unit_label, ")")) +
        coord_cartesian(clip = "off") +
        theme_ipsum(axis_title_size = 10, base_size = 9) +
        theme(legend.position = "right",
              legend.title    = element_text(size = 8),
              legend.text     = element_text(size = 8),
              plot.title      = element_text(size = 10, face = "bold"),
              plot.margin     = margin(30, 10, 10, 10))
    }
    # ── Build panels ─────────────────────────────────────────────────────────
    p_value <- build_mirror_plot(
      dat_left     = src_val,
      dat_right    = tgt_val,
      tot_left     = tot_src_val,
      tot_right    = tot_tgt_val,
      x_lim        = val_lim,
      unit_label   = val_scale$label,
      metric_name  = "Primary value",
      panel_letter = "(a)"
    )

    p_weight <- build_mirror_plot(
      dat_left     = src_wgt,
      dat_right    = tgt_wgt,
      tot_left     = tot_src_wgt,
      tot_right    = tot_tgt_wgt,
      x_lim        = wgt_lim,
      unit_label   = wgt_scale$label,
      metric_name  = "Net weight",
      panel_letter = "(b)"
    )


    # ── Composite: 2-panel landscape ─────────────────────────────────────────
    composite_desc <- p_value + p_weight +
      plot_layout(ncol = 2, guides = "collect") &
      theme(legend.position = "bottom",
            plot.margin     = margin(30, 5, 10, 5))

    # ── Save ─────────────────────────────────────────────────────────────────
    # Landscape A4 at 300 dpi: 3508 x 2480 px
    for (ext in snakemake@params$ext) {
      out_path <- file.path(
        "results", "network_analysis", agg_lvl, "plot",
        fao_division,
        paste0("network_mirrored_desc_stat.", ext)
      )
      ggsave(filename = out_path, plot = composite_desc, device = ext,
             create.dir = TRUE, width = 3508, height = 2480,
             units = "px", dpi = 300, bg = "white")
    }
  } # end fao_division loop
} # end agg_lvl loop

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 3 — Price descriptive statistics figure
# Two-panel composite: (a) supply (source) perspective, (b) demand (target)
# Each panel: one line per trader type, years on x-axis, $/tonne on y-axis
# ══════════════════════════════════════════════════════════════════════════════

for (input_file in snakemake@input) {

  agg_lvl <- basename(dirname(dirname(input_file)))

  data <- read.csv(file   = input_file,
                   header = TRUE,
                   sep    = ";")

  for (fao_division in snakemake@params$fao_divisions) {

    prod <- as.numeric(fao_division)

    dat <- data %>% # nolint
      filter(cmd == prod) %>% # nolint
      arrange(period)

    # ── Compute price (USD/tonne) — net_wgt in kg, divide by 1000 for tonnes ──
    # All trader types shown, including balanced.
    # Total average price added: tot_primary_value / tot_net_wgt (kg->t).
    price_dat <- tibble(period = dat$period) %>% # nolint
      mutate(
        src_main_exp = dat$src_main_exp_primary_value / (dat$src_main_exp_net_wgt / 1000), # nolint
        src_main_imp = dat$src_main_imp_primary_value / (dat$src_main_imp_net_wgt / 1000), # nolint
        src_balanced = dat$src_balanced_primary_value / (dat$src_balanced_net_wgt / 1000), # nolint
        src_total    = dat$tot_primary_value           / (dat$tot_net_wgt           / 1000), # nolint
        tgt_main_exp = dat$tgt_main_exp_primary_value / (dat$tgt_main_exp_net_wgt / 1000), # nolint
        tgt_main_imp = dat$tgt_main_imp_primary_value / (dat$tgt_main_imp_net_wgt / 1000), # nolint
        tgt_balanced = dat$tgt_balanced_primary_value / (dat$tgt_balanced_net_wgt / 1000), # nolint
        tgt_total    = dat$tot_primary_value           / (dat$tot_net_wgt           / 1000) # nolint
      )

    # Palette and labels: all trader types + total
    pal_price    <- c(pal_vol["main_exp"], pal_vol["balanced"], pal_vol["main_imp"],
                      total = "black")
    price_labels <- c(vol_labels["main_exp"], vol_labels["balanced"],
                      vol_labels["main_imp"], total = "Total")
    price_levels <- c("main_exp", "balanced", "main_imp", "total")

    # Long format per perspective — all trader types + total
    price_src <- price_dat %>% # nolint
      select(period, src_main_exp, src_balanced, src_main_imp, src_total) %>% # nolint
      setNames(c("period", "main_exp", "balanced", "main_imp", "total")) %>% # nolint
      melt(id.vars = "period", variable.name = "trader_type",
           value.name = "price") %>% # nolint
      mutate(trader_type = factor(trader_type, levels = price_levels))

    price_tgt <- price_dat %>% # nolint
      select(period, tgt_main_exp, tgt_balanced, tgt_main_imp, tgt_total) %>% # nolint
      setNames(c("period", "main_exp", "balanced", "main_imp", "total")) %>% # nolint
      melt(id.vars = "period", variable.name = "trader_type",
           value.name = "price") %>% # nolint
      mutate(trader_type = factor(trader_type, levels = price_levels))

    # Shared y-axis limit across both panels for direct comparison
    y_max <- max(c(price_src$price, price_tgt$price), na.rm = TRUE) * 1.06

    # ── Helper: build one price line panel ───────────────────────────────────
    build_price_panel <- function(price_long, panel_letter, panel_title) {

      last_year <- max(price_long$period)
      last_vals <- price_long %>% filter(period == last_year) # nolint

      ggplot(price_long,
             aes(x = period, y = price,
                 color = trader_type, group = trader_type)) +

        geom_vline(xintercept = benchmark_years,
                   color = "grey80", linewidth = 0.35, linetype = "dashed") +

        geom_line(linewidth = 0.85) +
        geom_point(size = 1.8, shape = 16) +

        # End-of-series labels — repel vertically to avoid overlap
        ggrepel::geom_text_repel(
                  data      = last_vals, # nolint
                  aes(x     = last_year + 0.3, # nolint
                      label = price_labels[as.character(trader_type)]), # nolint
                  hjust        = 0, # nolint
                  direction    = "y",
                  nudge_x      = 0.2,
                  segment.size = 0.3,
                  segment.color = "grey60",
                  box.padding  = 0.3,
                  size         = 2.8,
                  fontface     = "bold",
                  force        = 2) +

        scale_color_manual(values = pal_price, labels = price_labels,
                           name = "Trader type") +
        scale_x_continuous(
          breaks = c(min(price_long$period), benchmark_years, last_year),
          expand = expansion(mult = c(0.02, 0.20))
        ) +
        scale_y_continuous(
          limits = c(0, y_max),
          expand = c(0, 0),
          labels = scales::comma
        ) +
        coord_cartesian(clip = "off") +
        labs(x     = "Year",
             y     = "Price ($/tonne)",
             title = paste0(panel_letter, " ", panel_title)) +
        theme_ipsum(axis_title_size = 10, base_size = 9) +
        theme(legend.position = "none",
              plot.margin     = margin(10, 70, 10, 10),
              plot.title      = element_text(size = 10, face = "bold"))
    }

    p_src <- build_price_panel(price_src, "(a)", "Average price — supply (source) perspective") # nolint
    p_tgt <- build_price_panel(price_tgt, "(b)", "Average price — demand (target) perspective") # nolint

    composite_price <- p_src + p_tgt +
      plot_layout(ncol = 2)

    # ── Save ─────────────────────────────────────────────────────────────────
    # Landscape — same width as composition figure for visual consistency
    for (ext in snakemake@params$ext) {
      out_path <- file.path(
        "results", "network_analysis", agg_lvl, "plot",
        fao_division,
        paste0("network_price_desc_stat.", ext)
      )
      ggsave(filename = out_path, plot = composite_price, device = ext,
             create.dir = TRUE, width = 4960, height = 1860,
             units = "px", dpi = 300, bg = "white")
    }
  } # end fao_division loop
} # end agg_lvl loop