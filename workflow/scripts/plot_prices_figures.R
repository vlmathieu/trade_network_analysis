library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")
library("patchwork")
library("scales")
library("tibble")
library("ggrepel")

# ── Shared palette and labels ────────────────────────────────────────────────
# Duplicated from plot_network_composition.R — both scripts share the
# trader-type colour scheme and the benchmark years

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

benchmark_years <- c(2000, 2005, 2010, 2015, 2020)

# Named inputs: "composition" = the per-agg_lvl network_composition.csv files
# (LOOP 1 iterates over them; LOOP 2 reads the country_lvl one for the
# main-importer classification); "mirror_flows" = flow-level mirror_flows.csv
# (LOOP 2 only)
composition_inputs <- snakemake@input[["composition"]]

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 1 — Price descriptive statistics figure
# Two-panel composite: (a) supply (source) perspective, (b) demand (target)
# Each panel: one line per trader type, years on x-axis, $/tonne on y-axis
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

    # ── Compute price (USD/tonne) — net_wgt in kg, divide by 1000 for tonnes ──
    # All trader types shown, including balanced.
    # Natural pairing of attribution and report: supply-side prices from
    # exporter reports (FOB), demand-side prices from importer reports (CIF).
    # Mirror reports are never mixed within a ratio.
    # The numerator uses the *_primary_value_wp_* columns (value summed only over
    # flows whose same-side net weight is reported), so numerator and denominator
    # cover one flow set and the ratio is not inflated where net weight is missing
    # over the 2000-2006 gap (see network_composition.py; sec-netweight-gaps).
    # Total average price added per side: tot_primary_value_wp / tot_net_wgt.
    price_dat <- tibble(period = dat$period) %>% # nolint
      mutate(
        src_main_exp = dat$src_main_exp_primary_value_wp_exp / (dat$src_main_exp_net_wgt_exp / 1000), # nolint
        src_main_imp = dat$src_main_imp_primary_value_wp_exp / (dat$src_main_imp_net_wgt_exp / 1000), # nolint
        src_balanced = dat$src_balanced_primary_value_wp_exp / (dat$src_balanced_net_wgt_exp / 1000), # nolint
        src_total    = dat$tot_primary_value_wp_exp           / (dat$tot_net_wgt_exp           / 1000), # nolint
        tgt_main_exp = dat$tgt_main_exp_primary_value_wp_imp / (dat$tgt_main_exp_net_wgt_imp / 1000), # nolint
        tgt_main_imp = dat$tgt_main_imp_primary_value_wp_imp / (dat$tgt_main_imp_net_wgt_imp / 1000), # nolint
        tgt_balanced = dat$tgt_balanced_primary_value_wp_imp / (dat$tgt_balanced_net_wgt_imp / 1000), # nolint
        tgt_total    = dat$tot_primary_value_wp_imp           / (dat$tot_net_wgt_imp           / 1000) # nolint
      )

    # Palette and labels: all trader types + total
    pal_price    <- c(pal_vol["main_exp"], pal_vol["balanced"], pal_vol["main_imp"], # nolint
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
    # Layer order (bottom → top): shaded band < grid < curves. The theme grid
    # is blanked and redrawn as layers so it sits ABOVE the band but BELOW the
    # curves (panel.ontop would put it above the curves as well).
    build_price_panel <- function(price_long, panel_letter, panel_title) {

      last_year <- max(price_long$period)
      last_vals <- price_long %>% filter(period == last_year) # nolint

      # Grid positions — same breaks as the axes so lines and labels align
      y_breaks <- scales::breaks_extended()(c(0, y_max))
      y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_max]
      x_breaks <- c(min(price_long$period), benchmark_years, last_year)

      ggplot(price_long,
             aes(x = period, y = price,
                 color = trader_type, group = trader_type)) +

        # HS-revision discontinuity band — prices divide value by net weight,
        # so the net-weight reporting gap inflates BOTH panels; same windows
        # as plot_market_concentration.R / plot_network_contribution.R:
        # division 07: disruption 1996–1999 (+ transient 2007 spike);
        # divisions 01/05: disruption 2000–2006
        (if (fao_division == "07") annotate("rect",
          xmin = 1996, xmax = 1999, ymin = 0, ymax = y_max,
          fill = "grey85"
        ) else annotate("rect",
          xmin = 2000, xmax = 2006, ymin = 0, ymax = y_max,
          fill = "grey85"
        )) +

        # Grid redrawn above the band (theme grid blanked below)
        geom_hline(yintercept = y_breaks, color = "#cccccc", linewidth = 0.25) + # nolint
        geom_vline(xintercept = x_breaks, color = "#cccccc", linewidth = 0.25) + # nolint

        # 2007 transient spike (division 07 only)
        (if (fao_division == "07") geom_vline(
          xintercept = 2007, color = "grey70",
          linetype = "solid", linewidth = 0.5
        ) else NULL) +

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
          breaks = y_breaks,
          expand = c(0, 0),
          labels = scales::comma
        ) +
        coord_cartesian(clip = "off") +
        labs(x     = "Year",
             y     = "Price (US$/tonne)",
             title = paste0(panel_letter, " ", panel_title)) +
        theme_ipsum(axis_title_size = 10, base_size = 9) +
        theme(legend.position  = "none",
              plot.margin      = margin(10, 70, 10, 10),
              plot.title       = element_text(size = 10, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }

    p_src <- build_price_panel(price_src, "(a)", "Average price — supply (source) perspective, exporter reports (FOB)") # nolint
    p_tgt <- build_price_panel(price_tgt, "(b)", "Average price — demand (target) perspective, importer reports (CIF)") # nolint

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

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 2 — Price distribution figures (country_lvl, division 01 only)
# FIVE STANDALONE figures:
# network_price_partners: spaghetti small multiples (3 x 2) of the average
#     import price of the top-6 main importers — each panel highlights one
#     country against the main-importers group price (same series as the
#     network_price_desc_stat figure) and the five other selected countries
#     in light grey; dotted segments bridge years with no computable price.
# network_price_ridgeline_countries_{imp,exp}: all-years-pooled price
#     distribution per main trader (top 10 by traded value + the pooled
#     main-importers / main-exporters group), ridgeline style with nested
#     quantile strips, median dots and an inset legend.
# network_price_ridgeline_years_{imp,exp}: the transposed companions —
#     all-countries pooled price distribution per YEAR (chronological
#     rows), same ridgeline construction, sharing the per-side inset
#     legend with the countries figure.
# The _imp figures read importer reports (CIF, blue quantile strips); the
# _exp figures read exporter reports (FOB, red strips). Own reports on each
# side, mirror sides never mixed within a ratio. Country level: EU
# aggregation would fold the euro-dollar exchange-rate component into
# prices.
# ══════════════════════════════════════════════════════════════════════════════

flows <- read.csv(file       = snakemake@input[["mirror_flows"]],
                  header     = TRUE,
                  sep        = ";",
                  colClasses = c(fao_code_agg = "character"),
                  na.strings = c("", "NA"))

# Flow-level import unit prices: a price exists only where the importer
# reported both value and net weight strictly positive (a zero/missing
# report yields no price and the flow is excluded, not imputed)
flows01 <- flows %>% # nolint
  filter(fao_code_agg == "01",
         !is.na(primary_value_imp), primary_value_imp > 0, # nolint
         !is.na(net_wgt_imp),       net_wgt_imp > 0) %>% # nolint
  mutate(price = primary_value_imp / (net_wgt_imp / 1000))

# Flow-level export unit prices — same rule on the exporter-reported (FOB)
# value and net weight; feeds the _exp ridgeline figures. Display alias:
# "Russian Federation" -> "Russia" (as in plot_contributor_profiles.R)
flows01_exp <- flows %>% # nolint
  filter(fao_code_agg == "01",
         !is.na(primary_value_exp), primary_value_exp > 0, # nolint
         !is.na(net_wgt_exp),       net_wgt_exp > 0) %>% # nolint
  mutate(price = primary_value_exp / (net_wgt_exp / 1000),
         exporter_desc = ifelse(exporter_desc == "Russian Federation",
                                "Russia", exporter_desc))

yrs_dist    <- sort(unique(flows01$period))
x_breaks_yr <- c(min(yrs_dist), benchmark_years, max(yrs_dist))

# ── Figure 1: import price spaghetti — top-6 main importers ─────────────────
# Spaghetti small multiples (data-to-viz.com/caveat/spaghetti.html, 3 x 2):
# each panel highlights one country's price curve against the main-importers
# group reference and the five other selected countries in light grey.
# Selection: top-6 importers by CIF value pooled over the period (China,
# Japan, India, Austria, Rep. of Korea, Sweden) — each classified a main
# importer in most years. Each country-year price is the ratio of sums over
# that country's import flows, i.e. its value-weighted average import price
top6_importers <- flows01 %>% # nolint
  group_by(importer_desc) %>% # nolint
  summarise(tot_value = sum(primary_value_imp), .groups = "drop") %>% # nolint
  arrange(desc(tot_value)) %>% # nolint
  slice_head(n = 6) %>% # nolint
  pull(importer_desc)

partner_prices <- flows01 %>% # nolint
  filter(importer_desc %in% top6_importers) %>% # nolint
  group_by(importer_desc, period) %>% # nolint
  summarise(price = sum(primary_value_imp) / sum(net_wgt_imp / 1000),
            .groups = "drop")

# Country-level composition rows for division 01 — read once: the partners
# figure takes the main-importers group price from it, the ridgeline
# figures the per-year trader-type classification (list_main_imp /
# list_main_exp)
comp_file <- grep("country_lvl", composition_inputs, value = TRUE)
comp01 <- read.csv(comp_file, header = TRUE, sep = ";") %>% # nolint
  filter(cmd == 1)

# Main-importers group price (target attribution, importer reports) —
# the same series as panel (b) of the network_price_desc_stat figure,
# for direct comparison
group_prices <- comp01 %>% # nolint
  transmute(
    period, # nolint
    importer_desc = "Main importers",
    price = tgt_main_imp_primary_value_imp / (tgt_main_imp_net_wgt_imp / 1000) # nolint
  ) %>% # nolint
  filter(is.finite(price))

# Complete series x year grid: years with no computable price become NA
# rows, which BREAK the solid line there (notably China/Japan/India over
# 2000-2006, whose importer-side net weights are entirely zeroed)
partner_full <- expand.grid(importer_desc    = top6_importers,
                            period           = yrs_dist,
                            stringsAsFactors = FALSE) %>% # nolint
  left_join(partner_prices, by = c("importer_desc", "period"))

# Panel scaffolding: one facet per selected country, ordered by pooled
# import value ((a) China ... (f) Sweden). The grey background of a panel
# is the five OTHER countries; the blue reference and the shaded band
# replicate in every panel
panel_labels <- setNames(
  sprintf("(%s) %s", letters[seq_along(top6_importers)], top6_importers),
  top6_importers
)
as_panel <- function(country) {
  factor(panel_labels[country], levels = unname(panel_labels))
}
panels_df <- data.frame(panel_country = top6_importers,
                        stringsAsFactors = FALSE)

# Context spaghetti: complete grid (solid lines break at NA) plus a
# sparse frame feeding the dotted bridge over the missing years, same
# convention as the highlighted series
spag_grey <- merge(panels_df, partner_full, by = NULL) %>% # nolint
  filter(importer_desc != panel_country) %>% # nolint
  mutate(panel = as_panel(panel_country), role = "other")
spag_grey_sparse <- merge(panels_df, partner_prices, by = NULL) %>% # nolint
  filter(importer_desc != panel_country) %>% # nolint
  mutate(panel = as_panel(panel_country), role = "other")

# Main-importers group reference, replicated in every panel
spag_ref <- merge(panels_df, group_prices, by = NULL) %>% # nolint
  mutate(panel = as_panel(panel_country), role = "ref")

# Highlighted series: sparse frame (points + dotted bridge) and complete
# grid (solid line broken over the missing years)
spag_hl <- partner_prices %>% # nolint
  mutate(panel = as_panel(importer_desc), role = "sel")
spag_hl_full <- partner_full %>% # nolint
  mutate(panel = as_panel(importer_desc), role = "sel")

# Roles, not countries, carry the colour: the facet strip names the
# highlighted country. Selected importer = dark blue (imports stay in the
# blue family; red is reserved for the export side)
spag_pal <- c(sel = "#16418F", ref = "#3D85F7", other = "grey70")
spag_lab <- c(sel   = "Selected importer",
              ref   = "Main importers (reference)",
              other = "Other top-6 importers")

# Log scale here too: partner prices in artefact country-years (partial net
# weights) spike orders of magnitude above the clean band, and showing them
# as-is is the point — they mark reporting errors, not real prices
y_min_b    <- min(c(partner_prices$price, group_prices$price)) * 0.7
y_max_b    <- max(c(partner_prices$price, group_prices$price)) * 1.4
y_breaks_b <- 10^seq(ceiling(log10(y_min_b)), floor(log10(y_max_b)))

p_partner <- ggplot(mapping = aes(x = period, y = price,
                                  color = role,
                                  group = importer_desc)) +

  # Shaded band drawn first so the redrawn grid sits above it (annotate
  # replicates in every facet)
  annotate("rect", xmin = 2000, xmax = 2006,
           ymin = y_min_b, ymax = y_max_b, fill = "grey92") +

  # Grid redrawn above the band (theme grid blanked below)
  geom_hline(yintercept = y_breaks_b, color = "#cccccc", linewidth = 0.25) +
  geom_vline(xintercept = x_breaks_yr, color = "#cccccc", linewidth = 0.25) +

  geom_vline(xintercept = benchmark_years,
             color = "grey80", linewidth = 0.35, linetype = "dashed") +

  # Context spaghetti underneath (dotted bridge first, solid lines over
  # it), then the blue reference, then the highlighted country on top
  geom_line(data = spag_grey_sparse, linetype = "dotted",
            linewidth = 0.35) +
  geom_line(data = spag_grey, linewidth = 0.35, na.rm = TRUE) +
  geom_line(data = spag_ref, linewidth = 0.6) +

  # Highlighted country: dotted bridge connecting observed points straight
  # across missing years (drawn first; the solid line overplots it wherever
  # data exist — the NA rows of the complete grid break it over the years
  # with no computable price, na.rm only silences the warning)
  geom_line(data = spag_hl, linetype = "dotted", linewidth = 0.5) +
  geom_line(data = spag_hl_full, linewidth = 0.8, na.rm = TRUE) +
  geom_point(data = spag_hl, size = 1.1, shape = 16) +

  facet_wrap(~panel, ncol = 2) +

  scale_color_manual(values = spag_pal,
                     labels = spag_lab,
                     breaks = c("sel", "ref", "other")) +
  scale_x_continuous(breaks = x_breaks_yr,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_log10(limits = c(y_min_b, y_max_b),
                breaks = y_breaks_b,
                labels = scales::comma,
                expand = c(0, 0)) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.1))) +
  labs(x     = "Year",
       y     = "Import price (US$/tonne, log scale)",
       color = NULL) +
  theme_ipsum(axis_title_size = 10, base_size = 9) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 8),
        plot.margin      = margin(10, 15, 10, 10),
        strip.text       = element_text(size = 9, face = "bold", hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ── Figure 2: ridgeline of pooled price distributions by main importer ──────
# Adapted from r-graph-gallery web-ridgeline-plot-with-inside-plot-and-
# annotations (ggdist stat_halfeye + stat_interval + median point),
# reproduced with core ggplot2 + patchwork only: per-group density ribbon
# on log10 price, nested 50/80/95% quantile strips (blues instead of the
# tutorial's greens), black median dot, dotted reference line at the
# main-importers median, an "n = x" column between the country labels and
# the distributions (instead of the tutorial's avg-bedroom counts) and an
# inset legend explaining the encoding. All years pooled — the temporal
# view is network_price_ridgeline_years_imp.

top10_importers <- flows01 %>% # nolint
  group_by(importer_desc) %>% # nolint
  summarise(tot_value = sum(primary_value_imp), .groups = "drop") %>% # nolint
  arrange(desc(tot_value)) %>% # nolint
  slice_head(n = 10) %>% # nolint
  pull(importer_desc)

# Pooled "Main importers" reference group: flows whose importer is
# classified main importer in that flow's year (list_main_imp is
# pipe-separated — same convention as plot_network_contribution.R)
main_imp_long <- do.call(rbind, Map(function(p, l) {
  data.frame(period           = p,
             importer_desc    = strsplit(l, "|", fixed = TRUE)[[1]],
             stringsAsFactors = FALSE)
}, comp01$period, comp01$list_main_imp))

ridge_flows <- bind_rows(
  flows01 %>% # nolint
    filter(importer_desc %in% top10_importers) %>% # nolint
    transmute(grp = importer_desc, lp = log10(price)), # nolint
  flows01 %>% # nolint
    inner_join(main_imp_long, by = c("period", "importer_desc")) %>% # nolint
    transmute(grp = "Main importers", lp = log10(price)) # nolint
)

# One row per group, reading top -> bottom from smallest to highest
# median (row 1 = highest median = bottom line, so the smallest-median
# group sits on the first line); quantile pairs feed the interval strips
grp_stats <- ridge_flows %>% # nolint
  group_by(grp) %>% # nolint
  summarise(n    = n(),
            med  = median(lp), # nolint
            q025 = quantile(lp, 0.025), q10 = quantile(lp, 0.10), # nolint
            q25  = quantile(lp, 0.25),  q75 = quantile(lp, 0.75), # nolint
            q90  = quantile(lp, 0.90),  q975 = quantile(lp, 0.975), # nolint
            .groups = "drop") %>% # nolint
  arrange(desc(med)) %>% # nolint
  mutate(row = row_number())

n_rows <- nrow(grp_stats)

# Per-group density on log10 price, normalised to a common max height
# (as stat_halfeye does per group)
ridge_height <- 0.85
dens_df <- do.call(rbind, lapply(seq_len(n_rows), function(i) {
  lp <- ridge_flows$lp[ridge_flows$grp == grp_stats$grp[i]]
  d  <- density(lp, from = min(lp), to = max(lp))
  data.frame(row = i, x = d$x,
             ymax = i + d$y / max(d$y) * ridge_height)
}))

med_mi   <- grp_stats$med[grp_stats$grp == "Main importers"]
int_cols <- c(p95 = "#CBDCF9", p80 = "#8FB3F0", p50 = "#3D85F7")

# Horizontal range: pooled 0.1-99.9 percentiles (the extreme tails carry
# no visible density). The n column sits LEFT of the ridges, between the
# country labels and the distributions (as in the tutorial); it is
# right-aligned ending at x_lab_n, well short of the first x gridline
# (10 US$/t, kept — informative on a log scale)
x_lo    <- as.numeric(quantile(ridge_flows$lp, 0.001))
x_hi    <- as.numeric(quantile(ridge_flows$lp, 0.999))
x_brk   <- seq(ceiling(x_lo), floor(x_hi))
# n column ends just short of the first gridline (10 US$/t), and the
# panel's left edge sits just left of the column — tight on both sides
x_lab_n <- x_brk[1] - 0.15

# Trim the ridge tails to the plotted range so they cannot run under the
# n column (their density out there is ~0 anyway)
dens_df <- dens_df %>% filter(x >= x_lo - 0.02, x <= x_hi + 0.02) # nolint

p_ridge <- ggplot() +

  # Layer order (bottom -> top): density ridges < quantile strips <
  # median dots < dotted reference line (no horizontal grid)

  # Density ridges — grey, fill only (no outline)
  geom_ribbon(data = dens_df,
              aes(x = x, ymin = row, ymax = ymax, group = row),
              fill = "grey75", alpha = 0.5) +

  # Nested quantile strips on the baseline: 95% < 80% < 50% (light -> dark)
  geom_segment(data = grp_stats,
               aes(x = q025, xend = q975, y = row, yend = row),
               color = int_cols[["p95"]], linewidth = 2.4) +
  geom_segment(data = grp_stats,
               aes(x = q10, xend = q90, y = row, yend = row),
               color = int_cols[["p80"]], linewidth = 2.4) +
  geom_segment(data = grp_stats,
               aes(x = q25, xend = q75, y = row, yend = row),
               color = int_cols[["p50"]], linewidth = 2.4) +

  # Median dot — black, fill only, small enough to stay inside the bar
  geom_point(data = grp_stats, aes(x = med, y = row),
             shape = 16, color = "black", size = 1.3) +

  # Dotted reference line on top of everything: median price of the
  # pooled main-importers group
  geom_vline(xintercept = med_mi, linetype = "dotted",
             color = "grey30", linewidth = 0.45) +
  # Annotation right of the line, on the same level as the n-column header
  annotate("text", x = med_mi + 0.06, y = n_rows + 1.05,
           label = sprintf("Main importers median: US$%s/t",
                           scales::comma(round(10^med_mi))),
           hjust = 0, size = 2.8, color = "grey30") +

  # Number of price observations (the tutorial's avg nb of bedrooms) —
  # right-aligned column between the country labels and the ridges
  geom_text(data = grp_stats,
            aes(x = x_lab_n, y = row,
                label = paste0("n = ", scales::comma(n))),
            hjust = 1, size = 2.5, color = "grey40") +
  annotate("text", x = x_lab_n, y = n_rows + 1.05,
           label = "Nb. of price obs.", hjust = 1,
           size = 2.8, color = "grey30",
           lineheight = 0.95) +

  scale_x_continuous(breaks = x_brk,
                     labels = scales::comma(10^x_brk)) +
  scale_y_continuous(breaks = grp_stats$row,
                     labels = grp_stats$grp,
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(x_lab_n - 0.42, x_hi + 0.25),
                  ylim = c(0.55, n_rows + 1.9)) +
  labs(x = "Import unit price, CIF (US$/tonne, log scale)", y = NULL) +
  theme_ipsum(axis_title_size = 10, base_size = 9) +
  theme(legend.position    = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(10, 10, 10, 10))

# Inset legend (the tutorial's "inside plot"): one pooled example
# distribution with the encoding annotated. Its x-range is tightened to
# the bar + label room (NOT the main panel's full range), so the
# distribution and bar fill the box; labels follow the tutorial layout —
# "Median" just above the dot, "Distribution of prices" at the right of
# the density at the same level, the three interval labels on one shared
# horizontal line below the bar
leg_d  <- density(ridge_flows$lp, from = x_lo, to = x_hi)
leg_df <- data.frame(x = leg_d$x, ymax = leg_d$y / max(leg_d$y))
# Cut the distribution's tails: the negative crop margins push the panel
# edges outside the drawn frame, so the near-zero density queues rendered
# past the border. Keep only the contiguous range where the density
# reaches at least 5% of its peak — the quantile strips still mark the
# full 95% spread on the baseline.
leg_rng <- range(leg_df$x[leg_df$ymax >= 0.05])
leg_df  <- leg_df %>% filter(x >= leg_rng[1], x <= leg_rng[2]) # nolint
leg_q  <- quantile(ridge_flows$lp,
                   c(0.025, 0.10, 0.25, 0.5, 0.75, 0.90, 0.975))
# Height of the density curve near a given x (arrow target helper)
leg_y_at <- function(x0) leg_df$ymax[which.min(abs(leg_df$x - x0))]

# Label anchors, all on the same line below the bar: 50% pushed left of
# its strip, its arrow targeting the point 25% ALONG the dark bar
# (between the q25 frontier and the median, inside the dark blue); 80%
# under its own visible strip segment; 95% pushed right
leg_lab_y <- -0.42
lab50_x   <- leg_q[["25%"]] - 0.45
lab50_tgt <- leg_q[["25%"]] + 0.25 * (leg_q[["75%"]] - leg_q[["25%"]])
lab80_x   <- (leg_q[["75%"]] + leg_q[["90%"]]) / 2 - 0.5
lab95_x   <- leg_q[["97.5%"]] + 0.45 - 0.80
# Rightmost x where the density still reaches 0.35 of its peak — the
# "Distribution of prices" arrow lands on the curve there
leg_x35   <- max(leg_df$x[leg_df$ymax >= 0.8])

p_ridge_legend <- ggplot() +
  # Same encoding as the main panel: grey fill-only density
  geom_ribbon(data = leg_df, aes(x = x, ymin = 0, ymax = ymax),
              fill = "grey75", alpha = 0.5) +
  geom_segment(aes(x = leg_q[["2.5%"]], xend = leg_q[["97.5%"]],
                   y = 0, yend = 0),
               color = int_cols[["p95"]], linewidth = 3.4) +
  geom_segment(aes(x = leg_q[["10%"]], xend = leg_q[["90%"]],
                   y = 0, yend = 0),
               color = int_cols[["p80"]], linewidth = 3.4) +
  geom_segment(aes(x = leg_q[["25%"]], xend = leg_q[["75%"]],
                   y = 0, yend = 0),
               color = int_cols[["p50"]], linewidth = 3.4) +
  geom_point(aes(x = leg_q[["50%"]], y = 0),
             shape = 16, color = "black", size = 1.6) +

  # "Median" just above the median dot, short arrow down to it
  annotate("text", x = leg_q[["50%"]], y = 0.32, label = "Median",
           hjust = 0.5, size = 2.1, color = "grey20") +
  annotate("segment", x = leg_q[["50%"]], y = 0.23,
           xend = leg_q[["50%"]], yend = 0.10,
           linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +

  # "Distribution of prices" at the right of the density, same vertical
  # level as the Median label, arrow onto the curve's right slope
  annotate("text", x = leg_x35 + 0.85, y = 0.75,
           label = "Distribution\nof prices",
           hjust = 0.5, size = 2.1, color = "grey20") +
  annotate("curve", x = leg_x35 + 0.3, y = 0.75,
           xend = leg_x35 + 0.12, yend = 0.66,
           curvature = 0.15, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +

  # Interval labels — one shared horizontal line below the bar, arrows up
  # to the VISIBLE part of each strip (the portion no darker strip covers)
  annotate("text", x = lab50_x, y = leg_lab_y,
           label = "50% of prices\nfall within this range",
           hjust = 0.5, size = 2.1, color = "grey20") +
  annotate("curve", x = lab50_x + 0.15, y = leg_lab_y + 0.18,
           xend = lab50_tgt, yend = -0.12,
           curvature = -0.1, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +
  annotate("text", x = lab80_x, y = leg_lab_y, label = "80% of prices",
           hjust = 0, size = 2.1, color = "grey20") +
  annotate("curve", x = lab80_x + 0.50, y = leg_lab_y + 0.10,
           xend = (leg_q[["75%"]] + leg_q[["90%"]]) / 2, yend = -0.12,
           curvature = 0.2, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +
  annotate("text", x = lab95_x, y = leg_lab_y, label = "95% of prices",
           hjust = 0, size = 2.1, color = "grey20") +
  annotate("curve", x = lab95_x + 0.50, y = leg_lab_y + 0.10,
           xend = (leg_q[["90%"]] + leg_q[["97.5%"]]) / 2, yend = -0.12,
           curvature = 0.2, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +

  # Legend title — bold, top-left corner of the box (as in the tutorial),
  # nudged down and right off the frame edge
  annotate("text", x = x_lo + 0.45, y = 1.00, label = "Legend",
           hjust = 0, vjust = 1, size = 2.3, fontface = "bold",
           color = "grey20") +

  # x range data-driven and tight: just enough for the pushed-right 95%
  # label and the "Distribution of prices" label, no dead space.
  # clip = "on": anything crossing the panel boundary is cut instead of
  # spilling past the legend frame
  coord_cartesian(xlim = c(x_lo - 0.05,
                           max(lab95_x + 0.95, leg_x35 + 1.95)),
                  ylim = c(-0.82, 1.12), clip = "on") +
  theme_void() +
  # Square border around the legend box (transparent inside, so the grid
  # and ridge tails behind stay visible); slim margins keep the frame
  # close to the content
  theme(plot.margin     = margin(-10, -10, -20, -20),
        plot.background = element_rect(color = "grey60", fill = "white",
                                       linewidth = 0.4))

# Embed at the right of Austria's distribution: the box's vertical centre
# sits on Austria's row baseline. inset_element fractions are panel-relative
# (its default align_to), so the centre is converted from data coords
# using the main panel's ylim = c(0.55, n_rows + 1.9)
austria_row <- grp_stats$row[grp_stats$grp == "Austria"]
if (length(austria_row) == 0) austria_row <- ceiling(n_rows / 2)
leg_ctr <- (austria_row - 0.55) / ((n_rows + 1.9) - 0.55)
# Flatter but wider box, flush against the panel's right edge (right = 1)
leg_hgt <- 0.26
p_ridgeline_countries_imp <- p_ridge +
  inset_element(p_ridge_legend,
                left   = 0.72,
                bottom = max(0.02, leg_ctr - leg_hgt / 2),
                right  = 1,
                top    = min(0.98, leg_ctr + leg_hgt / 2),
                clip   = FALSE)

# ── Figure 3: ridgeline of price distributions by year ──────────────────────
# The transposed companion of Figure 2: rows are the 28 trade years, each
# pooling ALL countries' flow-level import prices for that year. Same
# construction (density ridge, nested 50/80/95% strips, median dot, n
# column, inset legend). Chronological order — earliest year on the top
# line, latest at the bottom — NOT median order: the point of this figure
# is the time reading. The dotted reference line is the overall pooled
# median across all years, and the 2000-2006 HS-revision net-weight break
# gets its usual grey band behind the affected rows (years sit on an axis
# again here, unlike the pooled Figure 2). Variables suffixed _ry;
# x_breaks_yr (calendar-year axis breaks of the partners figure) is
# unrelated.

ry_flows <- flows01 %>% # nolint
  transmute(period, lp = log10(price)) # nolint

# row 1 = latest year = bottom line, so reading top -> bottom runs
# forward in time
ry_stats <- ry_flows %>% # nolint
  group_by(period) %>% # nolint
  summarise(n    = n(),
            med  = median(lp), # nolint
            q025 = quantile(lp, 0.025), q10 = quantile(lp, 0.10), # nolint
            q25  = quantile(lp, 0.25),  q75 = quantile(lp, 0.75), # nolint
            q90  = quantile(lp, 0.90),  q975 = quantile(lp, 0.975), # nolint
            .groups = "drop") %>% # nolint
  arrange(desc(period)) %>% # nolint
  mutate(row = row_number())

ry_n_rows <- nrow(ry_stats)

# Per-year density on log10 price, same common-max normalisation and
# ridge height as Figure 2
ry_dens <- do.call(rbind, lapply(seq_len(ry_n_rows), function(i) {
  lp <- ry_flows$lp[ry_flows$period == ry_stats$period[i]]
  d  <- density(lp, from = min(lp), to = max(lp))
  data.frame(row = i, x = d$x,
             ymax = i + d$y / max(d$y) * ridge_height)
}))

# Reference line: pooled median over ALL flows and years (the analogue of
# Figure 2's main-importers median)
ry_med_all <- median(ry_flows$lp)

# Horizontal range, axis breaks and n-column anchor — same recipe as
# Figure 2, on this figure's own pooled distribution
ry_x_lo    <- as.numeric(quantile(ry_flows$lp, 0.001))
ry_x_hi    <- as.numeric(quantile(ry_flows$lp, 0.999))
ry_x_brk   <- seq(ceiling(ry_x_lo), floor(ry_x_hi))
ry_x_lab_n <- ry_x_brk[1] - 0.15

# Trim the ridge tails to the plotted range (no visible density out there)
ry_dens <- ry_dens %>% filter(x >= ry_x_lo - 0.02, x <= ry_x_hi + 0.02) # nolint

# Rows hit by the HS-revision net-weight break (2000-2006, division 01)
ry_break_rows <- ry_stats$row[ry_stats$period %in% 2000:2006]

p_ridge_ry <- ggplot() +

  # Layer order (bottom -> top): shaded band < grid < density ridges <
  # quantile strips < median dots < dotted reference line — the full
  # shade < grid < geoms convention shared with plot_network_composition.R (theme grid blanked,
  # vertical gridlines redrawn as a layer above the band)

  # HS-revision band behind the affected rows. Edges midway between the
  # in-band "n = x" labels and their neighbours on both sides (+/- 0.5);
  # left border starts at the n column's slot rather than the panel edge,
  # so the band stays clear of the year labels
  annotate("rect", xmin = ry_x_lab_n - 0.42, xmax = Inf,
           ymin = min(ry_break_rows) - 0.5,
           ymax = max(ry_break_rows) + 0.5,
           fill = "grey95") +

  # Grid redrawn above the band (theme grid blanked below)
  geom_vline(xintercept = ry_x_brk, color = "#cccccc", linewidth = 0.25) +

  # Density ridges — grey, fill only (no outline)
  geom_ribbon(data = ry_dens,
              aes(x = x, ymin = row, ymax = ymax, group = row),
              fill = "grey75", alpha = 0.5) +

  # Nested quantile strips on the baseline: 95% < 80% < 50% (light -> dark)
  geom_segment(data = ry_stats,
               aes(x = q025, xend = q975, y = row, yend = row),
               color = int_cols[["p95"]], linewidth = 2.4) +
  geom_segment(data = ry_stats,
               aes(x = q10, xend = q90, y = row, yend = row),
               color = int_cols[["p80"]], linewidth = 2.4) +
  geom_segment(data = ry_stats,
               aes(x = q25, xend = q75, y = row, yend = row),
               color = int_cols[["p50"]], linewidth = 2.4) +

  # Median dot — black, fill only, inside the bar
  geom_point(data = ry_stats, aes(x = med, y = row),
             shape = 16, color = "black", size = 1.3) +

  # Dotted reference line on top: pooled median across all years
  geom_vline(xintercept = ry_med_all, linetype = "dotted",
             color = "grey30", linewidth = 0.45) +
  annotate("text", x = ry_med_all + 0.06, y = ry_n_rows + 1.05,
           label = sprintf("All-years median: US$%s/t",
                           scales::comma(round(10^ry_med_all))),
           hjust = 0, size = 2.8, color = "grey30") +

  # Number of price observations per year — right-aligned column between
  # the year labels and the ridges
  geom_text(data = ry_stats,
            aes(x = ry_x_lab_n, y = row,
                label = paste0("n = ", scales::comma(n))),
            hjust = 1, size = 2.5, color = "grey40") +
  annotate("text", x = ry_x_lab_n, y = ry_n_rows + 1.05,
           label = "Nb. of price obs.", hjust = 1,
           size = 2.8, color = "grey30",
           lineheight = 0.95) +

  scale_x_continuous(breaks = ry_x_brk,
                     labels = scales::comma(10^ry_x_brk)) +
  scale_y_continuous(breaks = ry_stats$row,
                     labels = ry_stats$period,
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(ry_x_lab_n - 0.42, ry_x_hi + 0.25),
                  ylim = c(0.55, ry_n_rows + 1.9)) +
  labs(x = "Import unit price, CIF (US$/tonne, log scale)", y = NULL) +
  theme_ipsum(axis_title_size = 10, base_size = 9) +
  theme(legend.position    = "none",
        panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(10, 10, 10, 10))

# Inset legend: REUSED from Figure 2 (p_ridge_legend) — it explains the
# encoding, which is identical, and its pooled example distribution is
# near-identical in shape to this figure's. Height fraction scaled to the
# taller 2900 px page so the box keeps roughly the same absolute size;
# vertical centre on the 1997/1998 row midpoint (top of the figure, clear
# of the 2000-2006 band), converted from data coords through the panel's
# ylim = c(0.55, ry_n_rows + 1.9)
ry_leg_hgt  <- 0.20
ry_leg_rows <- ry_stats$row[ry_stats$period %in% c(1997, 1998)]
if (length(ry_leg_rows) == 0) ry_leg_rows <- ceiling(ry_n_rows / 2)
ry_leg_ctr  <- (mean(ry_leg_rows) - 0.55) / ((ry_n_rows + 1.9) - 0.55)
p_ridgeline_years_imp <- p_ridge_ry +
  inset_element(p_ridge_legend,
                left   = 0.72,
                bottom = max(0.02, ry_leg_ctr - ry_leg_hgt / 2),
                right  = 1,
                top    = min(0.98, ry_leg_ctr + ry_leg_hgt / 2),
                clip   = FALSE)

# ── Figure 4: ridgeline of pooled price distributions by main exporter ──────
# The exporter-side mirror of Figure 2, on exporter reports (FOB): same
# construction (density ridge, nested 50/80/95% quantile strips, median
# dot, n column, inset legend), red strips instead of blue (ladder anchored
# on the main-exporters palette red), top-10 exporters by exported value
# plus the pooled main-exporters group, dotted reference line at the
# main-exporters median. Variables suffixed _exp.

top10_exporters <- flows01_exp %>% # nolint
  group_by(exporter_desc) %>% # nolint
  summarise(tot_value = sum(primary_value_exp), .groups = "drop") %>% # nolint
  arrange(desc(tot_value)) %>% # nolint
  slice_head(n = 10) %>% # nolint
  pull(exporter_desc)

# Pooled "Main exporters" reference group: flows whose exporter is
# classified main exporter in that flow's year (list_main_exp is
# pipe-separated — same convention as list_main_imp above)
main_exp_long <- do.call(rbind, Map(function(p, l) {
  data.frame(period           = p,
             exporter_desc    = strsplit(l, "|", fixed = TRUE)[[1]],
             stringsAsFactors = FALSE)
}, comp01$period, comp01$list_main_exp))
# Same alias as flows01_exp so the group join still matches Russia's flows
main_exp_long$exporter_desc[
  main_exp_long$exporter_desc == "Russian Federation"] <- "Russia"

ridge_flows_exp <- bind_rows(
  flows01_exp %>% # nolint
    filter(exporter_desc %in% top10_exporters) %>% # nolint
    transmute(grp = exporter_desc, lp = log10(price)), # nolint
  flows01_exp %>% # nolint
    inner_join(main_exp_long, by = c("period", "exporter_desc")) %>% # nolint
    transmute(grp = "Main exporters", lp = log10(price)) # nolint
)

# One row per group, same reading order as Figure 2 (row 1 = highest
# median = bottom line)
grp_stats_exp <- ridge_flows_exp %>% # nolint
  group_by(grp) %>% # nolint
  summarise(n    = n(),
            med  = median(lp), # nolint
            q025 = quantile(lp, 0.025), q10 = quantile(lp, 0.10), # nolint
            q25  = quantile(lp, 0.25),  q75 = quantile(lp, 0.75), # nolint
            q90  = quantile(lp, 0.90),  q975 = quantile(lp, 0.975), # nolint
            .groups = "drop") %>% # nolint
  arrange(desc(med)) %>% # nolint
  mutate(row = row_number())

n_rows_exp <- nrow(grp_stats_exp)

dens_df_exp <- do.call(rbind, lapply(seq_len(n_rows_exp), function(i) {
  lp <- ridge_flows_exp$lp[ridge_flows_exp$grp == grp_stats_exp$grp[i]]
  d  <- density(lp, from = min(lp), to = max(lp))
  data.frame(row = i, x = d$x,
             ymax = i + d$y / max(d$y) * ridge_height)
}))

med_me       <- grp_stats_exp$med[grp_stats_exp$grp == "Main exporters"]
# Red ladder — same lightness steps as the blue int_cols, anchored on the
# shared main-exporters red (#C32E5A)
int_cols_exp <- c(p95 = "#F5CFDA", p80 = "#E289A3", p50 = "#C32E5A")

# Horizontal range, axis breaks and n-column anchor — same recipe as
# Figure 2, on the exporter-side pooled distribution
x_lo_exp    <- as.numeric(quantile(ridge_flows_exp$lp, 0.001))
x_hi_exp    <- as.numeric(quantile(ridge_flows_exp$lp, 0.999))
x_brk_exp   <- seq(ceiling(x_lo_exp), floor(x_hi_exp))
x_lab_n_exp <- x_brk_exp[1] - 0.15

dens_df_exp <- dens_df_exp %>% # nolint
  filter(x >= x_lo_exp - 0.02, x <= x_hi_exp + 0.02) # nolint

p_ridge_exp <- ggplot() +

  # Layer order (bottom -> top): density ridges < quantile strips <
  # median dots < dotted reference line (no horizontal grid)

  geom_ribbon(data = dens_df_exp,
              aes(x = x, ymin = row, ymax = ymax, group = row),
              fill = "grey75", alpha = 0.5) +

  geom_segment(data = grp_stats_exp,
               aes(x = q025, xend = q975, y = row, yend = row),
               color = int_cols_exp[["p95"]], linewidth = 2.4) +
  geom_segment(data = grp_stats_exp,
               aes(x = q10, xend = q90, y = row, yend = row),
               color = int_cols_exp[["p80"]], linewidth = 2.4) +
  geom_segment(data = grp_stats_exp,
               aes(x = q25, xend = q75, y = row, yend = row),
               color = int_cols_exp[["p50"]], linewidth = 2.4) +

  geom_point(data = grp_stats_exp, aes(x = med, y = row),
             shape = 16, color = "black", size = 1.3) +

  # Dotted reference line: median price of the pooled main-exporters group
  geom_vline(xintercept = med_me, linetype = "dotted",
             color = "grey30", linewidth = 0.45) +
  annotate("text", x = med_me + 0.06, y = n_rows_exp + 1.05,
           label = sprintf("Main exporters median: US$%s/t",
                           scales::comma(round(10^med_me))),
           hjust = 0, size = 2.8, color = "grey30") +

  geom_text(data = grp_stats_exp,
            aes(x = x_lab_n_exp, y = row,
                label = paste0("n = ", scales::comma(n))),
            hjust = 1, size = 2.5, color = "grey40") +
  annotate("text", x = x_lab_n_exp, y = n_rows_exp + 1.05,
           label = "Nb. of price obs.", hjust = 1,
           size = 2.8, color = "grey30",
           lineheight = 0.95) +

  scale_x_continuous(breaks = x_brk_exp,
                     labels = scales::comma(10^x_brk_exp)) +
  scale_y_continuous(breaks = grp_stats_exp$row,
                     labels = grp_stats_exp$grp,
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(x_lab_n_exp - 0.42, x_hi_exp + 0.25),
                  ylim = c(0.55, n_rows_exp + 1.9)) +
  labs(x = "Export unit price, FOB (US$/tonne, log scale)", y = NULL) +
  theme_ipsum(axis_title_size = 10, base_size = 9) +
  theme(legend.position    = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(10, 10, 10, 10))

# Inset legend for the _exp figures — same construction and hand-tuned
# offsets as the blue one, on the exporter-side pooled distribution and
# the red interval ladder
leg_d_exp  <- density(ridge_flows_exp$lp, from = x_lo_exp, to = x_hi_exp)
leg_df_exp <- data.frame(x = leg_d_exp$x, ymax = leg_d_exp$y / max(leg_d_exp$y)) # nolint
leg_rng_exp <- range(leg_df_exp$x[leg_df_exp$ymax >= 0.05])
leg_df_exp  <- leg_df_exp %>% # nolint
  filter(x >= leg_rng_exp[1], x <= leg_rng_exp[2]) # nolint
leg_q_exp  <- quantile(ridge_flows_exp$lp,
                       c(0.025, 0.10, 0.25, 0.5, 0.75, 0.90, 0.975))

leg_lab_y_exp <- -0.42
lab50_x_exp   <- leg_q_exp[["25%"]] - 0.45
lab50_tgt_exp <- leg_q_exp[["25%"]] +
  0.25 * (leg_q_exp[["75%"]] - leg_q_exp[["25%"]])
lab80_x_exp   <- (leg_q_exp[["75%"]] + leg_q_exp[["90%"]]) / 2 - 0.5
lab95_x_exp   <- leg_q_exp[["97.5%"]] + 0.45 - 0.80
leg_x35_exp   <- max(leg_df_exp$x[leg_df_exp$ymax >= 0.8])

p_ridge_legend_exp <- ggplot() +
  geom_ribbon(data = leg_df_exp, aes(x = x, ymin = 0, ymax = ymax),
              fill = "grey75", alpha = 0.5) +
  geom_segment(aes(x = leg_q_exp[["2.5%"]], xend = leg_q_exp[["97.5%"]],
                   y = 0, yend = 0),
               color = int_cols_exp[["p95"]], linewidth = 3.4) +
  geom_segment(aes(x = leg_q_exp[["10%"]], xend = leg_q_exp[["90%"]],
                   y = 0, yend = 0),
               color = int_cols_exp[["p80"]], linewidth = 3.4) +
  geom_segment(aes(x = leg_q_exp[["25%"]], xend = leg_q_exp[["75%"]],
                   y = 0, yend = 0),
               color = int_cols_exp[["p50"]], linewidth = 3.4) +
  geom_point(aes(x = leg_q_exp[["50%"]], y = 0),
             shape = 16, color = "black", size = 1.6) +

  annotate("text", x = leg_q_exp[["50%"]], y = 0.32, label = "Median",
           hjust = 0.5, size = 2.1, color = "grey20") +
  annotate("segment", x = leg_q_exp[["50%"]], y = 0.23,
           xend = leg_q_exp[["50%"]], yend = 0.10,
           linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +

  # Text pulled left, close to the arrow start (the wider box leaves it room)
  annotate("text", x = leg_x35_exp + 0.60, y = 0.75,
           label = "Distribution\nof prices",
           hjust = 0.5, size = 2.1, color = "grey20") +
  annotate("curve", x = leg_x35_exp + 0.3, y = 0.75,
           xend = leg_x35_exp + 0.12, yend = 0.66,
           curvature = 0.15, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +

  annotate("text", x = lab50_x_exp, y = leg_lab_y_exp,
           label = "50% of prices\nfall within this range",
           hjust = 0.5, size = 2.1, color = "grey20") +
  annotate("curve", x = lab50_x_exp + 0.15, y = leg_lab_y_exp + 0.18,
           xend = lab50_tgt_exp, yend = -0.12,
           curvature = -0.1, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +
  annotate("text", x = lab80_x_exp, y = leg_lab_y_exp,
           label = "80% of prices",
           hjust = 0, size = 2.1, color = "grey20") +
  annotate("curve", x = lab80_x_exp + 0.50, y = leg_lab_y_exp + 0.10,
           xend = (leg_q_exp[["75%"]] + leg_q_exp[["90%"]]) / 2, yend = -0.12,
           curvature = 0.2, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +
  annotate("text", x = lab95_x_exp, y = leg_lab_y_exp, label = "95% of prices",
           hjust = 0, size = 2.1, color = "grey20") +
  annotate("curve", x = lab95_x_exp + 0.50, y = leg_lab_y_exp + 0.10,
           xend = (leg_q_exp[["90%"]] + leg_q_exp[["97.5%"]]) / 2, yend = -0.12,
           curvature = 0.2, linewidth = 0.25, color = "grey40",
           arrow = arrow(length = unit(3, "pt"))) +

  # Title in the box's top-left corner (no left crop anymore — the panel
  # sits fully inside the frame, so a small pad suffices)
  annotate("text", x = x_lo_exp + 0.10, y = 1.00, label = "Legend",
           hjust = 0, vjust = 1, size = 2.3, fontface = "bold",
           color = "grey20") +

  # Right bound trimmed to just past the "95% of prices" label (the old
  # +0.95/+1.95 paddings, tuned on the imp shape and the pre-move text
  # position, left dead space right of the content on the exp distribution)
  coord_cartesian(xlim = c(x_lo_exp - 0.05,
                           max(lab95_x_exp + 0.80, leg_x35_exp + 1.30)),
                  ylim = c(-0.70, 1.08), clip = "on") +
  theme_void() +
  # UNLIKE the imp legend, no negative crop margins: those push the panel
  # edges outside the frame, so edge labels render past the border — and
  # the narrower the inset box, the more data-units the fixed crop eats.
  # Small positive padding + clip = "on" keeps everything inside the frame;
  # the ylim above is tightened accordingly (the crop no longer hides the
  # empty strips beyond the bottom labels and the title line)
  theme(plot.margin     = margin(2, 2, 2, 2),
        plot.background = element_rect(color = "grey60", fill = "white",
                                       linewidth = 0.4))

# Embed at the right of New Zealand's distribution (short right tail near
# the lower rows); fallback: middle row. Box narrowed together with the
# legend's trimmed internal x-range (left 0.63, was 0.60 — same ratio, so
# the drawn content keeps its scale), right edge flush
nz_row_exp <- grp_stats_exp$row[grp_stats_exp$grp == "New Zealand"]
if (length(nz_row_exp) == 0) nz_row_exp <- ceiling(n_rows_exp / 2)
leg_ctr_exp <- (nz_row_exp - 0.55) / ((n_rows_exp + 1.9) - 0.55)
p_ridgeline_countries_exp <- p_ridge_exp +
  inset_element(p_ridge_legend_exp,
                left   = 0.70,
                bottom = max(0.02, leg_ctr_exp - leg_hgt / 2),
                right  = 1,
                top    = min(0.98, leg_ctr_exp + leg_hgt / 2),
                clip   = FALSE)

# ── Figure 5: ridgeline of export price distributions by year ───────────────
# The transposed companion of Figure 4 and the exporter-side mirror of
# Figure 3: rows are the 28 trade years pooling ALL countries' flow-level
# export prices, chronological order, reference line at the pooled
# all-years median, 2000-2006 HS-revision band, red quantile strips,
# shared _exp inset legend. Variables suffixed _ry_exp.

ry_flows_exp <- flows01_exp %>% # nolint
  transmute(period, lp = log10(price)) # nolint

ry_stats_exp <- ry_flows_exp %>% # nolint
  group_by(period) %>% # nolint
  summarise(n    = n(),
            med  = median(lp), # nolint
            q025 = quantile(lp, 0.025), q10 = quantile(lp, 0.10), # nolint
            q25  = quantile(lp, 0.25),  q75 = quantile(lp, 0.75), # nolint
            q90  = quantile(lp, 0.90),  q975 = quantile(lp, 0.975), # nolint
            .groups = "drop") %>% # nolint
  arrange(desc(period)) %>% # nolint
  mutate(row = row_number())

ry_n_rows_exp <- nrow(ry_stats_exp)

ry_dens_exp <- do.call(rbind, lapply(seq_len(ry_n_rows_exp), function(i) {
  lp <- ry_flows_exp$lp[ry_flows_exp$period == ry_stats_exp$period[i]]
  d  <- density(lp, from = min(lp), to = max(lp))
  data.frame(row = i, x = d$x,
             ymax = i + d$y / max(d$y) * ridge_height)
}))

ry_med_all_exp <- median(ry_flows_exp$lp)

ry_x_lo_exp    <- as.numeric(quantile(ry_flows_exp$lp, 0.001))
ry_x_hi_exp    <- as.numeric(quantile(ry_flows_exp$lp, 0.999))
ry_x_brk_exp   <- seq(ceiling(ry_x_lo_exp), floor(ry_x_hi_exp))
ry_x_lab_n_exp <- ry_x_brk_exp[1] - 0.15

ry_dens_exp <- ry_dens_exp %>% # nolint
  filter(x >= ry_x_lo_exp - 0.02, x <= ry_x_hi_exp + 0.02) # nolint

# The 2000-2006 net-weight break zeroes BOTH mirror sides of the affected
# flows (bilateral signature), so the exporter-side figure carries the
# same band as the importer-side one
ry_break_rows_exp <- ry_stats_exp$row[ry_stats_exp$period %in% 2000:2006]

p_ridge_ry_exp <- ggplot() +

  # Layer order (bottom -> top): shaded band < grid < density ridges <
  # quantile strips < median dots < dotted reference line — the full
  # shade < grid < geoms convention shared with plot_network_composition.R

  annotate("rect", xmin = ry_x_lab_n_exp - 0.42, xmax = Inf,
           ymin = min(ry_break_rows_exp) - 0.5,
           ymax = max(ry_break_rows_exp) + 0.5,
           fill = "grey95") +

  geom_vline(xintercept = ry_x_brk_exp, color = "#cccccc", linewidth = 0.25) +

  geom_ribbon(data = ry_dens_exp,
              aes(x = x, ymin = row, ymax = ymax, group = row),
              fill = "grey75", alpha = 0.5) +

  geom_segment(data = ry_stats_exp,
               aes(x = q025, xend = q975, y = row, yend = row),
               color = int_cols_exp[["p95"]], linewidth = 2.4) +
  geom_segment(data = ry_stats_exp,
               aes(x = q10, xend = q90, y = row, yend = row),
               color = int_cols_exp[["p80"]], linewidth = 2.4) +
  geom_segment(data = ry_stats_exp,
               aes(x = q25, xend = q75, y = row, yend = row),
               color = int_cols_exp[["p50"]], linewidth = 2.4) +

  geom_point(data = ry_stats_exp, aes(x = med, y = row),
             shape = 16, color = "black", size = 1.3) +

  geom_vline(xintercept = ry_med_all_exp, linetype = "dotted",
             color = "grey30", linewidth = 0.45) +
  annotate("text", x = ry_med_all_exp + 0.06, y = ry_n_rows_exp + 1.05,
           label = sprintf("All-years median: US$%s/t",
                           scales::comma(round(10^ry_med_all_exp))),
           hjust = 0, size = 2.8, color = "grey30") +

  geom_text(data = ry_stats_exp,
            aes(x = ry_x_lab_n_exp, y = row,
                label = paste0("n = ", scales::comma(n))),
            hjust = 1, size = 2.5, color = "grey40") +
  annotate("text", x = ry_x_lab_n_exp, y = ry_n_rows_exp + 1.05,
           label = "Nb. of price obs.", hjust = 1,
           size = 2.8, color = "grey30",
           lineheight = 0.95) +

  scale_x_continuous(breaks = ry_x_brk_exp,
                     labels = scales::comma(10^ry_x_brk_exp)) +
  scale_y_continuous(breaks = ry_stats_exp$row,
                     labels = ry_stats_exp$period,
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(ry_x_lab_n_exp - 0.42, ry_x_hi_exp + 0.25),
                  ylim = c(0.55, ry_n_rows_exp + 1.9)) +
  labs(x = "Export unit price, FOB (US$/tonne, log scale)", y = NULL) +
  theme_ipsum(axis_title_size = 10, base_size = 9) +
  theme(legend.position    = "none",
        panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(10, 10, 10, 10))

# Inset legend: the _exp legend, height fraction scaled to the 2900 px
# page; vertical centre on the 2007 row (just below the 2000-2006 band),
# box narrowed with the trimmed internal x-range (left 0.63), right edge
# flush
ry_leg_rows_exp <- ry_stats_exp$row[ry_stats_exp$period == 2007]
if (length(ry_leg_rows_exp) == 0) ry_leg_rows_exp <- ceiling(ry_n_rows_exp / 2)
ry_leg_ctr_exp  <- (mean(ry_leg_rows_exp) - 0.55) /
  ((ry_n_rows_exp + 1.9) - 0.55)
p_ridgeline_years_exp <- p_ridge_ry_exp +
  inset_element(p_ridge_legend_exp,
                left   = 0.70,
                bottom = max(0.02, ry_leg_ctr_exp - ry_leg_hgt / 2),
                right  = 1,
                top    = min(0.98, ry_leg_ctr_exp + ry_leg_hgt / 2),
                clip   = FALSE)

# ── Save — five standalone figures ───────────────────────────────────────────
# partners is sized to fit an A4 page inside normal 2.54 cm margins
# (1880 px at 300 dpi, as network_mirrored_desc_stat) with the height of
# its 3 x 2 spaghetti panels; the countries ridgelines are taller for
# their 11 rows; the years ridgelines taller still for their 28 year rows
price_dist_figs    <- list(
  partners                = p_partner,
  ridgeline_countries_imp = p_ridgeline_countries_imp,
  ridgeline_countries_exp = p_ridgeline_countries_exp,
  ridgeline_years_imp     = p_ridgeline_years_imp,
  ridgeline_years_exp     = p_ridgeline_years_exp
)
price_dist_widths  <- c(partners = 1880,
                        ridgeline_countries_imp = 2480,
                        ridgeline_countries_exp = 2480,
                        ridgeline_years_imp = 2480,
                        ridgeline_years_exp = 2480)
price_dist_heights <- c(partners = 2600,
                        ridgeline_countries_imp = 2200,
                        ridgeline_countries_exp = 2200,
                        ridgeline_years_imp = 2900,
                        ridgeline_years_exp = 2900)
for (fig_key in names(price_dist_figs)) {
  for (ext in snakemake@params$ext) {
    out_path <- file.path(
      "results", "network_analysis", "country_lvl", "plot", "01",
      paste0("network_price_", fig_key, ".", ext)
    )
    ggsave(filename = out_path, plot = price_dist_figs[[fig_key]],
           device = ext, create.dir = TRUE,
           width = price_dist_widths[[fig_key]],
           height = price_dist_heights[[fig_key]],
           units = "px", dpi = 300, bg = "white")
  }
}