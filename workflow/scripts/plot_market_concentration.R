library("ggplot2")
library("hrbrthemes")
suppressPackageStartupMessages(library("dplyr"))
library("reshape2")
library("patchwork")

# Define color palette
pal <- c("#3D85F7", "#C32E5A")

# Define order for corresponding colors
order <- c("hhi_imp", "hhi_exp")

# Map weight values to human-readable titles
wgt_titles <- c(
  "primary_value" = "Market concentration based on primary value",
  "net_wgt"       = "Market concentration based on net weight"
)

# Function to build one panel for a given weight
build_panel <- function(data, prod, wgt, y_min, y_max, panel_letter, fao_division) {

  # Filter data for the given commodity and weight
  plot_data <- data |>
    filter(cmd == prod) |> # nolint
    filter(weight == wgt) |> # nolint
    # Select relevant columns
    select(c("period", "hhi_imp", "hhi_exp")) |>
    # Convert wide to long dataframe
    melt(id.vars = "period",
         variable.name = "trader_type",
         value.name = "hhi") |>
    arrange(period) |> # nolint
    # Convert trader_type to factor and order its values
    mutate(trader_type = factor(trader_type, levels = order)) # nolint

  # Temporary plot used only to extract LOESS smooth end-values for labels
  tmp_panel <- ggplot(plot_data,
                      aes(x = period, # nolint
                          y = hhi, # nolint
                          color = trader_type, # nolint
                          label = trader_type)) +
    geom_point(alpha = 0.6) +
    suppressMessages(
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                  linewidth = 0.75, span = 0.4)
    )

  # Extract LOESS smooth end-values for label placement
  smooth_data <- ggplot_build(tmp_panel)$data[[2]]
  smooth_imp  <- smooth_data[smooth_data$group == 1, ]
  smooth_exp  <- smooth_data[smooth_data$group == 2, ]
  y_imp       <- smooth_imp$y[nrow(smooth_imp)]
  y_exp       <- smooth_exp$y[nrow(smooth_exp)]

  # Minimum vertical gap between any two labels stacked in the right margin
  y_range <- y_max - y_min
  lab_gap <- 0.06 * y_range

  # 1. Separate the two end-of-curve labels from each other
  raw_gap <- abs(y_exp - y_imp)
  offset  <- if (raw_gap < lab_gap) (lab_gap - raw_gap) / 2 else 0
  y_imp_lbl <- if (y_imp >= y_exp) y_imp + offset else y_imp - offset
  y_exp_lbl <- if (y_exp >= y_imp) y_exp + offset else y_exp - offset

  # 2. Push each end-of-curve label clear of the fixed HHI threshold labels
  #    (0.25, 0.15, 0.01), which also sit in the right margin at x = max + 0.4
  thr_y <- c(0.25, 0.15, 0.01)
  push_from_thr <- function(y) {
    for (t in sort(thr_y)) {
      if (abs(y - t) < lab_gap) y <- if (y >= t) t + lab_gap else t - lab_gap
    }
    y
  }
  y_imp_lbl <- push_from_thr(y_imp_lbl)
  y_exp_lbl <- push_from_thr(y_exp_lbl)

  # 3. Re-separate the two end-of-curve labels if the threshold push brought
  #    them back together
  if (abs(y_imp_lbl - y_exp_lbl) < lab_gap) {
    mid <- (y_imp_lbl + y_exp_lbl) / 2
    if (y_imp_lbl >= y_exp_lbl) {
      y_imp_lbl <- mid + lab_gap / 2
      y_exp_lbl <- mid - lab_gap / 2
    } else {
      y_imp_lbl <- mid - lab_gap / 2
      y_exp_lbl <- mid + lab_gap / 2
    }
  }

  # Break vectors for the redrawn grid: shared y-range (round breaks within the
  # limits) and the fixed x benchmark years
  y_breaks <- scales::extended_breaks()(c(y_min, y_max))
  y_breaks <- y_breaks[y_breaks >= y_min & y_breaks <= y_max]
  x_breaks <- c(1996, 2000, 2005, 2010, 2015, 2020, max(data$period))

  # Build the actual plot — HS band first so it renders behind all other layers
  plot_panel <- ggplot(plot_data,
                       aes(x = period, # nolint
                           y = hhi, # nolint
                           color = trader_type, # nolint
                           label = trader_type)) +

    # HS-revision discontinuity band (net weight panel only)
    # Division 07: disruption 1996–1999; divisions 01/05: disruption 2000–2006
    (if (wgt == "net_wgt" && fao_division == "07") geom_rect(
      aes(xmin = 1996, xmax = 1999, ymin = y_min, ymax = y_max),
      fill = "grey92", alpha = 0.08, inherit.aes = FALSE
    ) else if (wgt == "net_wgt") geom_rect(
      aes(xmin = 2000, xmax = 2006, ymin = y_min, ymax = y_max),
      fill = "grey92", alpha = 0.08, inherit.aes = FALSE
    ) else NULL) +
    # 2007 transient spike (division 07 net weight only)
    (if (wgt == "net_wgt" && fao_division == "07") geom_vline(
      xintercept = 2007,
      color      = "grey85",
      linetype   = "solid",
      linewidth  = 0.5
    ) else NULL) +

    # Grid redrawn as explicit layers (theme grid blanked below) so it sits
    # above the HS band but below the points, curves, and labels
    geom_hline(yintercept = y_breaks, color = "grey90", linewidth = 0.3) +
    geom_vline(xintercept = x_breaks, color = "grey90", linewidth = 0.3) +

    # HHI threshold reference lines — drawn before the curves so the curves
    # render on top of them
    geom_segment(aes(x = min(data$period), y = 0.25,
                     xend = max(data$period), yend = 0.25),
                 color = "grey50", linetype = "dotted") +
    geom_segment(aes(x = min(data$period), y = 0.15,
                     xend = max(data$period), yend = 0.15),
                 color = "grey50", linetype = "dotted") +
    geom_segment(aes(x = min(data$period), y = 0.01,
                     xend = max(data$period), yend = 0.01),
                 color = "grey50", linetype = "dotted") +

    geom_point(alpha = 0.6, size = 0.3) +
    suppressMessages(
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                  linewidth = 0.4, span = 0.4)
    ) +

    # End-of-curve labels (right margin)
    annotate("text",
             x = max(data$period) + 0.4,
             y = y_imp_lbl,
             label = "Imports",
             hjust = 0,
             size = 2.6,
             lineheight = .8,
             fontface = "bold",
             color = pal[1]) +

    annotate("text",
             x = max(data$period) + 0.4,
             y = y_exp_lbl,
             label = "Exports",
             hjust = 0,
             size = 2.6,
             lineheight = .8,
             fontface = "bold",
             color = pal[2]) +

    # HHI threshold labels (right margin)
    annotate("text",
             x = max(data$period) + 0.4,
             y = 0.25,
             label = "HHI = 0.25",
             hjust = 0,
             size = 2.2,
             lineheight = .8,
             fontface = "bold",
             color = "grey50") +

    annotate("text",
             x = max(data$period) + 0.4,
             y = 0.15,
             label = "HHI = 0.15",
             hjust = 0,
             size = 2.2,
             lineheight = .8,
             fontface = "bold",
             color = "grey50") +

    annotate("text",
             x = max(data$period) + 0.4,
             y = 0.01,
             label = "HHI = 0.01",
             hjust = 0,
             size = 2.2,
             lineheight = .8,
             fontface = "bold",
             color = "grey50") +

    scale_color_manual(values = pal) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = c("1996", "2000", "2005", "2010", "2015", "2020", as.character(max(data$period))) # nolint
    ) +
    scale_y_continuous(breaks = y_breaks,
                       limits = c(y_min, y_max),
                       expand = expansion(mult = c(0.04, 0))) +
    labs(x = "Year",
         y = "Herfindahl-Hirschman Index",
         title = paste0(panel_letter, " ", wgt_titles[wgt])) +
    coord_cartesian(clip = "off") +
    theme_ipsum(base_size = 9, axis_title_size = 9) +
    theme(legend.position  = "none",
          plot.title       = element_text(size = 9, face = "bold"),
          plot.margin      = margin(3, 13, 3, 2, unit = "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  return(plot_panel) # nolint
}

for (input_file in snakemake@input) {

  # Derive aggregation level and output root from input path
  agg_lvl     <- basename(dirname(dirname(input_file)))
  output_root <- dirname(dirname(dirname(input_file)))

  # Load data
  data <- read.csv(file = input_file,
                   header = TRUE,
                   sep = ";")

  for (fao_division in snakemake@params$fao_divisions) {

    # Convert fao_division to numeric to filter data
    prod <- as.numeric(fao_division)

    # Compute shared y-axis limits across all weights (5% padding on max)
    data_prod <- data |>
      filter(cmd == prod)
    shared_y_max <- max(data_prod$hhi_imp, data_prod$hhi_exp)
    shared_y_max <- shared_y_max + 0.05 * abs(shared_y_max)
    shared_y_min <- 0

    # Build one panel per weight
    plot_lst <- lapply(seq_along(snakemake@params$wgt), function(i) {
      panel_letters <- c("(a)", "(b)")
      build_panel(data, prod, snakemake@params$wgt[[i]],
                  shared_y_min, shared_y_max, panel_letters[i], fao_division)
    })

    # Assemble composite figure (1x2 layout)
    patchwork_plot <- plot_lst[[1]] + plot_lst[[2]] +
      plot_layout(ncol = 2)

    # Save composite figure to all file extensions
    for (ext in snakemake@params$ext) {

      ggsave(
        filename  = file.path(output_root, agg_lvl,
                              "plot", fao_division,
                              paste("market_concentration", ext, sep = ".")),
        plot      = patchwork_plot,
        device    = ext,
        create.dir = TRUE,
        scale     = 1,
        width     = 190,
        height    = 85,
        units     = "mm",
        dpi       = 600,
        limitsize = TRUE,
        bg        = "white"
      )
    }
  }
}