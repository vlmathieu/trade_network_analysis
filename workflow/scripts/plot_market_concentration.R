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
build_panel <- function(data, prod, wgt, y_min, y_max, panel_letter) {

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

  # Conditional offset: nudge labels apart if closer than 5% of y range
  y_range <- y_max - y_min
  min_gap <- 0.05 * y_range
  raw_gap <- abs(y_exp - y_imp)
  offset  <- if (raw_gap < min_gap) (min_gap - raw_gap) / 2 else 0

  y_imp_lbl <- if (y_imp >= y_exp) y_imp + offset else y_imp - offset
  y_exp_lbl <- if (y_exp >= y_imp) y_exp + offset else y_exp - offset

  # Build the actual plot — geom_rect first so it renders behind all other layers
  plot_panel <- ggplot(plot_data,
                       aes(x = period, # nolint
                           y = hhi, # nolint
                           color = trader_type, # nolint
                           label = trader_type)) +

    # HS-revision discontinuity band (net weight panel only)
    (if (wgt == "net_wgt") geom_rect(
      aes(xmin = 2000, xmax = 2006, ymin = y_min, ymax = y_max),
      fill = "grey85", alpha = 0.08, inherit.aes = FALSE
    ) else NULL) +

    geom_point(alpha = 0.6) +
    suppressMessages(
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                  linewidth = 0.75, span = 0.4)
    ) +

    # End of curve labels
    annotate("text",
             x = max(data$period) + 0.4,
             y = y_imp_lbl,
             label = "Imports",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = pal[1]) +

    annotate("text",
             x = max(data$period) + 0.4,
             y = y_exp_lbl,
             label = "Exports",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = pal[2]) +

    # Horizontal lines to indicate HHI thresholds
    geom_segment(aes(x = min(data$period),
                     y = 0.25,
                     xend = max(data$period),
                     yend = 0.25),
                 color = "black",
                 linetype = "dotted") +
    annotate("text",
             x = max(data$period) + 0.4,
             y = 0.25,
             label = "HHI = 0.25",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = "black") +

    geom_segment(aes(x = min(data$period),
                     y = 0.15,
                     xend = max(data$period),
                     yend = 0.15),
                 color = "black",
                 linetype = "dotted") +
    annotate("text",
             x = max(data$period) + 0.4,
             y = 0.15,
             label = "HHI = 0.15",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = "black") +

    geom_segment(aes(x = min(data$period),
                     y = 0.01,
                     xend = max(data$period),
                     yend = 0.01),
                 color = "black",
                 linetype = "dotted") +
    annotate("text",
             x = max(data$period) + 0.4,
             y = 0.01,
             label = "HHI = 0.01",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = "black") +

    scale_color_manual(values = pal) +
    scale_x_continuous(
      breaks = c(1996, 2000, 2005, 2010, 2015, 2020, max(data$period)),
      labels = c("1996", "2000", "2005", "2010", "2015", "2020", as.character(max(data$period))) # nolint
    ) +
    scale_y_continuous(limits = c(y_min, y_max),
                       expand = c(0, 0)) +
    labs(x = "Year",
         y = "Herfindahl-Hirschman Index",
         title = paste0(panel_letter, " ", wgt_titles[wgt])) +
    coord_cartesian(clip = "off") +
    theme_ipsum(axis_title_size = 11) +
    theme(legend.position  = "none",
          plot.title       = element_text(size = 11, face = "bold"),
          plot.margin      = margin(40, 80, 20, 20),
          panel.ontop      = TRUE,
          panel.background = element_rect(fill = "transparent", colour = NA))

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
                  shared_y_min, shared_y_max, panel_letters[i])
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
        width     = 2480 * 2.2,
        height    = 1240 * 1.5,
        units     = "px",
        dpi       = 300,
        limitsize = TRUE,
        bg        = "white"
      )
    }
  }
}