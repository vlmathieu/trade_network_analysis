library("ggplot2")
library("hrbrthemes")
suppressPackageStartupMessages(library("dplyr"))
library("reshape2")
library("patchwork")

for (input_file in snakemake@input) {

  # Derive aggregation level from input path
  agg_lvl <- basename(dirname(dirname(input_file)))

  # Load data
  data <- read.csv(file = input_file,
                   header = TRUE,
                   sep = ";")

  # Define color palette
  pal <- c("#3D85F7", "#C32E5A")

  for (fao_division in snakemake@params$fao_divisions) {

    # Convert fao_division to numeric to filter data
    prod <- as.numeric(fao_division)

    # Collect list of variable names to plot
    lst <- gsub("_exp", "", names(data)[grep("exp", names(data))])
    metrics <- lst[lst != "nb_edge"]

    # Create emply list to store plot
    plot_lst <- list()

    # Create plot for every connectivity metric
    for (metric in metrics) {

      # Select columns for plot
      columns <- c("period", paste0(metric, c("_imp", "_exp")))

      # Define order for corresponding colors
      order <- paste0(metric, c("_imp", "_exp"))

      # Define name for y axis
      if ("mean" %in% metric) {
        axis_name <- "Mean nb. of trading partners"
      } else if ("var" %in% metric) {
        axis_name <- "Variance of the nb. of trading partners"
      } else if ("skew" %in% metric) {
        axis_name <- "Skewness of the nb. of trading partners"
      } else if ("kurt" %in% metric) {
        axis_name <- "Kurtosis of the nb. of trading partners"
      }

      # Long-format data for this metric (imports + exports series)
      base_df <- data |>
        filter(cmd == prod) |>
        # Select relevant columns
        select(all_of(columns)) |>
        # Convert wide to long dataframe
        melt(id.vars = "period",
             variable.name = "trade_flow",
             value.name = "metric") |>
        arrange(period) |>
        # Convert trader_type to factor and order its values
        mutate(trade_flow = factor(trade_flow, levels = order))

      # Temporary plot used only to extract the LOESS end-values (for the
      # end-of-curve labels) and the y-range (for axis limits and gridlines)
      tmp <- ggplot(base_df,
                    aes(x = period, y = metric,
                        color = trade_flow, label = trade_flow)) +
        geom_point(alpha = 0.6, size = 0.3) +
        suppressMessages(
          geom_smooth(method = "loess",
                      formula = y ~ x,
                      se = FALSE,
                      linewidth = 0.6)
        )

      # Collect final y values of geom_smooth to plot legend label
      y_imp <- filter(ggplot_build(tmp)$data[[2]],
                      x == max(data$period) & grepl("_imp", label))$y
      y_exp <- filter(ggplot_build(tmp)$data[[2]],
                      x == max(data$period) & grepl("_exp", label))$y

      # Collect y_min and y_max for setting axis limits
      y_min <- min(min(ggplot_build(tmp)$data[[2]]$y),
                   min(select(filter(data, cmd == prod), contains(metric))))
      y_max <- max(max(ggplot_build(tmp)$data[[2]]$y),
                   max(select(filter(data, cmd == prod), contains(metric))))

      # Conditionally offset end-of-curve labels if they are too close to each
      # other. The minimum gap is set to 9% of the y range (each label is itself
      # ~5-6% of the range tall at this font size); when the raw gap is smaller,
      # one label is nudged up and the other down by half the shortfall so they
      # remain symmetric around their natural positions.
      y_range   <- y_max - y_min
      min_gap   <- 0.09 * y_range
      raw_gap   <- abs(y_exp - y_imp)
      offset    <- if (raw_gap < min_gap) (min_gap - raw_gap) / 2 else 0

      y_imp_lbl <- if (y_imp >= y_exp) y_imp + offset else y_imp - offset
      y_exp_lbl <- if (y_exp >= y_imp) y_exp + offset else y_exp - offset

      # Metric-specific rounded axis boundaries: floor/ceil the data range to a
      # "nice" step so every point (and the smooth) falls inside a clean framed
      # range whose limits and gridlines are round numbers. The step is derived
      # from the data magnitude, so it adapts per metric (mean, variance,
      # skewness, kurtosis) without hard-coding.
      raw_breaks <- scales::extended_breaks(n = 5)(c(y_min, y_max))
      step       <- raw_breaks[2] - raw_breaks[1]
      y_lo       <- floor(y_min / step) * step
      y_hi       <- ceiling(y_max / step) * step
      y_breaks   <- seq(y_lo, y_hi, by = step)
      x_breaks   <- c(1996, 2000, 2005, 2010, 2015, 2020, max(data$period))

      # Build the final plot. The grid is redrawn as explicit geom_hline /
      # geom_vline layers (theme grid blanked below) so it sits under the
      # points, curves, and labels — as in plot_network_contribution.R. The
      # same break vectors feed the axes so ticks and gridlines align.
      plot_metric <- ggplot(base_df,
                            aes(x = period, y = metric,
                                color = trade_flow, label = trade_flow)) +
        geom_hline(yintercept = y_breaks, color = "grey90", linewidth = 0.3) +
        geom_vline(xintercept = x_breaks, color = "grey90", linewidth = 0.3) +
        geom_point(alpha = 0.6, size = 0.3) +
        suppressMessages(
          geom_smooth(method = "loess",
                      formula = y ~ x,
                      se = FALSE,
                      linewidth = 0.6)
        ) +
        # End of chart labels
        annotate("text",
                 x = max(data$period) + 0.4,
                 y = y_imp_lbl,
                 label = "Imports",
                 hjust = 0,
                 size = 2.8,
                 lineheight = .8,
                 fontface = "bold",
                 color = pal[1]) +

        annotate("text",
                 x = max(data$period) + 0.4,
                 y = y_exp_lbl,
                 label = "Exports",
                 hjust = 0,
                 size = 2.8,
                 lineheight = .8,
                 fontface = "bold",
                 color = pal[2]) +

        # Set up scale colors, breaks, and limits + themes
        scale_color_manual(values = pal) +
        scale_x_continuous(
          breaks = x_breaks,
          labels = c("1996", "2000", "2005", "2010", "2015", "2020", as.character(max(data$period))) # nolint
        ) +
        scale_y_continuous(breaks = y_breaks,
                           limits = c(y_lo, y_hi),
                           expand = c(0, 0)) +
        labs(x = "Year",
             y = axis_name) +
        coord_cartesian(clip = "off") +
        theme_ipsum(base_size = 9, axis_title_size = 9) +
        theme(legend.position  = "none",
              plot.margin      = margin(3, 10, 3, 2, unit = "mm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

      # Store plot in list
      plot_lst[[match(metric, metrics)]] <- plot_metric

    }

    # Produce composite plot using patchwork
    patchwork_plot <- plot_lst[[1]] + plot_lst[[2]] +
      plot_lst[[3]] + plot_lst[[4]] + plot_layout(ncol = 2)

    # Save plot to all file extensions
    for (ext in snakemake@params$ext) {

      ggsave(
        filename = file.path("results", "network_analysis", agg_lvl, "plot",
                             fao_division,
                             paste("network_connectivity", ext, sep = ".")),
        plot = patchwork_plot,
        device = ext,
        create.dir = TRUE,
        scale = 1,
        width = 190,
        height = 150,
        units = "mm",
        dpi = 600,
        limitsize = TRUE,
        bg = "white"
      )
    }
  }
}