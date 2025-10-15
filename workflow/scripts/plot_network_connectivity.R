library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")
library("patchwork")

# Load data
data <- read.csv(file = snakemake@input[[1]],
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

    # Produce basic plot
    plot_metric <- data %>%
      filter(cmd == prod) %>%
      # Select relevant columns
      select(all_of(columns)) %>%
      # Convert wide to long dataframe
      melt(id.vars = "period",
           variable.name = "trade_flow",
           value.name = "metric") %>%
      arrange(period) %>%
      # Convert trader_type to factor and order its values
      mutate(trade_flow = factor(trade_flow, levels = order)) %>%

      # Create ggplot line plot
      ggplot(aes(x = period,
                 y = metric,
                 color = trade_flow,
                 label = trade_flow)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = NULL,
                  se = FALSE,
                  linewidth = 1)

    # Collect final y values of geom_smooth to plot legend label
    y_imp <- filter(ggplot_build(plot_metric)$data[[2]],
                    x == max(data$period) & grepl("_imp", label))$y
    y_exp <- filter(ggplot_build(plot_metric)$data[[2]],
                    x == max(data$period) & grepl("_exp", label))$y

    # Collect y_min and y_max for setting axis limits
    y_min <- min(min(ggplot_build(plot_metric)$data[[2]]$y),
                 min(select(filter(data, cmd == prod), contains(metric))))
    y_max <- max(max(ggplot_build(plot_metric)$data[[2]]$y),
                 max(select(filter(data, cmd == prod), contains(metric))))

    # Complete plot with annotation and theme
    plot_metric <- plot_metric +
      # End of chart labels
      annotate("text",
               x = max(data$period) + 0.4,
               y = y_imp,
               label = "Imports",
               hjust = 0,
               size = 5,
               lineheight = .8,
               fontface = "bold",
               color = pal[1]) +

      annotate("text",
               x = max(data$period) + 0.4,
               y = y_exp,
               label = "Exports",
               hjust = 0,
               size = 5,
               lineheight = .8,
               fontface = "bold",
               color = pal[2]) +

      # Set up scale colors, breaks, and limits + themes
      scale_color_manual(values = pal) +
      scale_x_continuous(
        breaks = c(1996, 2000, 2005, 2010, 2015, 2020, max(data$period)),
        labels = c("1996", "2000", "2005", "2010", "2015", "2020", as.character(max(data$period))) # nolint
      ) +
      scale_y_continuous(limits = c(trunc(y_min - 1),
                                    round(y_max + 1)),
                         expand = c(0, 0)) +
      labs(x = "Year",
           y = axis_name) +
      coord_cartesian(clip = "off") +
      theme_ipsum(axis_title_size = 16) +
      theme(legend.position = "none",
            plot.margin = margin(10, 50, 10, 20))

    # Store plot in list
    plot_lst[[match(metric, metrics)]] <- plot_metric

  }

  # Produce composite plot using patchwork
  patchwork_plot <- plot_lst[[1]] + plot_lst[[2]] +
    plot_lst[[3]] + plot_lst[[4]] + plot_layout(ncol = 2)

  # Save plot to all file extensions
  for (ext in snakemake@params$ext) {

    ggsave(
      paste(sub(pattern = "(.*)\\..*$",
                replacement = "\\1",
                basename(snakemake@output[[1]])),
            ext,
            sep = "."),
      plot = patchwork_plot,
      device = ext,
      path = paste(dirname(dirname(snakemake@output[[1]])),
                   fao_division,
                   sep = "/"),
      create.dir = TRUE,
      scale = 1,
      width = 2480 * 2.2,
      height = 1240 * 3,
      units = c("px"),
      dpi = 300,
      limitsize = TRUE,
      bg = "white"
    )
  }
}
