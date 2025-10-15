library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")

# Load data
data <- read.csv(file = snakemake@input[[1]],
                 header = TRUE,
                 sep = ";")

# Define color palette
pal <- c("#3D85F7", "#C32E5A")

# Define order for corresponding colors
order <- c("hhi_imp", "hhi_exp")

for (fao_division in snakemake@params$fao_divisions) {

  # Convert fao_division to numeric to filter data
  prod <- as.numeric(fao_division)

  # Data processing for plotting
  plot_mkt <- data %>%
    # Filter data for the given commodity
    filter(cmd == prod) %>%
    # Select relevant columns
    select(c("period", "hhi_imp", "hhi_exp")) %>%
    # Convert wide to long dataframe
    melt(id.vars = "period",
         variable.name = "trader_type",
         value.name = "hhi") %>%
    arrange(period) %>%
    # Convert trader_type to factor and order its values
    mutate(trader_type = factor(trader_type, levels = order)) %>%

    # Create ggplot line plot
    ggplot(aes(x = period,
               y = hhi,
               color = trader_type,
               label = trader_type)) +
    geom_line() +

    # End of chart labels
    annotate("text",
             x = max(data$period) + 0.4,
             y = data[data$period == max(data$period) &
                        data$cmd == prod, ]$hhi_imp,
             label = "Imports",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = pal[1]) +

    annotate("text",
             x = max(data$period) + 0.4,
             y = data[data$period == max(data$period) &
                        data$cmd == prod, ]$hhi_exp,
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
    scale_y_continuous(limits = c(0,
                                  round(max(data[data$cmd == prod, ]$hhi_imp,
                                            data[data$cmd == prod, ]$hhi_exp) +
                                          0.05, 1)),
                       expand = c(0, 0)) +
    labs(x = "Year",
         y = "Herfindahl-Hirschman Index") +
    coord_cartesian(clip = "off") +
    theme_ipsum(axis_title_size = 11) +
    theme(legend.position = "none",
          plot.margin = margin(40, 80, 20, 20))

  # Save plot to all file extensions
  for (ext in snakemake@params$ext) {

    ggsave(
      paste(sub(pattern = "(.*)\\..*$",
                replacement = "\\1",
                basename(snakemake@output[[1]])),
            ext,
            sep = "."),
      plot = plot_mkt,
      device = ext,
      path = paste(dirname(dirname(snakemake@output[[1]])),
                   fao_division,
                   sep = "/"),
      create.dir = TRUE,
      scale = 1,
      width = 2480 * 1.1,
      height = 1240 * 1.5,
      units = c("px"),
      dpi = 300,
      limitsize = TRUE,
      bg = "white"
    )
  }
}
