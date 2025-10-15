library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")
library("patchwork")
library("tidyverse")

# Function to shape data for plotting
shape_data <- function(network_contribution,
                       prod,
                       time_span = 5,
                       min_contrib = .05,
                       span = .5) {

  # Get last year of the data
  max_year <- max(network_contribution$period)

  # Get the main contributors to trade in the last years
  # (timespan years before last year)
  main_contributors <- tibble::as_tibble(network_contribution) %>%
    dplyr::filter(
      cmd == prod,  # nolint
      period >= max_year - time_span, # nolint
      contrib_trade_value_imp > quantile(.$contrib_trade_value_imp, # nolint
                                         probs = 1 - min_contrib) |
      contrib_trade_value_exp > quantile(.$contrib_trade_value_imp, # nolint
                                         probs = 1 - min_contrib)
    ) %>%
    select("country") %>%
    {unique(c(.$country))} # nolint

  # Create value to plot: mean contribution, min contribution, max contribution
  # and smooth values for nicer plotting
  shape_data <- tibble::as_tibble(network_contribution) %>%
    dplyr::filter(cmd == prod,  # nolint
                  country %in% main_contributors) %>%  # nolint
    mutate(contrib = (contrib_trade_value_exp +  # nolint
                        contrib_trade_value_imp) / 2 * 100) %>%  # nolint
    rowwise() %>%
    mutate(contrib_min = (min(contrib_trade_value_imp,
                              contrib_trade_value_exp) * 100),
           contrib_max = (max(contrib_trade_value_imp,
                              contrib_trade_value_exp) * 100))

  shape_data <- shape_data %>%  # nolint
    group_by(country) %>%  # nolint
    arrange(country, period) %>%  # nolint
    nest() %>%  # nolint
    mutate(
      models = purrr::map(data,
                          loess,
                          formula = contrib ~ period,
                          span = span),
      contrib_pred = purrr::map(models, `[[`, "fitted")  # nolint
    ) %>%
    select(-models) %>%
    mutate(
      models = purrr::map(data,
                          loess,
                          formula = contrib_min ~ period,
                          span = span),
      contrib_min_pred = purrr::map(models, `[[`, "fitted")
    ) %>%
    select(-models) %>%
    mutate(
      models = purrr::map(data,
                          loess,
                          formula = contrib_max ~ period,
                          span = span),
      contrib_max_pred = purrr::map(models, `[[`, "fitted")
    ) %>%
    select(-models) %>%
    unnest(cols = c(data, contrib_pred, contrib_min_pred, contrib_max_pred))  # nolint

  return(shape_data)
}

# Function to separate country in groups according to their trade profile
grouping_countries <- function(network_contribution,
                               contributor_profiles,
                               prod,
                               edge_ratio_min = 0.6,
                               edge_ratio_max = 2) {

  # Extract first and last year of dataset
  min_year <- min(network_contribution$period)
  max_year <- max(network_contribution$period)

  # Get trading countries in data
  contributors <- network_contribution %>% distinct(country) # nolint

  # Get contributors profiles
  extract <- contributor_profiles %>%
    filter(
      cmd == prod,  # nolint
      period %in% (min_year:max_year), # nolint
      country %in% contributors$country  # nolint
    )

  # Group countries : producers, consumers, intermediates; based on their
  # connections to import/export
  grouping <- extract %>%
    mutate(edge_ratio = nb_edge_imp / nb_edge_exp) %>%  # nolint
    group_by(country) %>%  # nolint
    summarize(mean_edge_ratio = mean(edge_ratio, na.rm = TRUE)) %>%  # nolint
    mutate(group = ifelse(mean_edge_ratio < edge_ratio_min, "producers",  # nolint
                          ifelse(mean_edge_ratio > edge_ratio_max, "consumers",
                                 "intermediates"))) %>%
    group_by(group)  # nolint
  group_name <- group_keys(grouping)$group
  grouping <- grouping %>%
    group_by(group) %>%  # nolint
    group_split() %>%
    setNames(group_name)

  country_groups <- list()
  for (group in names(grouping)) {
    country_groups[[group]] <- grouping[[group]]$country
  }

  return(country_groups)
}

# Load data
network_contribution <- read.csv(file = snakemake@input[[1]],
                                 header = TRUE,
                                 sep = ";")
contributor_profiles <- read.csv(file = snakemake@input[[2]],
                                 header = TRUE,
                                 sep = ";")

for (fao_division in snakemake@params$fao_divisions) {

  # Format data
  plot_data <- shape_data(network_contribution,
                          prod = as.numeric(fao_division))

  # Group main contributors
  country_groups <- grouping_countries(plot_data,
                                       contributor_profiles,
                                       prod = as.numeric(fao_division))

  # Get max year and contribution to scale plot axis
  y_max <- round(max(plot_data$contrib_max_pred, plot_data$contrib_max))
  x_max <- max(plot_data$period)

  # Create emply list to store plot
  plot_lst <- list()

  # Create plot for every connectivity metric
  for (group in country_groups) {

    plot_contribution <-  plot_data %>%
      # Filter data to keep main contributors from country_group
      dplyr::filter(
        country %in% group
      ) %>%

      ggplot(aes(x = period,
                 y = contrib_pred,
                 ymin = contrib_min_pred,
                 ymax = contrib_max_pred,
                 group = country)
      ) +
      geom_line(aes(x = period,
                    y = contrib_pred,
                    linetype = country,
                    color = country)) +
      geom_ribbon(aes(x = period,
                      ymin = contrib_min_pred,
                      ymax = contrib_max_pred,
                      fill = country),
                  alpha = 0.5) +
      geom_point(aes(x = period,
                     y = contrib_pred,
                     shape = country,
                     color = country)) +
      geom_point(aes(x = period,
                     y = contrib_min,
                     shape = country,
                     color = country,
                     alpha = 0.2)) +
      geom_point(aes(x = period,
                     y = contrib_max,
                     shape = country,
                     color = country,
                     alpha = 0.2)) +

      # Set up scale colors, breaks, and limits + themes
      scale_color_grey() +
      scale_fill_grey() +
      geom_text(data = plot_data %>%
                  filter(period == x_max,
                         country %in% group),
                aes(label = country,
                    x = period + .5,
                    y = contrib_pred,
                    color = "black",
                    lineheight = .8,
                    fontface = "bold"),
                hjust = 0) +
      scale_x_continuous(
        breaks = c(1996, 2000, 2005, 2010, 2015, 2020, x_max),
        labels = c("1996", "2000", "2005", "2010", "2015", "2020", as.character(x_max)) # nolint
      ) +
      scale_y_continuous(limits = c(0, y_max),
                         expand = c(0, 0)) +
      labs(x = "Year",
           y = "Contribution to global traded value (in %)") +
      coord_cartesian(clip = "off") +
      theme_ipsum(axis_title_size = 11) +
      theme(legend.position = "none",
            plot.margin = margin(10, 100, 10, 20))

    # Store plot in list
    plot_lst[[which(sapply(country_groups,
                           identical,
                           group))]] <- plot_contribution

  }

  # Produce composite plot using patchwork
  if (length(plot_lst) == 3) {
    patchwork_plot <- plot_lst[[1]] + plot_lst[[2]] + plot_lst[[3]] +
      plot_layout(ncol = 1)
    scale <- 3.3
  } else if (length(plot_lst) == 2) {
    patchwork_plot <- plot_lst[[1]] + plot_lst[[2]] +
      plot_layout(ncol = 1)
    scale <- 2.2
  } else {
    patchwork_plot <- plot_lst[[1]] + plot_layout(ncol = 1)
    scale <- 1.1
  }

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
      width = 2480 * 1.1,
      height = 1240 * scale,
      units = c("px"),
      dpi = 300,
      limitsize = TRUE,
      bg = "white"
    )
  }
}
