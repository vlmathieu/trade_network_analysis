r_packages <- c("ggplot2",
                "ggh4x",
                "ggrepel",
                "ggpubr",
                "hrbrthemes",
                "scales",
                "patchwork",
                "poweRlaw",
                "maps",
                "geosphere",
                "CoordinateCleaner",
                "dplyr",
                "reshape2",
                "FAOSTAT",
                "utils")

for (pkg in r_packages){
  if (!require(pkg)) {
    install.packages(pkg)
  }
}

library("ggplot2")
library("dplyr")
library("hrbrthemes")
library("reshape2")
library("ggrepel")
library("tidyverse")

library("ggh4x")
library("ggrepel")
library("ggpubr")
library("scales")
library("patchwork")
library("grid")
library("tidyverse")
library("ggstream")
library("showtext")
library("ggtext")

network_contribution <- read.csv(
  file = "/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/output/network_contribution.csv", # nolint
  header = TRUE,
  sep = ";")
contributor_profiles <- read.csv(
  file = "/Users/valentinmathieu/Desktop/wd/trade_network_analysis/results/network_analysis/output/contributor_profiles.csv", # nolint
  header = TRUE,
  sep = ";")
head(network_contribution)

main_contributors <- network_contribution %>%
  filter(
    cmd == 12,
    period %in% (2017:2022),
    contrib_trade_value_imp > quantile(.$contrib_trade_value_imp,
                                       probs = 0.95) |
      contrib_trade_value_exp > quantile(.$contrib_trade_value_imp,
                                         probs = 0.95)
  ) %>%
  select("country") %>%
  {unique(c(.$country))}

extract <- contributor_profiles %>%
  filter(
    cmd == 12,
    period %in% (1996:2022),
    country %in% main_contributors
  )

grouping <- extract %>%
  mutate(edge_ratio = nb_edge_imp / nb_edge_exp) %>%
  group_by(country) %>%
  summarize(mean_edge_ratio = mean(edge_ratio, na.rm = TRUE)) %>%
  mutate(group = ifelse(mean_edge_ratio < 0.6, "producers",
                        ifelse(mean_edge_ratio > 2, "consumers",
                               "intermediates"))) %>%
  group_by(group)
group_name <- group_keys(grouping)$group
grouping <- grouping %>%
  group_by(group) %>%
  group_split() %>%
  setNames(group_name)

country_groups <- list()
for (group in names(grouping)) {
  country_groups[[group]] <- grouping[[group]]$country
}

plot_data <- network_contribution %>%
  filter(cmd == 12,
         country %in% country_groups[[3]]) %>%
  mutate(contrib = (contrib_trade_value_exp +
                      contrib_trade_value_imp) / 2 * 100) %>%
  rowwise() %>%
  mutate(contrib_min = (min(contrib_trade_value_imp,
                            contrib_trade_value_exp) * 100),
         contrib_max = (max(contrib_trade_value_imp,
                            contrib_trade_value_exp) * 100))

span <- .5
plot_data_pred <- plot_data %>%
  group_by(country) %>%
  arrange(country, period) %>%
  nest() %>%
  mutate(
    models = purrr::map(data,
                        loess,
                        formula = contrib ~ period,
                        span = span),
    contrib_pred = purrr::map(models, `[[`, "fitted")
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
  unnest(cols = c(data, contrib_pred, contrib_min_pred, contrib_max_pred))

y_max <- max(plot_data_pred$contrib_max_pred, plot_data_pred$contrib_max)
x_max <- max(plot_data_pred$period)

plot_contribution <-  plot_data_pred %>%

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
  # geom_point(aes(x = period,
  #                y = contrib_pred,
  #                shape = country,
  #                color = country)) +
  geom_point(aes(x = period,
                 y = contrib_min,
                 shape = country,
                 color = country),
                 alpha = 0.3) +
  geom_point(aes(x = period,
                 y = contrib_max,
                 shape = country,
                 color = country),
                 alpha = 0.3) +
  geom_text(data = plot_data_pred %>% filter(period == x_max),
            aes(label = country,
                x = period + 1,
                y = contrib_pred,
                color = "black",
                lineheight = .8,
                fontface = "bold")) +

  # Set up scale colors, breaks, and limits + themes
  scale_color_grey() +
  scale_fill_grey() +
  scale_x_continuous(
    breaks = c(1996, 2000, 2005, 2010, 2015, 2020, 2022),
    labels = c("1996", "2000", "2005", "2010", "2015", "2020", "2022")
  ) +
  scale_y_continuous(limits = c(0, y_max),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Contribution to global traded value (in %)") +
  coord_cartesian(clip = "off") +
  theme_ipsum(axis_title_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(20, 20, 20, 20))

ggsave(
  "contribution.png",
  plot = plot_contribution,
  device = "png",
  path = "/Users/valentinmathieu/Desktop/",
  scale = 1,
  width = 2480 * 1.5,
  height = 1240 * 2,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = "white"
)


#################
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
  main_contributors <- tibble::as.tibble(network_contribution) %>%
    filter(
      cmd == prod,  # nolint
      period >= max_year - 5, # nolint
      contrib_trade_value_imp > quantile(.$contrib_trade_value_imp, # nolint
                                         probs = 1 - min_contrib) |
      contrib_trade_value_exp > quantile(.$contrib_trade_value_imp, # nolint
                                         probs = 1 - min_contrib)
    ) %>%
    select("country") %>%
    {unique(c(.$country))} # nolint

  # Create value to plot: mean contribution, min contribution, max contribution
  # and smooth values for nicer plotting
  shape_data <- network_contribution %>%
    filter(cmd == prod,
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

test <- shape_data(network_contribution = network_contribution, prod = 12)
groups <- grouping_countries(test, contributor_profiles, prod = 12)

plot_lst <- list()
for (group in country_groups) {

  # Store plot in list
  plot_lst[[which(sapply(country_groups,
                         identical,
                         group))]] <- group

}
