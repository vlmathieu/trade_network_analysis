library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")
library("ggrepel")

# Function that return a bubble plot
bubble_plot <- function(data, prod, year, size, pal) {

  # Columns to keep
  axis_col <- c("country", "nb_edge_imp", "nb_edge_exp")
  to_match <- paste0(size, c("_imp", "_exp"))
  size_col <- colnames(data)[grep(paste(to_match, collapse = "|"),
                                  colnames(data))]
  col_to_keep <- c(axis_col, size_col)

  # Define size label
  if (grepl("value", size, fixed = TRUE)) {
    size_label <- "Traded value (million US$)"
  } else {
    size_label <- "Traded weight (million tons)"
  }

  # Build plot
  plot_profile <- data %>%
    # Filter for given commodity and period
    filter(cmd == prod, period == year) %>% # nolint
    # Select columns to keep
    select(all_of(col_to_keep)) %>%
    # Divide size column values by a million to scale
    mutate(across(all_of(size_col), ~ .x / 1000000)) %>%
    # Arrange by descending values of import size (imp bubbles on top)
    arrange(desc(size_col[grep("_imp", size_col)])) %>%
    # Change size column names by size for automation
    setNames(gsub(size, "size", names(.))) %>% # nolint
    # Mutate country column to factor type
    mutate(country = factor(country)) %>% # nolint

    # Create plot
    ggplot(aes(x = nb_edge_exp, # nolint
               y = nb_edge_imp, # nolint
               label = country)) +

    # Add segment
    geom_segment(aes(x = 0,
                     y = 0,
                     xend = ceiling(max(data$nb_edge_exp) / 10) * 10,
                     yend = ceiling(max(data$nb_edge_imp) / 10) * 10),
                 color = "grey",
                 linetype = "dotted") +

    # Add bubbles for exports and imports
    geom_point(aes(size = size_exp), alpha = 0.6, color = pal[2]) + # nolint
    geom_point(aes(size = size_imp), alpha = 0.6, color = pal[1]) + # nolint

    # Add country name to each point
    geom_text_repel(size = 2,
                    min.segment.length = 1.2,
                    position = position_nudge_repel(x = 3, y = 3),
                    fontface = "bold") +
    # Add year label to plot
    geom_label(label = as.character(year), x = 114, y = 4, size = 5) +

    # Scale size and axis
    scale_size(range = c(1, 20),
               breaks = c(1000, 2000, 5000, 10000),
               limits = c(0, 20000),
               name = size_label) +
    scale_x_continuous(name = "Number of export partners",
                       limits = c(0, ceiling(max(data$nb_edge_exp) / 10) * 10),
                       breaks = seq(0,
                                    ceiling(max(data$nb_edge_exp) / 10) * 10,
                                    by = 20),
                       expand = c(0, 0)) +
    scale_y_continuous(name = "Number of import partners",
                       limits = c(0, 120),
                       breaks = seq(20,
                                    ceiling(max(data$nb_edge_exp) / 10) * 10,
                                    by = 20),
                       expand = c(0, 0)) +
    theme_ipsum(axis_title_size = 11)
  return(plot_profile)
}

# Load data
data <- read.csv(file = snakemake@input[[1]],
                 header = TRUE,
                 sep = ";")

# Define color palette
pal <- c("#3D85F7", "#C32E5A")

# Get size variable
size <- snakemake@params$size

for (fao_division in snakemake@params$fao_divisions) {
  for (year in snakemake@params$years) {

    # Convert fao_division and year to numeric to filter data
    prod <- as.numeric(fao_division)
    year <- as.numeric(year)

    # Produce plot
    plot_profile <- bubble_plot(data, prod, year, size, pal)

    # Save plot to all file extensions
    for (ext in snakemake@params$ext) {

      ggsave(
        paste(sub(pattern = "(.*)\\_.*$",
                  replacement = "\\1",
                  basename(snakemake@output[[1]])),
              paste(year, ext, sep = "."),
              sep = "_"),
        plot = plot_profile,
        device = ext,
        path = paste(dirname(dirname(dirname(snakemake@output[[1]]))),
                     paste(fao_division, "contributor_profiles", sep = "/"),
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
}
