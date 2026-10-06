library("ggplot2")
library("hrbrthemes")
suppressPackageStartupMessages(library("dplyr"))
library("reshape2")
library("legendry")
library("ggrepel")
library("patchwork")

source(file.path(snakemake@scriptdir, "utils.R"))

# ── Label nudges ────────────────────────────────────────────────────
# One (x, y) pair per country, per panel. These tables hold what the main-text
# figure needs; a rule can shift or add entries for its own variant through the
# nudge_main / nudge_a / nudge_b params, merged on top of these (see
# profile_year_windows.smk). A country that appears only in a supplementary
# window belongs in that rule's override, not here, so these tables stay a
# faithful record of the published figure.
# Both axes are sqrt-transformed, so the same nudge value moves a label further
# in partner counts the further right or higher up the label already sits.
NUDGE_MAIN <- list(
  "Canada"         = c( 1.00, -0.20),  # right of its bubble, closer
  "Malaysia"       = c(-0.90, -0.30),  # bottom-left of its bubble
  "Thailand"       = c(-1.00,  0.00),  # left of its bubble
  "Switzerland"    = c( 1.30,  0.30),  # down-right, clear of zoom A box
  "Türkiye"        = c(-0.80, -0.15),  # a bit lower
  "European Union" = c( 0.00, -0.50),
  "USA"            = c( 0.50, -0.50),
  "Norway"         = c(-0.70,  0.30)
)

NUDGE_A <- list(
  "China"         = c(-0.70,  0.50),
  "India"         = c( 0.70,  0.00),  # right, closer to bubble
  "Hong Kong"     = c( 0.80,  0.00),  # right of its bubble (~east)
  "Japan"         = c( 0.60,  0.00),  # right, closer to bubble
  "Rep. of Korea" = c( 0.50,  0.70),  # top of its bubble, closer
  "Viet Nam"      = c(-0.70,  0.30)   # top-left of its bubble
)

NUDGE_B <- list(
  "New Zealand" = c( 0.35,  0.25),
  "Australia"   = c(-0.40,  0.05),
  "Russia"      = c( 0.00, -0.40),
  "Cameroon"    = c(-0.25, -0.20),
  "Papua N.G."  = c(-0.10,  0.30),  # on top, a bit more distance
  "Gabon"       = c( 0.25,  0.15),
  "Congo"       = c( 0.25,  0.15),
  "Uruguay"     = c(-0.20, -0.20),
  "Ukraine"     = c(-0.20, -0.20),
  "Brazil"      = c( 0.25, -0.15)
)

# Attach nudge_x / nudge_y to a label frame. A country absent from the merged
# table gets no nudge, which also covers the "" rows arrow_phantoms() adds.
apply_nudges <- function(d, base, override = NULL) {
  tbl <- if (is.null(override)) base else
    modifyList(base, lapply(override, function(v) as.numeric(unlist(v))))
  idx  <- match(d$country, names(tbl))
  pick <- function(k) vapply(idx, function(i) {
    if (is.na(i)) 0 else as.numeric(tbl[[i]])[k]
  }, numeric(1))
  d$nudge_x <- pick(1)
  d$nudge_y <- pick(2)
  d
}

# Long-format helper: two rows per country sorted so larger bubble renders first
to_long <- function(d) {
  bind_rows(
    d %>% mutate(trade_type = "Total exported value", size_val = size_exp), # nolint
    d %>% mutate(trade_type = "Total imported value", size_val = size_imp) # nolint
  ) %>% # nolint
    arrange(country, desc(size_val)) # nolint
}

arrow_phantoms <- function(d_arr, steps = c(0.25, 0.50, 0.75)) {
  if (nrow(d_arr) == 0) return(data.frame())
  do.call(bind_rows, lapply(steps, function(t) {
    d_arr %>% transmute( # nolint
      nb_edge_exp = (sqrt(x_start) + t * (sqrt(x_end) - sqrt(x_start)))^2, # nolint
      nb_edge_imp = (sqrt(y_start) + t * (sqrt(y_end) - sqrt(y_start)))^2, # nolint
      country = "", lbl_color = "black", lbl_face = "plain"
    )
  }))
}

trajectory_plot <- function(
  data, prod, year_start, year_end, size, threshold, pal
) {

  if (grepl("deflated", size, fixed = TRUE)) {
    size_label <- "Traded value (million 2015 US$)"
  } else if (grepl("value", size, fixed = TRUE)) {
    size_label <- "Traded value (million US$)"
  } else {
    size_label <- "Traded weight (million tons)"
  }

  size_exp_col <- paste0(size, "_exp")
  size_imp_col <- paste0(size, "_imp")

  prod_data <- data %>% # nolint
    filter(cmd == prod) %>% # nolint
    mutate(country = case_match(
      country, # nolint
      "Russian Federation"   ~ "Russia",
      "China, Hong Kong SAR" ~ "Hong Kong",
      "Papua New Guinea"     ~ "Papua N.G.",
      .default = country
    ))

  # Per-year threshold selection
  year_start_data  <- prod_data %>% filter(period == year_start) # nolint
  global_exp_start <- sum(year_start_data[[size_exp_col]], na.rm = TRUE)
  global_imp_start <- sum(year_start_data[[size_imp_col]], na.rm = TRUE)
  qualifying_start <- year_start_data %>% # nolint
    filter(.data[[size_exp_col]] >= threshold * global_exp_start | # nolint
           .data[[size_imp_col]] >= threshold * global_imp_start) %>% # nolint
    pull(country) # nolint

  year_end_data  <- prod_data %>% filter(period == year_end) # nolint
  global_exp_end <- sum(year_end_data[[size_exp_col]], na.rm = TRUE)
  global_imp_end <- sum(year_end_data[[size_imp_col]], na.rm = TRUE)
  qualifying_end <- year_end_data %>% # nolint
    filter(.data[[size_exp_col]] >= threshold * global_exp_end | # nolint
           .data[[size_imp_col]] >= threshold * global_imp_end) %>% # nolint
    pull(country) # nolint

  # All axis limits from qualifying countries, year_start + year_end only
  qualifying_union     <- union(qualifying_start, qualifying_end)
  qualifying_prod_data <- prod_data %>% # nolint
    filter(period %in% c(year_start, year_end), # nolint
           country %in% qualifying_union) # nolint

  xy_max   <- ceiling(max(qualifying_prod_data$nb_edge_exp,
                          qualifying_prod_data$nb_edge_imp,
                          na.rm = TRUE) * 1.15 / 10) * 10
  xy_min   <- floor(min(qualifying_prod_data$nb_edge_exp,
                        qualifying_prod_data$nb_edge_imp,
                        na.rm = TRUE) / 5) * 5
  # Bubble-area maximum. Overridable per rule so a set of figures meant to be
  # compared with each other can share one scale (profile_year_windows.smk
  # pins the three alternative comparison windows to a common maximum).
  size_max <- if (is.null(snakemake@params$size_max)) {
    ceiling(max(qualifying_prod_data[[size_exp_col]],
                qualifying_prod_data[[size_imp_col]],
                na.rm = TRUE) / 1e6)
  } else {
    as.numeric(snakemake@params$size_max)
  }

  # Zoom A: Asia high-importers (China, Viet Nam, Korea, Japan, India + neighbors) # nolint
  # Verified from data: China (18->44, 55->76), India (14->20, 54), Korea (6->5, 34->31) # nolint
  # y-range overridable per rule (robustness_deflated.smk widens it to bring
  # the 2000 positions of Viet Nam and Türkiye inside the inset)
  zoom_a_x <- if (is.null(snakemake@params$zoom_a_x)) c(3, 52) else unlist(snakemake@params$zoom_a_x) # nolint
  zoom_a_y <- if (is.null(snakemake@params$zoom_a_y)) c(27, 90) else unlist(snakemake@params$zoom_a_y) # nolint

  # Zoom B: low-import exporters (Australia, NZ, PNG, Russia, Cameroon, Brazil, Uruguay) # nolint
  # Verified from data: Russia (33->26, 7->6), Brazil (18->29, 8->5), NZ (21->23, 9->7) # nolint
  zoom_b_x <- if (is.null(snakemake@params$zoom_b_x)) c(10, 36) else unlist(snakemake@params$zoom_b_x) # nolint
  zoom_b_y <- if (is.null(snakemake@params$zoom_b_y)) c(0.5, 10) else unlist(snakemake@params$zoom_b_y) # nolint

  # Dynamic size legend breaks as multiples of 1000 (3 breaks). Fractions are
  # widely spaced (1/8, 1/2, 1) so the nested legend circles — whose radius
  # scales with sqrt(value) — don't overlap (0.2/0.6/1.0 put the small and mid
  # circles too close, e.g. 2000/5000/8000 collided).
  # Overridable per rule (robustness_deflated.smk).
  size_breaks <- if (is.null(snakemake@params$size_breaks)) {
    unique(round(c(size_max * 0.125, size_max * 0.5, size_max) / 1000) * 1000)
  } else {
    unlist(snakemake@params$size_breaks)
  }
  # guide_circles() anchors each label at the top of its circle, so closely
  # spaced breaks give touching labels. When breaks are set by the rule, pad
  # the outer labels with a newline to nudge them down / up (the padding line
  # is half-height, see the guide theme below).
  size_labels <- if (is.null(snakemake@params$size_breaks)) {
    scales::label_comma()
  } else {
    function(b) {
      l <- scales::label_comma()(b)
      n <- length(l)
      if (n >= 3) l <- c(paste0("\n", l[1]), l[2:(n - 1)], paste0(l[n], "\n"))
      l
    }
  }

  # Year slices
  d_start <- prod_data %>% # nolint
    filter(period == year_start, country %in% qualifying_start) %>% # nolint
    mutate(
      size_exp = .data[[size_exp_col]] / 1e6,
      size_imp = .data[[size_imp_col]] / 1e6
    )

  d_end <- prod_data %>% # nolint
    filter(period == year_end, country %in% qualifying_end) %>% # nolint
    mutate(
      size_exp = .data[[size_exp_col]] / 1e6,
      size_imp = .data[[size_imp_col]] / 1e6
    )

  d_start <- d_start %>% mutate(year_group = as.character(year_start)) # nolint
  d_end   <- d_end   %>% mutate(year_group = as.character(year_end))   # nolint

  # Long-format data: within each country, larger bubble row comes first
  d_start_long <- to_long(d_start)
  d_end_long   <- to_long(d_end)

  # Filtered label data for zoom panels (only countries within each zoom region)
  label_buf <- 2
  d_end_a   <- d_end %>% # nolint
    filter(nb_edge_exp >= zoom_a_x[1], nb_edge_exp <= zoom_a_x[2], # nolint
           nb_edge_imp >= zoom_a_y[1], nb_edge_imp <= zoom_a_y[2]) # nolint
  d_start_a <- d_start %>% # nolint
    filter(!country %in% qualifying_end, # nolint
           nb_edge_exp >= zoom_a_x[1], nb_edge_exp <= zoom_a_x[2], # nolint
           nb_edge_imp >= zoom_a_y[1], nb_edge_imp <= zoom_a_y[2]) # nolint
  d_end_b   <- d_end %>% # nolint
    filter(nb_edge_exp >= zoom_b_x[1], nb_edge_exp <= zoom_b_x[2], # nolint
           nb_edge_imp >= zoom_b_y[1], nb_edge_imp <= zoom_b_y[2]) # nolint
  d_start_b <- d_start %>% # nolint
    filter(!country %in% qualifying_end, # nolint
           nb_edge_exp >= zoom_b_x[1], nb_edge_exp <= zoom_b_x[2], # nolint
           nb_edge_imp >= zoom_b_y[1], nb_edge_imp <= zoom_b_y[2]) # nolint

  # Alpha scale values / labels for year legend
  alpha_vals <- setNames(c(0.3, 0.6),
                         c(as.character(year_start), as.character(year_end)))
  alpha_labs <- setNames(
    c(paste0("← ", year_start), as.character(year_end)),
    c(as.character(year_start), as.character(year_end))
  )

  # Arrow segments — only for countries present in both years
  # Midpoint computed in sqrt-display space so arrowhead sits visually centred
  d_arrows <- inner_join(
    d_start %>% select(country, x_start = nb_edge_exp, y_start = nb_edge_imp), # nolint
    d_end   %>% select(country, x_end   = nb_edge_exp, y_end   = nb_edge_imp), # nolint
    by = "country"
  ) %>% mutate( # nolint
    x_mid = ((sqrt(x_start) + sqrt(x_end)) / 2)^2, # nolint
    y_mid = ((sqrt(y_start) + sqrt(y_end)) / 2)^2 # nolint
  )

  # Filtered point + arrow data for zoom panels (± label_buf so near-edge bubbles included) # nolint
  # Necessary because ggrepel runs on ALL geom_point data including coord_cartesian-clipped points # nolint
  d_start_long_a <- d_start_long %>% # nolint
    filter(nb_edge_exp >= zoom_a_x[1] - label_buf, nb_edge_exp <= zoom_a_x[2] + label_buf, # nolint
           nb_edge_imp >= zoom_a_y[1] - label_buf, nb_edge_imp <= zoom_a_y[2] + label_buf) # nolint
  d_end_long_a <- d_end_long %>% # nolint
    filter(nb_edge_exp >= zoom_a_x[1] - label_buf, nb_edge_exp <= zoom_a_x[2] + label_buf, # nolint
           nb_edge_imp >= zoom_a_y[1] - label_buf, nb_edge_imp <= zoom_a_y[2] + label_buf) # nolint
  d_start_long_b <- d_start_long %>% # nolint
    filter(nb_edge_exp >= zoom_b_x[1] - label_buf, nb_edge_exp <= zoom_b_x[2] + label_buf, # nolint
           nb_edge_imp >= zoom_b_y[1] - label_buf, nb_edge_imp <= zoom_b_y[2] + label_buf) # nolint
  d_end_long_b <- d_end_long %>% # nolint
    filter(nb_edge_exp >= zoom_b_x[1] - label_buf, nb_edge_exp <= zoom_b_x[2] + label_buf, # nolint
           nb_edge_imp >= zoom_b_y[1] - label_buf, nb_edge_imp <= zoom_b_y[2] + label_buf) # nolint
  d_arrows_a <- d_arrows %>% # nolint
    filter((x_start >= zoom_a_x[1] - label_buf & x_start <= zoom_a_x[2] + label_buf & # nolint
              y_start >= zoom_a_y[1] - label_buf & y_start <= zoom_a_y[2] + label_buf) | # nolint
             (x_end >= zoom_a_x[1] - label_buf & x_end <= zoom_a_x[2] + label_buf & # nolint
                y_end >= zoom_a_y[1] - label_buf & y_end <= zoom_a_y[2] + label_buf)) # nolint
  d_arrows_b <- d_arrows %>% # nolint
    filter((x_start >= zoom_b_x[1] - label_buf & x_start <= zoom_b_x[2] + label_buf & # nolint
              y_start >= zoom_b_y[1] - label_buf & y_start <= zoom_b_y[2] + label_buf) | # nolint
             (x_end >= zoom_b_x[1] - label_buf & x_end <= zoom_b_x[2] + label_buf & # nolint
                y_end >= zoom_b_y[1] - label_buf & y_end <= zoom_b_y[2] + label_buf)) # nolint

  # p_main label data: only countries whose year_end position is outside both zoom boxes # nolint
  in_zoom_a <- function(d) {
    d$nb_edge_exp >= zoom_a_x[1] & d$nb_edge_exp <= zoom_a_x[2] &
      d$nb_edge_imp >= zoom_a_y[1] & d$nb_edge_imp <= zoom_a_y[2]
  }
  in_zoom_b <- function(d) {
    d$nb_edge_exp >= zoom_b_x[1] & d$nb_edge_exp <= zoom_b_x[2] &
      d$nb_edge_imp >= zoom_b_y[1] & d$nb_edge_imp <= zoom_b_y[2]
  }
  d_end_main   <- d_end   %>% filter(!in_zoom_a(.) & !in_zoom_b(.)) # nolint
  d_start_main <- d_start %>% filter(!country %in% qualifying_end,  # nolint
                                     !in_zoom_a(.) & !in_zoom_b(.)) # nolint

  d_labels_main <- bind_rows(
    d_end_main   %>% mutate(lbl_color = "black",  lbl_face = "bold"),  # nolint
    d_start_main %>% mutate(lbl_color = "grey60", lbl_face = "plain"), # nolint
    arrow_phantoms(d_arrows)
  ) %>%
    apply_nudges(NUDGE_MAIN, snakemake@params$nudge_main)
  d_labels_a <- bind_rows(
    d_end_a   %>% mutate(lbl_color = "black",  lbl_face = "bold"),  # nolint
    d_start_a %>% mutate(lbl_color = "grey60", lbl_face = "plain"), # nolint
    arrow_phantoms(d_arrows_a)
  ) %>%
    apply_nudges(NUDGE_A, snakemake@params$nudge_a)
  d_labels_b <- bind_rows(
    d_end_b   %>% mutate(lbl_color = "black",  lbl_face = "bold"),  # nolint
    d_start_b %>% mutate(lbl_color = "grey60", lbl_face = "plain"), # nolint
    arrow_phantoms(d_arrows_b)
  ) %>%
    apply_nudges(NUDGE_B, snakemake@params$nudge_b)

  # ── Custom draw_key for year legend: nested circles at data$alpha ────────
  draw_key_alpha_circles <- function(data, params, size) {
    al <- data$alpha
    grid::grobTree(
      grid::pointsGrob(0.5, 0.5, pch = 19,
                       gp = grid::gpar(col = pal[1], alpha = al, fontsize = 14)), # nolint
      grid::pointsGrob(0.5, 0.5, pch = 19,
                       gp = grid::gpar(col = pal[2], alpha = al, fontsize = 9))
    )
  }

  # ── Main plot: full range, zoom rectangles A and B, no legend ────────────
  p_main <- ggplot() +
    geom_segment(aes(x = xy_min, y = xy_min, xend = xy_max, yend = xy_max), # nolint
                 color = "grey80", linetype = "dotted") +
    annotate("rect",
             xmin = zoom_a_x[1], xmax = zoom_a_x[2],
             ymin = zoom_a_y[1], ymax = zoom_a_y[2],
             fill = NA, color = "black", linetype = "dashed",
             linewidth = 0.35) +
    annotate("label",
             x = zoom_a_x[1], y = zoom_a_y[2],
             label = "A", hjust = 0, vjust = 1,
             size = 3, fontface = "bold", color = "black",
             fill = "white", linewidth = 0.6,
             label.r = unit(0, "lines"),
             label.padding = unit(0.3, "lines")) +
    annotate("rect",
             xmin = zoom_b_x[1], xmax = zoom_b_x[2],
             ymin = zoom_b_y[1], ymax = zoom_b_y[2],
             fill = NA, color = "black", linetype = "dashed",
             linewidth = 0.35) +
    annotate("label",
             x = zoom_b_x[1], y = zoom_b_y[2],
             label = "B", hjust = 0, vjust = 1,
             size = 3, fontface = "bold", color = "black",
             fill = "white", linewidth = 0.6,
             label.r = unit(0, "lines"),
             label.padding = unit(0.3, "lines")) +
    geom_segment(data = d_arrows,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end), # nolint
                 color = "grey40", linewidth = 0.25) +
    geom_segment(data = d_arrows,
                 aes(x = x_start, y = y_start, xend = x_mid, yend = y_mid), # nolint
                 arrow = arrow(length = unit(0.08, "cm"), type = "open"),
                 color = "grey40", linewidth = 0.25) +
    geom_point(data = d_start_long,
               aes(x = nb_edge_exp, y = nb_edge_imp, # nolint
                   size = size_val, color = trade_type, alpha = year_group), # nolint
               show.legend = FALSE) +
    geom_point(data = d_end_long,
               aes(x = nb_edge_exp, y = nb_edge_imp, # nolint
                   size = size_val, color = trade_type, alpha = year_group), # nolint
               show.legend = FALSE) +
    geom_text_repel(data = d_labels_main,
                    aes(x = nb_edge_exp, y = nb_edge_imp, label = country, # nolint
                        color = lbl_color, fontface = lbl_face),            # nolint
                    nudge_x = d_labels_main$nudge_x,
                    nudge_y = d_labels_main$nudge_y,
                    size = 2.2, bg.color = "white", bg.r = 0.2,
                    segment.colour = "grey55", segment.size = 0.25,
                    min.segment.length = 0,
                    # Low repulsion + strong pull toward the nudged target, so
                    # manual nudge_x/nudge_y placement dominates (labels barely
                    # push each other apart or away from bubbles).
                    box.padding = 0.1, point.padding = 0.15,
                    force = 0.2, force_pull = 3, seed = 42,
                    direction = "both", max.overlaps = Inf,
                    show.legend = FALSE) +
    scale_color_manual(
      values = c(
        "Total exported value" = pal[2],
        "Total imported value" = pal[1],
        "black"  = "black",
        "grey60" = "grey60"
      ),
      breaks = c("Total exported value", "Total imported value"),
      name = NULL
    ) +
    scale_discrete_manual(aesthetics = "fontface",
                          values = c("bold" = "bold", "plain" = "plain"),
                          guide = "none") +
    scale_alpha_manual(values = alpha_vals, labels = alpha_labs,
                       name = NULL, guide = "none") +
    scale_size_continuous(range  = c(0.5, 10),
                          breaks = size_breaks,
                          limits = c(0, size_max),
                          labels = scales::label_comma(),
                          name   = size_label,
                          guide  = "none") +
    scale_x_sqrt(name   = "Number of export partners",
                 limits = c(xy_min, xy_max),
                 breaks = c(0, 5, 10, 20, 40, 60, 80, 100, 120),
                 expand = c(0, 0)) +
    scale_y_sqrt(name   = "Number of import partners",
                 limits = c(xy_min, xy_max),
                 breaks = c(0, 5, 10, 20, 40, 60, 80, 100, 120),
                 expand = c(0, 0)) +
    theme_ipsum(base_size = 8, axis_title_size = 8) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none", aspect.ratio = 1,
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.margin = margin(3, 5, 3, 0, "pt"))

  # ── Zoom A: Asia high-importers ───────────────────────────────────────────
  p_zoom_a <- ggplot() +
    geom_segment(aes(x = xy_min, y = xy_min, xend = xy_max, yend = xy_max), # nolint
                 color = "grey80", linetype = "dotted") +
    annotate("label",
             x = zoom_a_x[1], y = zoom_a_y[2],
             label = "A", hjust = 0, vjust = 1,
             size = 3, fontface = "bold", color = "black",
             fill = "white", linewidth = 0.5,
             label.r = unit(0, "lines"),
             label.padding = unit(0.3, "lines")) +
    geom_segment(data = d_arrows_a,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end), # nolint
                 color = "grey40", linewidth = 0.25) +
    geom_segment(data = d_arrows_a,
                 aes(x = x_start, y = y_start, xend = x_mid, yend = y_mid), # nolint
                 arrow = arrow(length = unit(0.07, "cm"), type = "open"),
                 color = "grey40", linewidth = 0.25) +
    geom_point(data = d_start_long_a,
               aes(x = nb_edge_exp, y = nb_edge_imp, # nolint
                   size = size_val, color = trade_type, alpha = year_group), # nolint
               key_glyph = draw_key_alpha_circles,
               show.legend = c(colour = FALSE, size = FALSE, alpha = TRUE)) +
    geom_point(data = d_end_long_a,
               aes(x = nb_edge_exp, y = nb_edge_imp, # nolint
                   size = size_val, color = trade_type, alpha = year_group), # nolint
               show.legend = c(colour = TRUE, size = TRUE, alpha = FALSE)) +
    geom_text_repel(data = d_labels_a,
                    aes(x = nb_edge_exp, y = nb_edge_imp, label = country, # nolint
                        color = lbl_color, fontface = lbl_face),            # nolint
                    nudge_x = d_labels_a$nudge_x,
                    nudge_y = d_labels_a$nudge_y,
                    size = 2.2, bg.color = "white", bg.r = 0.2,
                    segment.colour = "grey55", segment.size = 0.25,
                    min.segment.length = 0,
                    # Low repulsion + strong pull toward the nudged target, so
                    # manual nudge_x/nudge_y placement dominates (labels barely
                    # push each other apart or away from bubbles).
                    box.padding = 0.1, point.padding = 0.15,
                    force = 0.2, force_pull = 3, seed = 42,
                    direction = "both", max.overlaps = Inf,
                    show.legend = FALSE) +
    scale_color_manual(
      values = c(
        "Total exported value" = pal[2],
        "Total imported value" = pal[1],
        "black"  = "black",
        "grey60" = "grey60"
      ),
      breaks = c("Total exported value", "Total imported value"),
      name = NULL
    ) +
    scale_discrete_manual(aesthetics = "fontface",
                          values = c("bold" = "bold", "plain" = "plain"),
                          guide = "none") +
    scale_alpha_manual(
      values = alpha_vals,
      labels = setNames(
        c(paste0(year_start, "  →"), as.character(year_end)),
        c(as.character(year_start), as.character(year_end))
      ),
      name  = "Year",
      guide = guide_legend(nrow = 1, order = -1)
    ) +
    scale_size_continuous(range  = c(0.5, 10),
                          breaks = size_breaks,
                          limits = c(0, size_max),
                          labels = size_labels,
                          name   = size_label) +
    scale_x_sqrt(name   = "Number of export partners",
                 breaks = c(0, 5, 10, 20, 40, 60, 80, 100, 120),
                 expand = c(0, 0)) +
    scale_y_sqrt(name   = "Number of import partners",
                 breaks = c(0, 5, 10, 20, 40, 60, 80, 100, 120),
                 expand = c(0, 0)) +
    theme_ipsum(base_size = 8, axis_title_size = 8) +
    coord_cartesian(xlim = zoom_a_x, ylim = zoom_a_y) +
    guides(
      colour = guide_legend(override.aes = list(size = 5, alpha = 0.7), order = -2), # nolint
      size   = guide_circles(
        text_position = "right",
        # Half-height padding line: the newline padding in size_labels then
        # nudges the outer labels by a quarter line instead of half a line.
        # No effect on single-line labels.
        theme         = theme(legend.text = element_text(lineheight = 0.45)),
        override.aes  = aes(colour = "grey50", alpha = 0.5)
      )
    ) +
    theme(panel.border  = element_rect(fill = NA, linewidth = 0.5,
                                       color = "grey50"),
          aspect.ratio  = 0.60,
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.margin      = margin(3, 0, 3, 5, "pt"),
          legend.margin    = margin(5, 0, 5, 0, "pt"),
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 8),
          legend.key.size  = unit(0.38, "cm"),
          legend.spacing.y = unit(0.15, "cm"))

  # ── Zoom B: low-import exporters ─────────────────────────────────────────
  p_zoom_b <- ggplot() +
    geom_segment(aes(x = xy_min, y = xy_min, xend = xy_max, yend = xy_max), # nolint
                 color = "grey80", linetype = "dotted") +
    annotate("label",
             x = zoom_b_x[1], y = zoom_b_y[2],
             label = "B", hjust = 0, vjust = 1,
             size = 3, fontface = "bold", color = "black",
             fill = "white", linewidth = 0.5,
             label.r = unit(0, "lines"),
             label.padding = unit(0.3, "lines")) +
    geom_segment(data = d_arrows_b,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end), # nolint
                 color = "grey40", linewidth = 0.25) +
    geom_segment(data = d_arrows_b,
                 aes(x = x_start, y = y_start, xend = x_mid, yend = y_mid), # nolint
                 arrow = arrow(length = unit(0.07, "cm"), type = "open"),
                 color = "grey40", linewidth = 0.25) +
    geom_point(data = d_start_long_b,
               aes(x = nb_edge_exp, y = nb_edge_imp, # nolint
                   size = size_val, color = trade_type, alpha = year_group), # nolint
               show.legend = FALSE) +
    geom_point(data = d_end_long_b,
               aes(x = nb_edge_exp, y = nb_edge_imp, # nolint
                   size = size_val, color = trade_type, alpha = year_group), # nolint
               show.legend = FALSE) +
    geom_text_repel(data = d_labels_b,
                    aes(x = nb_edge_exp, y = nb_edge_imp, label = country, # nolint
                        color = lbl_color, fontface = lbl_face),            # nolint
                    nudge_x = d_labels_b$nudge_x,
                    nudge_y = d_labels_b$nudge_y,
                    size = 2.2, bg.color = "white", bg.r = 0.2,
                    segment.colour = "grey55", segment.size = 0.25,
                    min.segment.length = 0,
                    # Low repulsion + strong pull toward the nudged target, so
                    # manual nudge_x/nudge_y placement dominates (labels barely
                    # push each other apart or away from bubbles).
                    box.padding = 0.1, point.padding = 0.15,
                    force = 0.2, force_pull = 3, seed = 42,
                    direction = "both", max.overlaps = Inf,
                    show.legend = FALSE) +
    scale_color_manual(
      values = c(
        "Total exported value" = pal[2],
        "Total imported value" = pal[1],
        "black"  = "black",
        "grey60" = "grey60"
      ),
      breaks = c("Total exported value", "Total imported value"),
      name = NULL
    ) +
    scale_discrete_manual(aesthetics = "fontface",
                          values = c("bold" = "bold", "plain" = "plain"),
                          guide = "none") +
    scale_alpha_manual(values = alpha_vals, labels = alpha_labs,
                       name = NULL, guide = "none") +
    scale_size_continuous(range  = c(0.5, 10),
                          limits = c(0, size_max),
                          guide  = "none") +
    scale_x_sqrt(name   = "Number of export partners",
                 breaks = c(10, 15, 20, 25, 30, 35),
                 expand = c(0, 0)) +
    scale_y_sqrt(name   = "Number of import partners",
                 breaks = c(0, 5, 10, 20, 40, 60, 80, 100, 120),
                 expand = c(0, 0)) +
    theme_ipsum(base_size = 8, axis_title_size = 8) +
    coord_cartesian(xlim = zoom_b_x, ylim = zoom_b_y) +
    theme(legend.position = "none",
          panel.border    = element_rect(fill = NA, linewidth = 0.5,
                                         color = "grey50"),
          aspect.ratio    = 0.60,
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.margin     = margin(3, 0, 3, 5, "pt"))

  (p_main | (p_zoom_a / p_zoom_b)) +
    plot_layout(guides = "collect", widths = c(1.3, 1))
}

pal        <- unname(PAL_DIRECTION)   # (imp, exp)
size       <- snakemake@params$size
threshold  <- snakemake@params$threshold
year_start <- snakemake@params$year_start
year_end   <- snakemake@params$year_end
out_stem   <- if (is.null(snakemake@params$out_stem)) "contributor_profiles" else snakemake@params$out_stem # nolint

for (input_file in snakemake@input) {

  agg_lvl     <- basename(dirname(dirname(input_file)))
  output_root <- dirname(dirname(dirname(input_file)))

  data <- read.csv(file = input_file, header = TRUE, sep = ";")

  for (fao_division in snakemake@params$fao_divisions) {
    prod <- as.numeric(fao_division)

    plot_profile <- trajectory_plot(
      data, prod, year_start, year_end, size, threshold, pal
    )

    for (ext in snakemake@params$ext) {
      ggsave(
        filename   = paste0(out_stem, ".", ext),
        plot       = plot_profile,
        device     = ext,
        path       = file.path(output_root, agg_lvl, "plot", fao_division),
        create.dir = TRUE,
        scale      = 1,
        width      = 190,
        height     = 100,
        units      = "mm",
        dpi        = 600,
        limitsize  = TRUE,
        bg         = "white"
      )
    }
  }
}
