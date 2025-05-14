#' @title Enhanced Dot Plot for Single-Cell Data Visualization
#' @description Creates an enhanced dot plot for visualizing gene expression across different cell types or clusters in single-cell data, with support for split visualization.
#' @param seu A Seurat object.
#' @param features A vector of gene names or a list of named vectors for grouped features.
#' @param group.by Column name in seu@meta.data for grouping cells. Default: NULL (uses current Idents).
#' @param split.by Column name in seu@meta.data for splitting the groups. Default: NULL.
#' @param split.by.method Method for visualizing the split groups. Options are "border" or "color". Default: "border".
#' @param nudge_factor Factor to adjust the spacing between split groups. Default: 0.35.
#' @param color_scheme Color scheme for the plot (named vector with "low" and "high"). Default: c(low="white", high="darkblue").
#' @param split.by.colors Colors for split groups. Default: "default".
#' @param center_color Center color for diverging schemes. Default: NULL.
#' @param angle Angle of x-axis labels. Default: NULL (auto).
#' @param hjust Horizontal justification of x-axis labels. Default: NULL (auto).
#' @param vjust Vertical justification of x-axis labels. Default: NULL (auto).
#' @param legend_position Position of the legend. Default: "right".
#' @param plot.margin Margins around the plot. Default: margin(5.5,5.5,5.5,5.5).
#' @param panel.spacing Spacing between facet panels. Default: unit(5, "pt").
#' @param strip.placement Placement of facet labels. Default: "outside".
#' @param border Whether to draw borders around points when split.by is NULL. Default: TRUE.
#' @param border.width Width of point borders. Default: 0.6.
#' @param flip Whether to flip coordinates. Default: FALSE.
#' @param free_space Whether to allow free space in facets. Default: TRUE.
#' @param show_grid Whether to show grid lines. Default: TRUE.
#' @param method Which summary statistic to use: "zscore" or "mean". Default: "zscore".
#' @param min.value Minimum value for the color scale. Default: NULL.
#' @param max.value Maximum value for the color scale. Default: NULL.
#' @param ... Additional arguments passed to theme().
#' @return A ggplot2 object.

DotPlot2 <- function(
  seu,
  features,
  group.by = NULL,
  split.by = NULL,
  split.by.method = "border",
  nudge_factor = 0.35,
  color_scheme = c(low = "white", high = "darkblue"),
  split.by.colors = "default",
  center_color = NULL,
  angle = NULL,
  hjust = NULL,
  vjust = NULL,
  legend_position = "right",
  plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5),
  panel.spacing = unit(5, "pt"),
  strip.placement = "outside",
  border = TRUE,
  border.width = 0.6,
  flip = FALSE,
  free_space = TRUE,
  show_grid = TRUE,
  method = "zscore",
  min.value = NULL,
  max.value = NULL,
  ...
) {
  library(ggplot2); library(reshape2); library(dplyr)

  # Validate features
  features <- validate_features(features, seu)

  # Prepare grouping
  if (is.list(features)) {
    feature_groups <- unlist(lapply(names(features), function(g) rep(g, length(features[[g]]))))
    tp <- unlist(features)
  } else {
    tp <- features
    feature_groups <- NULL
  }
  if (is.null(group.by)) {
    group_levels <- levels(Idents(seu))
  } else {
    group_levels <- levels(factor(seu@meta.data[[group.by]]))
  }
  if (!is.null(split.by)) {
    split_levels <- levels(factor(seu@meta.data[[split.by]]))
    group_values <- if (is.null(group.by)) as.character(Idents(seu)) else seu@meta.data[[group.by]]
    combined <- paste(group_values, seu@meta.data[[split.by]], sep = "___")
    seu[["__internal_combined_group__"]] <- combined
    calc_group.by <- "__internal_combined_group__"
  } else {
    calc_group.by <- group.by
  }

  # Compute percent and statistic
  pct <- feature_percent(seu, tp, group.by = calc_group.by)
  pct.m <- melt(pct, value.name = "pct")
  if (method == "mean" || ncol(pct) == 1) {
    if (ncol(pct) == 1) warning("Only one identity present; using mean expression")
    z <- CalcStats(seu, tp, group.by = calc_group.by, method = "mean") %>%
         as.matrix() %>% melt(value.name = "zscore")
    lab_value <- "Average Expression"
  } else {
    z <- CalcStats(seu, tp, group.by = calc_group.by) %>%
         as.matrix() %>% melt(value.name = "zscore")
    lab_value <- "zscore"
  }

  # Join and split back out
  ToPlot <- inner_join(pct.m, z, by = c("Var1","Var2"))
  if (!is.null(split.by)) {
    ToPlot <- ToPlot %>%
      tidyr::separate(Var2, into = c("group","split"), sep = "___")
    ToPlot$group <- factor(ToPlot$group, levels = group_levels)
    ToPlot$split <- factor(ToPlot$split, levels = split_levels)
    n_splits <- length(split_levels)
    nudge_vals <- seq(-nudge_factor/2, nudge_factor/2, length.out = n_splits)
    ToPlot$nudge <- nudge_vals[ToPlot$split]
  } else {
    ToPlot$group <- ToPlot$Var2; ToPlot$split <- NA; ToPlot$nudge <- 0
  }
  if (!is.null(feature_groups)) {
    ToPlot$FeatureGroup <- factor(feature_groups, levels = unique(feature_groups))
  }
  if (flip) {
    ToPlot$group <- factor(ToPlot$group, levels = rev(unique(ToPlot$group)))
    ToPlot$Var1  <- factor(ToPlot$Var1, levels  = unique(ToPlot$Var1))
  } else {
    ToPlot$group <- factor(ToPlot$group, levels = unique(ToPlot$group))
    ToPlot$Var1  <- factor(ToPlot$Var1, levels  = rev(unique(ToPlot$Var1)))
  }

  # Set color scale limits
  value_range <- range(ToPlot$zscore, na.rm = TRUE)
  if (!is.null(min.value)) value_range[1] <- min.value
  if (!is.null(max.value)) value_range[2] <- max.value

  # Determine labels
  if (is.null(angle)) {
    max_len <- if (flip) max(nchar(levels(ToPlot$Var1))) else max(nchar(levels(ToPlot$group)))
    angle <- if (max_len <= 2) 0 else 45
  }
  if (is.null(hjust)) hjust <- ifelse(angle > 0, 1, ifelse(angle < 0, 0, 0.5))
  if (is.null(vjust)) vjust <- ifelse(abs(angle) == 90, 0.5, 1)

  # Build plot
  if (!is.null(split.by) && split.by.method == "border") {
    p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, fill = zscore, color = split)) +
         geom_point(shape = 21, stroke = border.width,
                    position = position_nudge(x = ToPlot$nudge)) +
         scale_fill_gradient(low = color_scheme["low"], high = color_scheme["high"],
                             limits = value_range, oob = scales::squish) +
         scale_color_disc_auto(split.by.colors, n_splits) +
         labs(color = split.by)
  } else if (!is.null(split.by)) {
    p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, color = split, alpha = zscore)) +
         geom_point(position = position_nudge(x = ToPlot$nudge)) +
         scale_color_disc_auto(split.by.colors, length(split_levels)) +
         scale_alpha(range = c(0.1, 1))
  } else if (border) {
    p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, fill = zscore)) +
         geom_point(shape = 21, color = "black", stroke = border.width) +
         scale_fill_gradient(low = color_scheme["low"], high = color_scheme["high"],
                             limits = value_range, oob = scales::squish)
  } else {
    p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, color = zscore)) +
         geom_point() +
         scale_color_gradient(low = color_scheme["low"], high = color_scheme["high"],
                              limits = value_range, oob = scales::squish)
  }

  p <- p +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
      strip.background = element_rect(fill = NA, size = 0),
      panel.spacing = panel.spacing,
      strip.placement = strip.placement,
      legend.position = legend_position,
      ...
    ) +
    labs(size = "Percent\nexpressed", color = lab_value, fill = lab_value)

  if (!show_grid) p <- p + theme(panel.grid = element_blank())
  if (flip)      p <- p + coord_flip()

  if (!is.null(feature_groups)) {
    facet_scales <- ifelse(flip, "free_x", "free_y")
    facet_space  <- ifelse(free_space, "free", "fixed")
    if (flip) {
      p <- p + facet_grid(cols = vars(FeatureGroup),
                          scales = facet_scales, space = facet_space)
    } else {
      p <- p + facet_grid(rows = vars(FeatureGroup),
                          scales = facet_scales, space = facet_space)
    }
  }

  return(p)
}

# validate_features remains unchanged from your first version