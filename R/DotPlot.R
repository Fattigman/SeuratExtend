#' @title Enhanced Dot Plot for Single-Cell Data Visualization (with manual scale limits)
#' @description Creates an enhanced dot plot for visualizing gene expression across different cell types or clusters in single-cell data, with support for split visualization and manual color scale cutoffs.
#' @param seu A Seurat object.
#' @param features A vector of gene names or a list of named vectors for grouped features.
#' @param group.by Column name in seu@meta.data for grouping cells. Default: NULL (uses current Idents).
#' @param split.by Column name in seu@meta.data for splitting the groups. Default: NULL.
#' @param split.by.method Method for visualizing the split groups. Options are "border" or "color". Default: "border".
#' @param nudge_factor Factor to adjust the spacing between split groups. Default: 0.35.
#' @param color_scheme Color scheme for the plot. Accepts viridis options, Brewer palettes, or custom low/high vectors.
#' @param split.by.colors Colors for split groups when using "border" or "color" methods.
#' @param center_color Center color for diverging schemes. Default: NULL.
#' @param angle Angle of x-axis labels. Default: auto.
#' @param hjust Horizontal justification. Default: auto.
#' @param vjust Vertical justification. Default: auto.
#' @param legend_position Legend placement. Default: "right".
#' @param plot.margin Plot margins.
#' @param panel.spacing Spacing between facets.
#' @param strip.placement Placement of facet strips.
#' @param border Draw borders around points (when no split.by). Default: TRUE.
#' @param border.width Border width. Default: 0.6.
#' @param flip Flip coordinates. Default: FALSE.
#' @param free_space Allow free space in facets. Default: TRUE.
#' @param show_grid Show grid lines. Default: TRUE.
#' @param method Statistic for expression ("zscore" or "mean"). Default: "zscore".
#' @param min.value Minimum value for color scale cutoff. Default: NULL (auto by data).
#' @param max.value Maximum value for color scale cutoff. Default: NULL (auto by data).
#' @param ... Additional theme arguments.
#' @return A ggplot2 object.
#' @export
DotPlot2 <- function(
  seu,
  features,
  group.by = NULL,
  split.by = NULL,
  split.by.method = "border",
  nudge_factor = 0.35,
  color_scheme = c(low="white", high="darkblue"),
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
  features <- validate_features(features, seu)
  if (is.list(features)) {
    feature_groups <- unlist(lapply(names(features), function(g) rep(g, length(features[[g]]))))
    tp <- unlist(features)
  } else {
    tp <- features; feature_groups <- NULL
  }

  if (is.null(group.by)) {
    groups <- Idents(seu); group_levels <- levels(groups)
  } else {
    group_levels <- levels(factor(seu@meta.data[[group.by]]))
  }

  if (!is.null(split.by)) {
    split_levels <- levels(factor(seu@meta.data[[split.by]]))
    group_values <- if (is.null(group.by)) as.character(Idents(seu)) else seu@meta.data[[group.by]]
    combined_group <- paste(group_values, seu@meta.data[[split.by]], sep = "___")
    seu[["__internal_combined_group__"]] <- combined_group
    calc_group.by <- "__internal_combined_group__"
  } else {
    calc_group.by <- group.by
  }

  pct  <- feature_percent(seu, tp, group.by = calc_group.by)
  pct.m <- reshape2::melt(pct, value.name = "pct")

  if (method == "mean") {
    z <- CalcStats(seu, tp, group.by = calc_group.by, method = "mean") %>% as.matrix() %>% reshape2::melt(value.name = "zscore")
    lab_value <- "Average Expression"
  } else {
    z <- CalcStats(seu, tp, group.by = calc_group.by) %>% as.matrix() %>% reshape2::melt(value.name = "zscore")
    lab_value <- "zscore"
  }

  ToPlot <- dplyr::inner_join(pct.m, z, by = c("Var1","Var2"))

  if (!is.null(split.by)) {
    ToPlot <- tidyr::separate(ToPlot, Var2, into = c("group","split"), sep = "___")
    ToPlot$group <- factor(ToPlot$group, levels = group_levels)
    ToPlot$split <- factor(ToPlot$split, levels = split_levels)
    n_splits <- length(split_levels)
    nudge_values <- seq(-nudge_factor/2, nudge_factor/2, length.out = n_splits)
    ToPlot$nudge <- nudge_values[ToPlot$split]
  } else {
    ToPlot$group <- ToPlot$Var2; ToPlot$split <- NA; ToPlot$nudge <- 0
  }

  if (!is.null(feature_groups)) {
    ToPlot$FeatureGroup <- rep(feature_groups, times = ncol(pct))
    ToPlot$FeatureGroup <- factor(ToPlot$FeatureGroup, unique(ToPlot$FeatureGroup))
  }

  if (flip) {
    ToPlot$group <- factor(ToPlot$group, rev(unique(ToPlot$group)))
    ToPlot$Var1  <- factor(ToPlot$Var1, unique(ToPlot$Var1))
  } else {
    ToPlot$group <- factor(ToPlot$group, unique(ToPlot$group))
    ToPlot$Var1  <- factor(ToPlot$Var1, rev(unique(ToPlot$Var1)))
  }

  lims <- range(ToPlot$zscore, na.rm = TRUE)
  if (!is.null(min.value)) lims[1] <- min.value
  if (!is.null(max.value)) lims[2] <- max.value

  if (!is.null(split.by) && split.by.method == "border") {
    p <- ggplot2::ggplot(ToPlot, aes(group, Var1, size=pct, fill=zscore, color=split)) +
      ggplot2::geom_point(shape=21, stroke=border.width, position=ggplot2::position_nudge(x=ToPlot$nudge))
    fill_scale <- ggplot2::scale_fill_gradient(low = color_scheme["low"], high = color_scheme["high"], limits = lims, oob = scales::squish)
    p <- p + fill_scale + scale_color_disc_auto(split.by.colors, n_splits) + ggplot2::labs(color=split.by)
  } else if (!is.null(split.by)) {
    p <- ggplot2::ggplot(ToPlot, aes(group, Var1, size=pct, color=split, alpha=zscore)) +
      ggplot2::geom_point(position=ggplot2::position_nudge(x=ToPlot$nudge))
    p <- p + scale_color_disc_auto(split.by.colors, length(split_levels)) + ggplot2::scale_alpha(range=c(0.1,1))
  } else if (border) {
    p <- ggplot2::ggplot(ToPlot, aes(group, Var1, size=pct, fill=zscore)) + ggplot2::geom_point(shape=21, color="black", stroke=border.width)
    fill_scale <- ggplot2::scale_fill_gradient(low = color_scheme["low"], high = color_scheme["high"], limits = lims, oob = scales::squish)
    p <- p + fill_scale
  } else {
    p <- ggplot2::ggplot(ToPlot, aes(group, Var1, size=pct, color=zscore)) + ggplot2::geom_point()
    p <- p + ggplot2::scale_color_gradient(low = color_scheme["low"], high = color_scheme["high"], limits = lims, oob = scales::squish)
  }

  p <- p + ggplot2::theme_bw() + ggplot2::theme(axis.title=ggplot2::element_blank(), axis.text.x=ggplot2::element_text(angle=angle %||% 45, hjust=hjust %||% 1, vjust=vjust %||% 1), strip.background=ggplot2::element_rect(fill=NA, size=0), panel.spacing=panel.spacing, strip.placement=strip.placement, legend.position=legend_position) + ggplot2::labs(size="Percent\nexpressed", color=lab_value, fill=lab_value) + ggplot2::theme(...)

  if (!show_grid) p <- p + ggplot2::theme(panel.grid=ggplot2::element_blank())
  if (flip) p <- p + ggplot2::coord_flip()

  if (!is.null(feature_groups)) {
    facet_args <- list(scales=ifelse(flip,"free_x","free_y"), space=ifelse(free_space,"free","fixed"))
    if (flip) {
      p <- p + ggplot2::facet_grid(cols=ggplot2::vars(FeatureGroup), scales=facet_args$scales, space=facet_args$space)
    } else {
      p <- p + ggplot2::facet_grid(rows=ggplot2::vars(FeatureGroup), scales=facet_args$scales, space=facet_args$space)
    }
  }

  return(p)
}

validate_features <- function(features, seu) {
  if (is.list(features)) {
    # Process each group
    validated_features <- lapply(features, function(f) {
      existing <- intersect(f, rownames(seu))
      if (length(existing) == 0) return(NULL)

      # Warn about missing features
      missing <- setdiff(f, rownames(seu))
      if (length(missing) > 0) {
        warning(sprintf("The following requested variables were not found: %s",
                        paste(missing, collapse = ", ")))
      }
      return(existing)
    })

    # Remove empty groups
    empty_groups <- names(validated_features)[sapply(validated_features, is.null)]
    if (length(empty_groups) > 0) {
      warning(sprintf("The following groups had no valid features and were removed: %s",
                      paste(empty_groups, collapse = ", ")))
      validated_features <- validated_features[!sapply(validated_features, is.null)]
    }

    if (length(validated_features) == 0) {
      stop("No valid features found in any group")
    }

    # Check for duplicates across all groups
    all_features <- unlist(validated_features)
    duplicates <- all_features[duplicated(all_features)]
    if (length(duplicates) > 0) {
      warning(sprintf("Removing duplicate features (keeping first occurrence): %s",
                      paste(unique(duplicates), collapse = ", ")))
      # Keep first occurrence of each feature
      seen <- c()
      validated_features <- lapply(validated_features, function(f) {
        unique_f <- setdiff(f, seen)
        seen <<- c(seen, unique_f)
        return(unique_f)
      })
      # Remove any groups that became empty after duplicate removal
      validated_features <- validated_features[sapply(validated_features, length) > 0]
    }

    return(validated_features)

  } else {
    # Process single vector of features
    existing <- intersect(features, rownames(seu))

    if (length(existing) == 0) {
      stop("No valid features found")
    }

    # Warn about missing features
    missing <- setdiff(features, rownames(seu))
    if (length(missing) > 0) {
      warning(sprintf("The following requested variables were not found: %s",
                      paste(missing, collapse = ", ")))
    }

    # Handle duplicates
    duplicates <- features[duplicated(features)]
    if (length(duplicates) > 0) {
      warning(sprintf("Removing duplicate features (keeping first occurrence): %s",
                      paste(unique(duplicates), collapse = ", ")))
      existing <- unique(existing)
    }

    return(existing)
  }
}