canonical_gene <- function(bundle, gene, fallback = NULL) {
  universe <- bundle_gene_names(bundle)
  match_index <- match(casefold_key(gene), casefold_key(universe))
  match_index <- match_index[!is.na(match_index)]
  if (length(match_index) > 0L) {
    universe[[match_index[[1L]]]]
  } else {
    fallback %||% universe[[1L]]
  }
}

default_active_gene <- function(bundle) {
  featured <- unlist(bundle$featured_gene_sets, use.names = FALSE)
  preferred <- c("Glul", "EGFP", featured, bundle_gene_names(bundle))
  canonical_gene(bundle, preferred, bundle_gene_names(bundle)[[1L]])
}

new_app_state <- function(bundle) {
  list(
    active_gene = shiny::reactiveVal(default_active_gene(bundle)),
    gene_set = shiny::reactiveVal(character()),
    selected_cells = shiny::reactiveVal(character()),
    active_pathway = shiny::reactiveVal(
      if (nrow(bundle$pathways) > 0L) {
        as.character(bundle$pathways$pathway_id[[1L]])
      } else {
        NULL
      }
    )
  )
}

state_analysis_genes <- function(state) {
  shiny::reactive({
    genes <- state$gene_set()
    if (length(genes) == 0L) state$active_gene() else genes
  })
}

set_state_gene <- function(state, bundle, gene, add_to_set = FALSE) {
  matched <- parse_gene_input(gene, bundle_gene_names(bundle))$genes
  if (length(matched) == 0L) {
    return(invisible(FALSE))
  }
  gene <- matched[[1L]]
  state$active_gene(gene)
  if (isTRUE(add_to_set)) {
    current <- state$gene_set()
    state$gene_set(c(current, gene)[
      !duplicated(casefold_key(c(current, gene)))
    ])
  }
  invisible(TRUE)
}

replace_state_gene_set <- function(state, bundle, genes) {
  parsed <- parse_gene_input(genes, bundle_gene_names(bundle))
  state$gene_set(parsed$genes)
  if (length(parsed$genes) > 0L) {
    state$active_gene(parsed$genes[[1L]])
  }
  invisible(parsed)
}

append_state_gene_set <- function(state, bundle, genes) {
  parsed <- parse_gene_input(genes, bundle_gene_names(bundle))
  combined <- c(state$gene_set(), parsed$genes)
  combined <- combined[!duplicated(casefold_key(combined))]
  state$gene_set(combined)
  invisible(parsed)
}

expression_palette <- function(bundle) {
  palette <- unlist(
    bundle$palette$expression %||% bundle$palette$gene,
    use.names = TRUE
  )
  if (length(palette) >= 2L) {
    unname(palette)
  } else {
    c("#d9e4ed", "#8f7ca8", "#b52865")
  }
}

expression_color_limit <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  limit <- if (length(values) > 0L) max(abs(values)) else 0
  if (!is.finite(limit) || limit <= 0) 1 else limit
}

discrete_palette <- function(bundle, field, values) {
  configured <- unlist(bundle$palette[[field]] %||% list(), use.names = TRUE)
  values <- unique(as.character(values))
  fallback <- grDevices::hcl.colors(max(3L, length(values)), "Dark 3")
  names(fallback) <- values
  if (length(configured) == 0L) {
    return(fallback[values])
  }
  if (is.null(names(configured)) || any(!nzchar(names(configured)))) {
    configured <- rep(configured, length.out = length(values))
    names(configured) <- values
  }
  missing <- setdiff(values, names(configured))
  if (length(missing) > 0L) {
    configured[missing] <- fallback[missing]
  }
  configured[values]
}

umap_plot_data <- function(bundle, color_by, gene) {
  cells <- bundle$cells
  if (identical(color_by, "expression")) {
    cells$color_value <- as.numeric(expression_matrix(bundle, gene)[, 1L])
    cells$color_label <- gene
    cells <- cells[order(cells$color_value, na.last = TRUE), , drop = FALSE]
  } else if (identical(color_by, "detection")) {
    matched_gene <- match_bundle_genes(bundle, gene)$genes[[1L]]
    detected <- as.numeric(bundle$counts[, matched_gene]) > 0
    cells$color_value <- factor(
      ifelse(detected, "Detected", "Not detected"),
      levels = c("Not detected", "Detected")
    )
    cells$color_label <- paste(matched_gene, "raw-count detection")
    cells <- cells[order(cells$color_value), , drop = FALSE]
  } else if (identical(color_by, "cluster")) {
    cells$color_value <- factor(
      as.character(cells$cluster),
      levels = sort(unique(as.character(cells$cluster)))
    )
    cells$color_label <- "Cluster"
  } else {
    cells$color_value <- factor(
      as.character(cells$condition),
      levels = unique(as.character(cells$condition))
    )
    cells$color_label <- "Condition"
  }
  cells
}

umap_point_style <- function(cell_count) {
  cell_count <- max(0L, as.integer(cell_count %||% 0L))
  if (cell_count <= 500L) {
    return(list(
      interactive_size = 5.2,
      interactive_opacity = 0.78,
      static_size = 1.25,
      static_opacity = 0.78
    ))
  }
  if (cell_count <= 2000L) {
    return(list(
      interactive_size = 3.8,
      interactive_opacity = 0.66,
      static_size = 1,
      static_opacity = 0.68
    ))
  }
  list(
    interactive_size = 3,
    interactive_opacity = 0.56,
    static_size = 0.75,
    static_opacity = 0.6
  )
}

selection_overlay_style <- function(cell_count, selected_count) {
  point_style <- umap_point_style(cell_count)
  selected_count <- max(0L, as.integer(selected_count %||% 0L))
  if (selected_count <= 25L) {
    return(list(
      interactive_size = point_style$interactive_size + 3.5,
      interactive_line_width = 1.4,
      interactive_line_color = "#111820",
      static_size = point_style$static_size + 1.35,
      static_stroke = 0.55,
      static_color = "#111820"
    ))
  }
  if (selected_count <= 200L) {
    return(list(
      interactive_size = point_style$interactive_size + 2,
      interactive_line_width = 0.85,
      interactive_line_color = "rgba(17,24,32,0.78)",
      static_size = point_style$static_size + 0.8,
      static_stroke = 0.4,
      static_color = grDevices::adjustcolor("#111820", alpha.f = 0.78)
    ))
  }
  list(
    interactive_size = point_style$interactive_size + 0.8,
    interactive_line_width = 0.45,
    interactive_line_color = "rgba(17,24,32,0.48)",
    static_size = point_style$static_size + 0.3,
    static_stroke = 0.2,
    static_color = grDevices::adjustcolor("#111820", alpha.f = 0.48)
  )
}

cluster_label_positions <- function(cells) {
  if (nrow(cells) == 0L) {
    return(data.frame(
      cluster = character(),
      umap_1 = numeric(),
      umap_2 = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  clusters <- sort(unique(as.character(cells$cluster)))
  labels <- lapply(clusters, function(cluster) {
    candidates <- cells[as.character(cells$cluster) == cluster, , drop = FALSE]
    center_1 <- stats::median(candidates$umap_1)
    center_2 <- stats::median(candidates$umap_2)
    distance <- (candidates$umap_1 - center_1)^2 +
      (candidates$umap_2 - center_2)^2
    nearest <- which.min(distance)
    data.frame(
      cluster = cluster,
      umap_1 = candidates$umap_1[[nearest]],
      umap_2 = candidates$umap_2[[nearest]],
      stringsAsFactors = FALSE
    )
  })
  labels <- do.call(rbind, labels)
  rownames(labels) <- NULL
  labels
}

cell_ids_for_cluster <- function(bundle, cluster) {
  cluster <- compact_character(cluster)
  if (length(cluster) == 0L) {
    return(character())
  }
  bundle$cells$cell_id[
    as.character(bundle$cells$cluster) == as.character(cluster[[1L]])
  ]
}

make_umap_plotly <- function(
  bundle,
  color_by,
  gene,
  selected_cell_ids,
  source,
  gene_data = NULL
) {
  genes <- normalize_plot_genes(bundle, gene, default_active_gene(bundle))
  expression_blend_mode <- identical(color_by, "expression") &&
    length(genes) == 2L
  detection_blend_mode <- identical(color_by, "detection") &&
    length(genes) == 2L
  two_gene_color_mode <- expression_blend_mode || detection_blend_mode
  if (two_gene_color_mode) {
    if (is.null(gene_data)) {
      gene_data <- prepare_plot_gene_data(bundle, genes, selected_cell_ids)
    }
    data <- if (expression_blend_mode) {
      prepare_umap_expression_blend_data(gene_data)
    } else {
      prepare_umap_detection_data(gene_data)
    }
  } else {
    data <- umap_plot_data(bundle, color_by, genes[[1L]])
  }
  point_style <- umap_point_style(nrow(data))
  hover <- paste0(
    "Cell: ",
    data$cell_id,
    "<br>Cluster: ",
    data$cluster,
    "<br>Condition: ",
    data$condition
  )
  if (expression_blend_mode) {
    hover <- paste0(
      hover,
      "<br>",
      genes[[1L]],
      " Log normalized expression: ",
      formatC(data$expression_1, digits = 4L, format = "fg"),
      "<br>",
      genes[[2L]],
      " Log normalized expression: ",
      formatC(data$expression_2, digits = 4L, format = "fg")
    )
  } else if (detection_blend_mode) {
    hover <- paste0(
      hover,
      "<br>",
      genes[[1L]],
      " Raw detected: ",
      ifelse(data$detected_1, "yes", "no"),
      "<br>",
      genes[[2L]],
      " Raw detected: ",
      ifelse(data$detected_2, "yes", "no"),
      "<br>Detection class: ",
      data$blend_class
    )
  } else if (identical(color_by, "detection")) {
    hover <- paste0(
      hover,
      "<br>",
      genes[[1L]],
      " Raw detected: ",
      ifelse(as.character(data$color_value) == "Detected", "yes", "no")
    )
  }
  if (two_gene_color_mode) {
    plot <- plotly::plot_ly(
      data = data,
      x = ~umap_1,
      y = ~umap_2,
      key = ~cell_id,
      text = hover,
      hoverinfo = "text",
      type = "scattergl",
      mode = "markers",
      marker = list(
        size = point_style$interactive_size,
        opacity = point_style$interactive_opacity,
        color = data$blend_color,
        showscale = FALSE,
        line = list(width = 0)
      ),
      source = source,
      showlegend = FALSE
    )
  } else if (identical(color_by, "expression")) {
    color_limit <- expression_color_limit(data$color_value)
    hover <- paste0(
      hover,
      "<br>",
      genes[[1L]],
      " Log normalized expression: ",
      formatC(data$color_value, digits = 4L, format = "fg")
    )
    plot <- plotly::plot_ly(
      data = data,
      x = ~umap_1,
      y = ~umap_2,
      key = ~cell_id,
      text = hover,
      hoverinfo = "text",
      type = "scattergl",
      mode = "markers",
      color = ~color_value,
      colors = expression_palette(bundle),
      marker = list(
        size = point_style$interactive_size,
        opacity = point_style$interactive_opacity,
        cmin = -color_limit,
        cmid = 0,
        cmax = color_limit,
        line = list(width = 0)
      ),
      source = source,
      showlegend = FALSE
    )
    plot <- plotly::colorbar(
      plot,
      title = paste(genes[[1L]], "log normalized expression"),
      limits = c(-color_limit, color_limit)
    )
  } else {
    palette <- if (identical(color_by, "detection")) {
      c(
        "Not detected" = blend_palette()[["Neither detected"]],
        "Detected" = blend_palette()[["Gene 1"]]
      )
    } else {
      discrete_palette(bundle, color_by, data$color_value)
    }
    plot <- plotly::plot_ly(
      data = data,
      x = ~umap_1,
      y = ~umap_2,
      key = ~cell_id,
      text = hover,
      hoverinfo = "text",
      type = "scattergl",
      mode = "markers",
      color = ~color_value,
      colors = palette,
      marker = list(
        size = point_style$interactive_size,
        opacity = point_style$interactive_opacity,
        line = list(width = 0)
      ),
      source = source,
      showlegend = !identical(color_by, "cluster")
    )
  }

  if (identical(color_by, "cluster")) {
    labels <- cluster_label_positions(data)
    plot <- plotly::add_trace(
      plot,
      data = labels,
      x = ~umap_1,
      y = ~umap_2,
      text = ~cluster,
      type = "scatter",
      mode = "markers+text",
      inherit = FALSE,
      marker = list(
        size = 24,
        color = "rgba(255,254,251,0.90)",
        line = list(color = "#17232b", width = 1)
      ),
      textfont = list(color = "#17232b", size = 12),
      textposition = "middle center",
      hoverinfo = "skip",
      showlegend = FALSE,
      name = ""
    )
  }

  selected <- data[data$cell_id %in% selected_cell_ids, , drop = FALSE]
  if (nrow(selected) > 0L) {
    selection_style <- selection_overlay_style(nrow(data), nrow(selected))
    plot <- plotly::add_trace(
      plot,
      data = selected,
      x = ~umap_1,
      y = ~umap_2,
      key = ~cell_id,
      type = "scattergl",
      mode = "markers",
      inherit = FALSE,
      marker = list(
        size = selection_style$interactive_size,
        color = "rgba(255,255,255,0)",
        line = list(
          color = selection_style$interactive_line_color,
          width = selection_style$interactive_line_width
        )
      ),
      text = paste0("Selected cell: ", selected$cell_id),
      hoverinfo = "text",
      showlegend = FALSE
    )
  }
  plot <- plotly::layout(
    plot,
    dragmode = "lasso",
    xaxis = list(title = "UMAP 1", zeroline = FALSE, showgrid = FALSE),
    yaxis = list(
      title = "UMAP 2",
      zeroline = FALSE,
      showgrid = FALSE,
      scaleanchor = "x",
      scaleratio = 1
    ),
    margin = list(l = 54, r = 18, t = 20, b = 48),
    hovermode = "closest",
    legend = list(orientation = "h", y = -0.15)
  )
  plot <- plotly::config(
    plot,
    displaylogo = FALSE,
    modeBarButtonsToRemove = c(
      "autoScale2d",
      "toggleSpikelines",
      "hoverClosestCartesian",
      "hoverCompareCartesian"
    ),
    toImageButtonOptions = list(
      filename = paste0(
        "espiviz-umap-",
        safe_filename(paste(genes, collapse = "-"))
      )
    )
  )
  plot <- plotly::event_register(plot, "plotly_selected")
  plot <- plotly::event_register(plot, "plotly_click")
  plotly::event_register(plot, "plotly_deselect")
}

make_umap_ggplot <- function(
  bundle,
  color_by,
  gene,
  selected_cell_ids,
  gene_data = NULL
) {
  genes <- normalize_plot_genes(bundle, gene, default_active_gene(bundle))
  expression_blend_mode <- identical(color_by, "expression") &&
    length(genes) == 2L
  detection_blend_mode <- identical(color_by, "detection") &&
    length(genes) == 2L
  two_gene_color_mode <- expression_blend_mode || detection_blend_mode
  if (two_gene_color_mode) {
    if (is.null(gene_data)) {
      gene_data <- prepare_plot_gene_data(bundle, genes, selected_cell_ids)
    }
    data <- if (expression_blend_mode) {
      prepare_umap_expression_blend_data(gene_data)
    } else {
      prepare_umap_detection_data(gene_data)
    }
  } else {
    data <- umap_plot_data(bundle, color_by, genes[[1L]])
  }
  point_style <- umap_point_style(nrow(data))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = umap_1, y = umap_2))
  caption <- NULL
  if (two_gene_color_mode) {
    plot <- plot +
      ggplot2::geom_point(
        ggplot2::aes(color = blend_color),
        size = point_style$static_size,
        alpha = point_style$static_opacity
      ) +
      ggplot2::scale_color_identity()
    caption <- if (expression_blend_mode) {
      paste0(
        "Expression blend — gray: low both; orange: high ",
        genes[[1L]],
        "; blue: high ",
        genes[[2L]],
        "; purple: high both. Each gene's log normalized expression is ",
        "scaled independently from its minimum to maximum across cells."
      )
    } else {
      paste0(
        "Detection legend — gray: Neither detected; orange: ",
        genes[[1L]],
        " detected; blue: ",
        genes[[2L]],
        " detected; purple: Both detected. Detection uses raw counts."
      )
    }
  } else if (identical(color_by, "expression")) {
    colors <- expression_palette(bundle)
    color_limit <- expression_color_limit(data$color_value)
    middle_color <- colors[[ceiling(length(colors) / 2)]]
    plot <- plot +
      ggplot2::geom_point(
        ggplot2::aes(color = color_value),
        size = point_style$static_size,
        alpha = point_style$static_opacity
      ) +
      ggplot2::scale_color_gradient2(
        low = colors[[1L]],
        mid = middle_color,
        high = colors[[length(colors)]],
        midpoint = 0,
        limits = c(-color_limit, color_limit),
        name = paste(genes[[1L]], "log normalized expression")
      )
  } else {
    palette <- if (identical(color_by, "detection")) {
      c(
        "Not detected" = blend_palette()[["Neither detected"]],
        "Detected" = blend_palette()[["Gene 1"]]
      )
    } else {
      discrete_palette(bundle, color_by, data$color_value)
    }
    plot <- plot +
      ggplot2::geom_point(
        ggplot2::aes(color = color_value),
        size = point_style$static_size,
        alpha = point_style$static_opacity
      ) +
      ggplot2::scale_color_manual(
        values = palette,
        name = data$color_label[[1L]]
      )
  }
  if (identical(color_by, "cluster")) {
    labels <- cluster_label_positions(data)
    plot <- plot +
      ggplot2::geom_label(
        data = labels,
        ggplot2::aes(x = umap_1, y = umap_2, label = cluster),
        inherit.aes = FALSE,
        size = 3.4,
        fontface = "bold",
        color = "#17232b",
        fill = "#fffefb",
        linewidth = 0.25,
        label.padding = grid::unit(0.18, "lines")
      )
  }
  if (length(selected_cell_ids) > 0L) {
    selected <- data[data$cell_id %in% selected_cell_ids, , drop = FALSE]
    selection_style <- selection_overlay_style(nrow(data), nrow(selected))
    plot <- plot +
      ggplot2::geom_point(
        data = selected,
        shape = 21,
        size = selection_style$static_size,
        stroke = selection_style$static_stroke,
        color = selection_style$static_color,
        fill = NA
      )
  }
  plot +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2", caption = caption) +
    ggplot2::theme_minimal(base_family = "sans", base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = if (
        identical(color_by, "cluster") || two_gene_color_mode
      ) {
        "none"
      } else {
        "bottom"
      },
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

comparison_plot_data <- function(comparison) {
  if (nrow(comparison) == 0L) {
    return(data.frame(
      gene = character(),
      group = character(),
      mean_expression = numeric(),
      detected_pct = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  selected_group <- ifelse(
    comparison$remaining_n > 0L,
    "Selected cells",
    "All cells"
  )
  selected <- data.frame(
    gene = comparison$gene,
    group = selected_group,
    mean_expression = comparison$selected_mean,
    detected_pct = comparison$selected_detected_pct,
    stringsAsFactors = FALSE
  )

  has_remaining <- !is.na(comparison$remaining_n) & comparison$remaining_n > 0L
  remaining <- data.frame(
    gene = comparison$gene[has_remaining],
    group = rep("Remaining cells", sum(has_remaining)),
    mean_expression = comparison$remaining_mean[has_remaining],
    detected_pct = comparison$remaining_detected_pct[has_remaining],
    stringsAsFactors = FALSE
  )

  rbind(selected, remaining)
}

group_summary_plot_data <- function(summary, group_column) {
  required <- c("gene", group_column, "mean_expression", "detected_pct")
  missing <- setdiff(required, names(summary))
  if (length(missing) > 0L) {
    stop(
      "Grouped expression summary is missing: ",
      paste(missing, collapse = ", ")
    )
  }
  group_values <- as.character(summary[[group_column]])
  group_levels <- if (identical(group_column, "cluster")) {
    numeric_clusters <- suppressWarnings(as.numeric(unique(group_values)))
    if (all(!is.na(numeric_clusters))) {
      as.character(sort(numeric_clusters))
    } else {
      sort(unique(group_values))
    }
  } else if (is.factor(summary[[group_column]])) {
    intersect(levels(summary[[group_column]]), unique(group_values))
  } else {
    unique(group_values)
  }
  data.frame(
    gene = as.character(summary$gene),
    group = factor(group_values, levels = group_levels),
    mean_expression = as.numeric(summary$mean_expression),
    detected_pct = as.numeric(summary$detected_pct),
    stringsAsFactors = FALSE
  )
}

comparison_plot_height <- function(gene_count) {
  gene_count <- max(0L, as.integer(gene_count %||% 0L))
  as.integer(max(300L, min(680L, 180L + 20L * gene_count)))
}

make_gene_comparison_plot <- function(comparison) {
  long <- comparison_plot_data(comparison)
  if (nrow(long) == 0L) {
    return(NULL)
  }
  long$gene <- factor(long$gene, levels = rev(unique(comparison$gene)))
  long$group <- factor(
    long$group,
    levels = c("All cells", "Selected cells", "Remaining cells")
  )

  make_expression_dot_plot(long, group_label = NULL)
}

make_expression_dot_plot <- function(data, group_label, rotate_x = FALSE) {
  if (nrow(data) == 0L) {
    return(NULL)
  }
  if (!is.factor(data$gene)) {
    data$gene <- factor(data$gene, levels = rev(unique(data$gene)))
  }
  if (!is.factor(data$group)) {
    data$group <- factor(data$group, levels = unique(data$group))
  }
  color_limit <- expression_color_limit(data$mean_expression)

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = group,
      y = gene,
      color = mean_expression,
      size = detected_pct
    )
  ) +
    ggplot2::geom_point(alpha = 0.92, na.rm = TRUE) +
    ggplot2::scale_color_gradient2(
      low = "#2166AC",
      mid = "#BDBDBD",
      high = "#E31A8C",
      midpoint = 0,
      limits = c(-color_limit, color_limit)
    ) +
    ggplot2::scale_size_continuous(
      range = c(2, 9),
      limits = c(0, 100),
      breaks = c(0, 50, 100)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(4, "cm")
      ),
      size = ggplot2::guide_legend(title.position = "top", nrow = 1)
    ) +
    ggplot2::labs(
      x = group_label,
      y = NULL,
      color = "Mean log normalized expression",
      size = "Detected (%)"
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1L, nrow = 1L),
      color = ggplot2::guide_colorbar(
        order = 2L,
        direction = "horizontal",
        title.position = "top",
        barwidth = grid::unit(4, "cm"),
        barheight = grid::unit(0.25, "cm")
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#e8ecef"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "left",
      axis.text.x = ggplot2::element_text(
        face = "bold",
        color = "#17232b",
        angle = if (isTRUE(rotate_x)) 30 else 0,
        hjust = if (isTRUE(rotate_x)) 1 else 0.5
      )
    )
}

make_group_summary_plot <- function(
  summary,
  group_column,
  group_label,
  rotate_x = FALSE
) {
  data <- group_summary_plot_data(summary, group_column)
  make_expression_dot_plot(data, group_label, rotate_x = rotate_x)
}

prepare_marker_overview <- function(
  bundle,
  marker_cluster,
  top_n = 8L
) {
  top_n <- max(1L, as.integer(top_n))
  marker_cluster <- compact_character(marker_cluster)
  empty <- data.frame(
    gene = character(),
    cluster = character(),
    cell_count = integer(),
    mean_expression = numeric(),
    median_expression = numeric(),
    detected_n = integer(),
    detected_pct = numeric(),
    marker_rank = integer(),
    within_gene_mean = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(marker_cluster) == 0L) {
    return(empty)
  }

  markers <- bundle$markers[
    as.character(bundle$markers$cluster) == marker_cluster[[1L]],
    ,
    drop = FALSE
  ]
  markers <- markers[order(markers$rank, markers$gene), , drop = FALSE]
  markers <- markers[!duplicated(casefold_key(markers$gene)), , drop = FALSE]
  if (nrow(markers) == 0L) {
    return(empty)
  }
  markers <- utils::head(markers, top_n)
  genes <- as.character(markers$gene)
  gene_index <- match(genes, bundle_gene_names(bundle))
  expression <- expression_matrix(bundle, genes)
  counts <- bundle$counts[, gene_index, drop = FALSE]
  colnames(counts) <- genes

  result <- group_expression_summary(
    expression,
    counts,
    bundle$cells$cluster,
    "cluster"
  )
  result$marker_rank <- as.integer(markers$rank[
    match(result$gene, markers$gene)
  ])
  result$within_gene_mean <- ave(
    result$mean_expression,
    result$gene,
    FUN = function(values) {
      finite <- is.finite(values)
      scaled <- rep(NA_real_, length(values))
      if (!any(finite)) {
        return(scaled)
      }
      limits <- range(values[finite])
      if (diff(limits) <= .Machine$double.eps^0.5) {
        scaled[finite] <- 0.5
      } else {
        scaled[finite] <- (values[finite] - limits[[1L]]) / diff(limits)
      }
      scaled
    }
  )
  result$cluster <- as.character(result$cluster)
  result <- result[order(
    result$marker_rank,
    match(result$cluster, sort(unique(result$cluster)))
  ), , drop = FALSE]
  rownames(result) <- NULL
  result
}

marker_overview_height <- function(gene_count) {
  gene_count <- max(0L, as.integer(gene_count %||% 0L))
  as.integer(max(360L, min(620L, 210L + 32L * gene_count)))
}

make_marker_overview_plot <- function(data) {
  if (nrow(data) == 0L) {
    return(NULL)
  }
  cluster_values <- unique(as.character(data$cluster))
  numeric_clusters <- suppressWarnings(as.numeric(cluster_values))
  cluster_levels <- if (all(!is.na(numeric_clusters))) {
    as.character(sort(numeric_clusters))
  } else {
    sort(cluster_values)
  }
  gene_order <- unique(data$gene[order(data$marker_rank)])
  data$cluster <- factor(as.character(data$cluster), levels = cluster_levels)
  data$gene <- factor(as.character(data$gene), levels = rev(gene_order))

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = cluster,
      y = gene,
      color = within_gene_mean,
      size = detected_pct
    )
  ) +
    ggplot2::geom_point(alpha = 0.94, na.rm = TRUE) +
    ggplot2::scale_color_gradientn(
      colors = c("#e4ebef", "#8b86aa", "#9d2857"),
      limits = c(0, 1)
    ) +
    ggplot2::scale_size_continuous(
      range = c(2, 9),
      limits = c(0, 100),
      breaks = c(0, 50, 100)
    ) +
    ggplot2::labs(
      x = "Final cluster",
      y = NULL,
      color = "Relative mean\n(within gene)",
      size = "Detected (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#e8ecef"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.box = "vertical",
      axis.text.x = ggplot2::element_text(face = "bold", color = "#17232b")
    )
}

prepare_cluster_composition <- function(bundle) {
  samples <- unique(as.character(bundle$cells$sample))
  clusters <- unique(as.character(bundle$cells$cluster))
  numeric_clusters <- suppressWarnings(as.numeric(clusters))
  if (all(!is.na(numeric_clusters))) {
    clusters <- as.character(sort(numeric_clusters))
  } else {
    clusters <- sort(clusters)
  }
  grid <- expand.grid(
    sample = samples,
    cluster = clusters,
    stringsAsFactors = FALSE
  )
  observed <- table(
    factor(as.character(bundle$cells$sample), levels = samples),
    factor(as.character(bundle$cells$cluster), levels = clusters)
  )
  grid$cell_count <- as.integer(observed[cbind(
    match(grid$sample, rownames(observed)),
    match(grid$cluster, colnames(observed))
  )])
  sample_totals <- table(
    factor(as.character(bundle$cells$sample), levels = samples)
  )
  grid$sample_total <- as.integer(sample_totals[grid$sample])
  grid$cell_pct <- 100 * safe_ratio(grid$cell_count, grid$sample_total)

  sample_lookup <- unique(bundle$cells[c("sample", "condition")])
  grid$condition <- as.character(sample_lookup$condition[
    match(grid$sample, as.character(sample_lookup$sample))
  ])
  condition_levels <- unique(as.character(bundle$cells$condition))
  grid$condition <- factor(grid$condition, levels = condition_levels)
  grid$sample <- factor(grid$sample, levels = samples)
  grid$cluster <- factor(grid$cluster, levels = rev(clusters))
  grid
}

make_cluster_composition_plot <- function(data, bundle) {
  if (nrow(data) == 0L) {
    return(NULL)
  }
  condition_levels <- unique(as.character(bundle$cells$condition))
  data$condition <- factor(as.character(data$condition), levels = condition_levels)
  data$tile_label <- paste0(
    format(data$cell_count, big.mark = ","),
    "\n",
    formatC(data$cell_pct, digits = 1L, format = "f"),
    "%"
  )
  threshold <- 0.62 * max(data$cell_pct, na.rm = TRUE)
  data$text_color <- ifelse(data$cell_pct >= threshold, "white", "#17232b")

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = sample, y = cluster, fill = cell_pct)
  ) +
    ggplot2::geom_tile(color = "#fffefb", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = tile_label, color = text_color),
      lineheight = 0.92,
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradientn(
      colors = c("#eef2f4", "#8c91ad", "#9d2857"),
      name = "% of sample"
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(condition),
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::labs(
      x = "Biological sample",
      y = "Final cluster",
      caption = paste(
        "Counts and percentages describe recovered cells only;\n",
        "no differential-abundance test is performed.",
        sep = ""
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(0.8, "lines"),
      strip.text = ggplot2::element_text(face = "bold", color = "#17232b"),
      axis.text.x = ggplot2::element_text(
        angle = 32,
        hjust = 1,
        vjust = 1,
        color = "#17232b"
      ),
      plot.caption = ggplot2::element_text(
        color = "#526b7b",
        hjust = 0
      ),
      legend.position = "top"
    )
}
