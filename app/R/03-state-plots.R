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
      if (nrow(bundle$pathways) > 0L) as.character(bundle$pathways$pathway_id[[1L]]) else NULL
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
  if (length(matched) == 0L) return(invisible(FALSE))
  gene <- matched[[1L]]
  state$active_gene(gene)
  if (isTRUE(add_to_set)) {
    current <- state$gene_set()
    state$gene_set(c(current, gene)[!duplicated(casefold_key(c(current, gene)))])
  }
  invisible(TRUE)
}

replace_state_gene_set <- function(state, bundle, genes) {
  parsed <- parse_gene_input(genes, bundle_gene_names(bundle))
  state$gene_set(parsed$genes)
  if (length(parsed$genes) > 0L) state$active_gene(parsed$genes[[1L]])
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
  palette <- unlist(bundle$palette$expression %||% bundle$palette$gene, use.names = TRUE)
  if (length(palette) >= 2L) unname(palette) else c("#d9e4ed", "#8f7ca8", "#b52865")
}

discrete_palette <- function(bundle, field, values) {
  configured <- unlist(bundle$palette[[field]] %||% list(), use.names = TRUE)
  values <- unique(as.character(values))
  fallback <- grDevices::hcl.colors(max(3L, length(values)), "Dark 3")
  names(fallback) <- values
  if (length(configured) == 0L) return(fallback[values])
  if (is.null(names(configured)) || any(!nzchar(names(configured)))) {
    configured <- rep(configured, length.out = length(values))
    names(configured) <- values
  }
  missing <- setdiff(values, names(configured))
  if (length(missing) > 0L) configured[missing] <- fallback[missing]
  configured[values]
}

umap_plot_data <- function(bundle, color_by, gene) {
  cells <- bundle$cells
  if (identical(color_by, "expression")) {
    cells$color_value <- as.numeric(expression_matrix(bundle, gene)[, 1L])
    cells$color_label <- gene
  } else if (identical(color_by, "cluster")) {
    cells$color_value <- factor(as.character(cells$cluster), levels = sort(unique(as.character(cells$cluster))))
    cells$color_label <- "Cluster"
  } else {
    cells$color_value <- factor(as.character(cells$condition), levels = unique(as.character(cells$condition)))
    cells$color_label <- "Condition"
  }
  cells
}

make_umap_plotly <- function(bundle, color_by, gene, selected_cell_ids, source) {
  data <- umap_plot_data(bundle, color_by, gene)
  hover <- paste0(
    "Cell: ", data$cell_id,
    "<br>Cluster: ", data$cluster,
    "<br>Condition: ", data$condition
  )
  if (identical(color_by, "expression")) {
    hover <- paste0(
      hover,
      "<br>", gene, ": ", formatC(data$color_value, digits = 4L, format = "fg")
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
      marker = list(size = 5, opacity = 0.76, line = list(width = 0)),
      source = source,
      showlegend = FALSE
    )
    plot <- plotly::colorbar(plot, title = gene)
  } else {
    palette <- discrete_palette(bundle, color_by, data$color_value)
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
      marker = list(size = 5, opacity = 0.76, line = list(width = 0)),
      source = source
    )
  }

  selected <- data[data$cell_id %in% selected_cell_ids, , drop = FALSE]
  if (nrow(selected) > 0L) {
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
        size = 8,
        color = "rgba(255,255,255,0)",
        line = list(color = "#111820", width = 1.4)
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
      title = "UMAP 2", zeroline = FALSE, showgrid = FALSE,
      scaleanchor = "x", scaleratio = 1
    ),
    margin = list(l = 54, r = 18, t = 20, b = 48),
    hovermode = "closest",
    legend = list(orientation = "h", y = -0.15)
  )
  plot <- plotly::config(
    plot,
    displaylogo = FALSE,
    modeBarButtonsToRemove = c(
      "autoScale2d", "toggleSpikelines", "hoverClosestCartesian",
      "hoverCompareCartesian"
    ),
    toImageButtonOptions = list(filename = paste0("espiviz-umap-", safe_filename(gene)))
  )
  plot <- plotly::event_register(plot, "plotly_selected")
  plot <- plotly::event_register(plot, "plotly_click")
  plotly::event_register(plot, "plotly_deselect")
}

make_umap_ggplot <- function(bundle, color_by, gene, selected_cell_ids) {
  data <- umap_plot_data(bundle, color_by, gene)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = umap_1, y = umap_2))
  if (identical(color_by, "expression")) {
    colors <- expression_palette(bundle)
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(color = color_value), size = 1.2, alpha = 0.78) +
      ggplot2::scale_color_gradientn(colors = colors, name = gene)
  } else {
    palette <- discrete_palette(bundle, color_by, data$color_value)
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(color = color_value), size = 1.2, alpha = 0.78) +
      ggplot2::scale_color_manual(values = palette, name = data$color_label[[1L]])
  }
  if (length(selected_cell_ids) > 0L) {
    selected <- data[data$cell_id %in% selected_cell_ids, , drop = FALSE]
    plot <- plot + ggplot2::geom_point(
      data = selected,
      shape = 21,
      size = 2.1,
      stroke = 0.55,
      color = "#111820",
      fill = NA
    )
  }
  plot +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_minimal(base_family = "sans", base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

make_gene_comparison_plot <- function(comparison) {
  if (nrow(comparison) == 0L) return(NULL)
  long <- rbind(
    data.frame(
      gene = comparison$gene,
      group = "Selected cells",
      mean_expression = comparison$selected_mean,
      stringsAsFactors = FALSE
    ),
    data.frame(
      gene = comparison$gene,
      group = "Remaining cells",
      mean_expression = comparison$remaining_mean,
      stringsAsFactors = FALSE
    )
  )
  long$gene <- factor(long$gene, levels = rev(comparison$gene))
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = mean_expression, y = gene, color = group)
  ) +
    ggplot2::geom_point(size = 2.4, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(
      "Selected cells" = "#b52865",
      "Remaining cells" = "#48677c"
    )) +
    ggplot2::labs(x = "Mean normalized expression", y = NULL, color = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#e8ecef"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )
}
