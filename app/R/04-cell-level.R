normalize_plot_genes <- function(bundle, genes, fallback) {
  parsed <- parse_gene_input(
    genes,
    bundle_gene_names(bundle),
    max_genes = 2L
  )
  if (length(parsed$genes) > 0L) {
    return(parsed$genes)
  }
  canonical_gene(bundle, fallback)
}

prepare_plot_gene_data <- function(
  bundle,
  genes,
  selected_cell_ids = character()
) {
  genes <- normalize_plot_genes(bundle, genes, default_active_gene(bundle))
  expression <- expression_matrix(bundle, genes)
  gene_index <- match(genes, bundle_gene_names(bundle))
  counts <- as.matrix(bundle$counts[, gene_index, drop = FALSE])
  selected_cell_ids <- intersect(
    as.character(bundle$cells$cell_id),
    as.character(selected_cell_ids %||% character())
  )

  data <- bundle$cells[, c(
    "cell_id",
    "umap_1",
    "umap_2",
    "cluster",
    "condition",
    "sample"
  ), drop = FALSE]
  data$expression_1 <- as.numeric(expression[, 1L])
  data$detected_1 <- as.numeric(counts[, 1L]) > 0
  if (length(genes) == 2L) {
    data$expression_2 <- as.numeric(expression[, 2L])
    data$detected_2 <- as.numeric(counts[, 2L]) > 0
  } else {
    data$expression_2 <- rep(NA_real_, nrow(data))
    data$detected_2 <- rep(FALSE, nrow(data))
  }
  data$selected <- data$cell_id %in% selected_cell_ids
  attr(data, "genes") <- genes
  data
}

blend_palette <- function() {
  c(
    "Neither detected" = "#D9DEE2",
    "Gene 1" = "#D55E00",
    "Gene 2" = "#0072B2",
    "Both detected" = "#6F4C9B"
  )
}

scale_detected_strength <- function(values, detected, minimum = 0.28) {
  values <- as.numeric(values)
  detected <- as.logical(detected)
  result <- numeric(length(values))
  index <- which(detected & is.finite(values))
  if (length(index) == 0L) {
    return(result)
  }
  limits <- range(values[index])
  if (diff(limits) <= .Machine$double.eps^0.5) {
    result[index] <- 1
    return(result)
  }
  scaled <- (values[index] - limits[[1L]]) / diff(limits)
  result[index] <- minimum + (1 - minimum) * pmin(1, pmax(0, scaled))
  result
}

blend_rgb <- function(strength_1, strength_2, palette = blend_palette()) {
  neutral <- grDevices::col2rgb(palette[["Neither detected"]])[, 1L]
  gene_1 <- grDevices::col2rgb(palette[["Gene 1"]])[, 1L]
  gene_2 <- grDevices::col2rgb(palette[["Gene 2"]])[, 1L]
  both <- grDevices::col2rgb(palette[["Both detected"]])[, 1L]
  colors <- vapply(seq_along(strength_1), function(index) {
    first <- strength_1[[index]]
    second <- strength_2[[index]]
    total <- first + second
    if (total <= 0) {
      return(unname(palette[["Neither detected"]]))
    }
    target <- (gene_1 * first + gene_2 * second) / total
    overlap <- min(first, second)
    target <- target * (1 - overlap) + both * overlap
    opacity <- max(first, second)
    mixed <- neutral * (1 - opacity) + target * opacity
    grDevices::rgb(
      mixed[[1L]],
      mixed[[2L]],
      mixed[[3L]],
      maxColorValue = 255
    )
  }, character(1L))
  toupper(colors)
}

prepare_umap_blend_data <- function(gene_data) {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    stop("Two plot genes are required for a blended UMAP.", call. = FALSE)
  }
  data <- gene_data
  data$strength_1 <- scale_detected_strength(
    data$expression_1,
    data$detected_1
  )
  data$strength_2 <- scale_detected_strength(
    data$expression_2,
    data$detected_2
  )
  data$blend_strength <- data$strength_1 + data$strength_2
  data$blend_class <- ifelse(
    data$detected_1 & data$detected_2,
    "Both detected",
    ifelse(
      data$detected_1,
      paste(genes[[1L]], "detected"),
      ifelse(
        data$detected_2,
        paste(genes[[2L]], "detected"),
        "Neither detected"
      )
    )
  )
  data$blend_color <- blend_rgb(data$strength_1, data$strength_2)
  data <- data[order(data$blend_strength, data$cell_id), , drop = FALSE]
  rownames(data) <- NULL
  attr(data, "genes") <- genes
  data
}

blend_legend_ui <- function(genes) {
  genes <- compact_character(genes)
  if (length(genes) != 2L) {
    return(NULL)
  }
  palette <- blend_palette()
  labels <- c(
    "Neither detected",
    paste(genes[[1L]], "detected"),
    paste(genes[[2L]], "detected"),
    "Both detected"
  )
  colors <- unname(palette[c(
    "Neither detected",
    "Gene 1",
    "Gene 2",
    "Both detected"
  )])
  htmltools::div(
    class = "blend-legend",
    role = "list",
    `aria-label` = paste("Two-gene expression blend for", paste(genes, collapse = " and ")),
    lapply(seq_along(labels), function(index) {
      htmltools::div(
        class = "blend-legend-item",
        role = "listitem",
        htmltools::span(
          class = "blend-swatch",
          style = paste0("background-color:", colors[[index]], ";"),
          `aria-hidden` = "true"
        ),
        htmltools::span(labels[[index]])
      )
    })
  )
}

plot_gene_data_long <- function(gene_data) {
  genes <- attr(gene_data, "genes")
  rows <- lapply(seq_along(genes), function(index) {
    data.frame(
      cell_id = as.character(gene_data$cell_id),
      cluster = gene_data$cluster,
      condition = gene_data$condition,
      sample = as.character(gene_data$sample),
      gene = genes[[index]],
      expression = as.numeric(gene_data[[paste0("expression_", index)]]),
      detected = as.logical(gene_data[[paste0("detected_", index)]]),
      selected = as.logical(gene_data$selected),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  data$gene <- factor(data$gene, levels = genes)
  rownames(data) <- NULL
  data
}

make_violin_plot <- function(gene_data, bundle) {
  data <- plot_gene_data_long(gene_data)
  if (nrow(data) == 0L) {
    return(NULL)
  }
  clusters <- sort(unique(as.character(data$cluster)))
  data$cluster <- factor(as.character(data$cluster), levels = clusters)
  palette <- discrete_palette(bundle, "cluster", clusters)
  selected <- data[data$selected, , drop = FALSE]

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = cluster, y = expression, fill = cluster)
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = FALSE,
      linewidth = 0.3,
      color = "#42515a",
      alpha = 0.72
    ) +
    ggplot2::geom_boxplot(
      width = 0.13,
      outlier.shape = NA,
      color = "#17232b",
      fill = "#fffefb",
      alpha = 0.78,
      linewidth = 0.32
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "#6b747a",
      linewidth = 0.35,
      linetype = "dashed"
    )
  if (nrow(selected) > 0L) {
    plot <- plot + ggplot2::geom_point(
      data = selected,
      shape = 21,
      size = 2,
      stroke = 0.7,
      color = "#111820",
      fill = "#fffefb",
      position = ggplot2::position_jitter(width = 0.08, height = 0)
    )
  }
  plot +
    ggplot2::facet_wrap(~gene, ncol = 1L) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::labs(
      x = "Final cluster",
      y = "Log normalized expression",
      caption = paste(
        "Log normalized expression is centered and can be negative;",
        "detection uses raw counts.",
        "Outlined points are explicitly selected cells."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", color = "#17232b"),
      plot.caption = ggplot2::element_text(
        color = "#526b7b",
        hjust = 0,
        size = 8.5
      ),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

gene_pair_scope <- function(gene_data) {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    return(NULL)
  }
  included <- as.logical(gene_data$detected_1) |
    as.logical(gene_data$detected_2)
  total_n <- nrow(gene_data)
  included_n <- sum(included)
  list(
    included = included,
    total_n = total_n,
    included_n = included_n,
    excluded_n = total_n - included_n,
    excluded_pct = if (total_n > 0L) {
      100 * (total_n - included_n) / total_n
    } else {
      NA_real_
    }
  )
}

gene_pair_scope_ui <- function(scope) {
  if (is.null(scope)) {
    return(NULL)
  }
  excluded_pct <- trimws(formatC(
    scope$excluded_pct,
    digits = 3L,
    format = "fg"
  ))
  excluded_label <- if (scope$excluded_n == 1L) "cell" else "cells"
  htmltools::p(
    paste0(
      "Showing ",
      format(scope$included_n, big.mark = ","),
      " of ",
      format(scope$total_n, big.mark = ","),
      " cells detected by raw count for at least one gene. Excluded ",
      format(scope$excluded_n, big.mark = ","),
      " double-negative ",
      excluded_label,
      " (",
      excluded_pct,
      "%) to avoid the artificial diagonal created by their shared per-cell ",
      "log normalized expression centering offset."
    ),
    class = "supporting-copy gene-pair-scope",
    role = "status"
  )
}

make_gene_pair_plotly <- function(gene_data, bundle, source) {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    return(NULL)
  }
  scope <- gene_pair_scope(gene_data)
  data <- gene_data[scope$included, , drop = FALSE]
  if (nrow(data) == 0L) {
    return(NULL)
  }
  clusters <- sort(unique(as.character(data$cluster)))
  palette <- discrete_palette(bundle, "cluster", clusters)
  plot <- plotly::plot_ly(source = source)
  for (cluster in clusters) {
    trace_data <- data[as.character(data$cluster) == cluster, , drop = FALSE]
    hover <- paste0(
      "Cell: ",
      trace_data$cell_id,
      "<br>Cluster: ",
      cluster,
      "<br>",
      genes[[1L]],
      " Log normalized expression: ",
      formatC(trace_data$expression_1, digits = 4L, format = "fg"),
      " (Raw detected: ",
      ifelse(trace_data$detected_1, "yes", "no"),
      ")<br>",
      genes[[2L]],
      " Log normalized expression: ",
      formatC(trace_data$expression_2, digits = 4L, format = "fg"),
      " (Raw detected: ",
      ifelse(trace_data$detected_2, "yes", "no"),
      ")<br>Explicitly selected: ",
      ifelse(trace_data$selected, "yes", "no")
    )
    plot <- plotly::add_trace(
      plot,
      data = trace_data,
      x = ~expression_1,
      y = ~expression_2,
      key = ~cell_id,
      text = hover,
      hoverinfo = "text",
      type = "scattergl",
      mode = "markers",
      name = paste("Cluster", cluster),
      marker = list(
        color = unname(palette[[cluster]]),
        size = ifelse(trace_data$selected, 8, 5),
        symbol = ifelse(trace_data$selected, "diamond", "circle"),
        opacity = 0.72,
        line = list(
          color = "#111820",
          width = ifelse(trace_data$selected, 1.2, 0)
        )
      ),
      inherit = FALSE,
      showlegend = TRUE
    )
  }
  plot |>
    plotly::layout(
      xaxis = list(
        title = list(text = paste(genes[[1L]], "log normalized expression")),
        zeroline = TRUE,
        zerolinecolor = "#aab2b7"
      ),
      yaxis = list(
        title = list(text = paste(genes[[2L]], "log normalized expression")),
        zeroline = TRUE,
        zerolinecolor = "#aab2b7"
      ),
      legend = list(title = list(text = "Final cluster")),
      hovermode = "closest",
      margin = list(l = 62, r = 20, t = 18, b = 58)
    ) |>
    plotly::config(
      displaylogo = FALSE,
      toImageButtonOptions = list(filename = "espiviz-gene-pair-pflog")
    )
}

prepare_selection_snapshot <- function(gene_data) {
  genes <- attr(gene_data, "genes")
  explicit <- any(gene_data$selected)
  index <- if (explicit) which(gene_data$selected) else seq_len(nrow(gene_data))
  selected_n <- length(index)
  total_n <- nrow(gene_data)
  gene_rows <- lapply(seq_along(genes), function(gene_index) {
    expression <- gene_data[[paste0("expression_", gene_index)]][index]
    detected <- gene_data[[paste0("detected_", gene_index)]][index]
    data.frame(
      gene = genes[[gene_index]],
      mean_pflog = if (selected_n > 0L) mean(expression) else NA_real_,
      detected_n = sum(detected),
      detected_pct = if (selected_n > 0L) 100 * mean(detected) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  list(
    explicit_selection = explicit,
    selected_n = selected_n,
    total_n = total_n,
    selected_pct = if (total_n > 0L) 100 * selected_n / total_n else NA_real_,
    cluster_n = length(unique(as.character(gene_data$cluster[index]))),
    sample_n = length(unique(as.character(gene_data$sample[index]))),
    condition_n = length(unique(as.character(gene_data$condition[index]))),
    genes = do.call(rbind, gene_rows)
  )
}

selection_snapshot_ui <- function(snapshot) {
  cell_label <- if (isTRUE(snapshot$explicit_selection)) {
    paste0(
      format(snapshot$selected_n, big.mark = ","),
      " of ",
      format(snapshot$total_n, big.mark = ","),
      " cells"
    )
  } else {
    paste("All", format(snapshot$total_n, big.mark = ","), "cells")
  }
  htmltools::div(
    class = "selection-snapshot",
    `aria-live` = "polite",
    htmltools::div(
      class = "selection-snapshot-heading",
      htmltools::strong(cell_label),
      htmltools::span(
        paste0(
          trimws(formatC(snapshot$selected_pct, digits = 3L, format = "fg")),
          "%"
        )
      )
    ),
    htmltools::div(
      class = "selection-snapshot-meta",
      paste(
        snapshot$cluster_n,
        if (snapshot$cluster_n == 1L) "cluster" else "clusters",
        "·",
        snapshot$sample_n,
        if (snapshot$sample_n == 1L) "sample" else "samples"
      )
    ),
    htmltools::div(
      class = "selection-snapshot-genes",
      lapply(seq_len(nrow(snapshot$genes)), function(index) {
        row <- snapshot$genes[index, , drop = FALSE]
        htmltools::div(
          class = "selection-snapshot-gene",
          htmltools::strong(row$gene),
          htmltools::span(
            paste0(
              "Mean log normalized expression: ",
              trimws(formatC(row$mean_pflog, digits = 3L, format = "fg"))
            )
          ),
          htmltools::span(paste0(
            format(row$detected_n, big.mark = ","),
            " of ",
            format(snapshot$selected_n, big.mark = ","),
            " detected by raw count (",
            trimws(formatC(row$detected_pct, digits = 3L, format = "fg")),
            "%)"
          ))
        )
      })
    )
  )
}

prepare_comparison_table <- function(comparison, explicit_selection) {
  has_reference <- isTRUE(explicit_selection) &&
    nrow(comparison) > 0L &&
    any(comparison$remaining_n > 0L, na.rm = TRUE)
  if (has_reference) {
    return(comparison[, c(
      "gene",
      "selected_mean",
      "selected_detected_pct",
      "remaining_mean",
      "remaining_detected_pct",
      "mean_difference",
      "detection_pp_difference"
    ), drop = FALSE])
  }
  data.frame(
    gene = as.character(comparison$gene),
    mean_expression = as.numeric(comparison$selected_mean),
    detected_pct = as.numeric(comparison$selected_detected_pct),
    stringsAsFactors = FALSE
  )
}
