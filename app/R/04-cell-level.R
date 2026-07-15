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

scale_expression_strength <- function(values) {
  values <- as.numeric(values)
  result <- numeric(length(values))
  index <- which(is.finite(values))
  if (length(index) == 0L) {
    return(result)
  }
  limits <- range(values[index])
  if (diff(limits) <= .Machine$double.eps^0.5) {
    result[index] <- as.numeric(limits[[1L]] > 0)
    return(result)
  }
  result[index] <- (values[index] - limits[[1L]]) / diff(limits)
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

prepare_umap_expression_blend_data <- function(gene_data) {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    stop(
      "Two plot genes are required for an expression-blended UMAP.",
      call. = FALSE
    )
  }
  data <- gene_data
  data$strength_1 <- scale_expression_strength(data$expression_1)
  data$strength_2 <- scale_expression_strength(data$expression_2)
  data$blend_strength <- data$strength_1 + data$strength_2
  data$blend_color <- blend_rgb(data$strength_1, data$strength_2)
  data <- data[order(data$blend_strength, data$cell_id), , drop = FALSE]
  rownames(data) <- NULL
  attr(data, "genes") <- genes
  data
}

prepare_umap_detection_data <- function(gene_data) {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    stop(
      "Two plot genes are required for a two-gene detection UMAP.",
      call. = FALSE
    )
  }
  data <- gene_data
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
  color_key <- ifelse(
    data$blend_class == "Neither detected",
    "Neither detected",
    ifelse(
      data$blend_class == paste(genes[[1L]], "detected"),
      "Gene 1",
      ifelse(
        data$blend_class == paste(genes[[2L]], "detected"),
        "Gene 2",
        "Both detected"
      )
    )
  )
  data$blend_color <- unname(blend_palette()[color_key])
  data$blend_strength <- as.integer(data$detected_1) +
    as.integer(data$detected_2)
  data <- data[order(data$blend_strength, data$cell_id), , drop = FALSE]
  rownames(data) <- NULL
  attr(data, "genes") <- genes
  data
}

blend_legend_ui <- function(genes, mode = c("expression", "detection")) {
  genes <- compact_character(genes)
  if (length(genes) != 2L) {
    return(NULL)
  }
  mode <- match.arg(mode)
  palette <- blend_palette()
  labels <- if (identical(mode, "expression")) {
    c(
      paste("Low", genes[[1L]], "+ low", genes[[2L]]),
      paste("High", genes[[1L]]),
      paste("High", genes[[2L]]),
      paste("High", genes[[1L]], "+ high", genes[[2L]])
    )
  } else {
    c(
      "Neither detected",
      paste(genes[[1L]], "detected"),
      paste(genes[[2L]], "detected"),
      "Both detected"
    )
  }
  colors <- unname(palette[c(
    "Neither detected",
    "Gene 1",
    "Gene 2",
    "Both detected"
  )])
  htmltools::div(
    class = "blend-legend",
    role = "list",
    `aria-label` = paste(
      "Two-gene",
      mode,
      "blend for",
      paste(genes, collapse = " and ")
    ),
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

gene_pair_scope <- function(gene_data, group = "all") {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    return(NULL)
  }
  group <- compact_character(group)
  group <- if (length(group) == 0L) "all" else group[[1L]]
  group_label <- "All cells"
  group_rows <- rep(TRUE, nrow(gene_data))
  if (!identical(group, "all")) {
    group_rows <- as.character(gene_data$cluster) == group
    group_label <- paste("Cluster", group)
  }
  included <- group_rows & (
    as.logical(gene_data$detected_1) |
      as.logical(gene_data$detected_2)
  )
  total_n <- sum(group_rows)
  included_n <- sum(included)
  list(
    included = included,
    group = group,
    group_label = group_label,
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
  group_label <- scope$group_label %||% "All cells"
  group_suffix <- if (identical(group_label, "All cells")) {
    ""
  } else {
    paste0(" in ", group_label)
  }
  htmltools::p(
    paste0(
      "Showing ",
      format(scope$included_n, big.mark = ","),
      " of ",
      format(scope$total_n, big.mark = ","),
      " cells",
      group_suffix,
      " detected by raw count for at least one gene. Excluded ",
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

gene_pair_group_data <- function(gene_data, group = "all") {
  scope <- gene_pair_scope(gene_data, group = group)
  if (is.null(scope)) {
    return(list(
      data = gene_data[0L, , drop = FALSE],
      group = "all",
      group_label = "All cells"
    ))
  }
  data <- gene_data[scope$included, , drop = FALSE]
  data <- data[
    is.finite(data$expression_1) & is.finite(data$expression_2),
    ,
    drop = FALSE
  ]
  list(
    data = data,
    group = scope$group,
    group_label = scope$group_label
  )
}

empty_gene_pair_loess <- function() {
  data.frame(
    x = numeric(),
    fit = numeric(),
    lower = numeric(),
    upper = numeric(),
    group_label = character(),
    cell_count = integer(),
    stringsAsFactors = FALSE
  )
}

prepare_gene_pair_loess <- function(
  gene_data,
  group = "all",
  span = 0.75,
  level = 0.95,
  grid_n = 100L
) {
  grouped <- gene_pair_group_data(gene_data, group)
  data <- grouped$data
  if (
    nrow(data) < 8L ||
      length(unique(data$expression_1)) < 4L ||
      length(unique(data$expression_2)) < 2L
  ) {
    return(empty_gene_pair_loess())
  }
  grid_n <- suppressWarnings(as.integer(grid_n[[1L]]))
  if (is.na(grid_n) || grid_n < 20L) {
    grid_n <- 100L
  }
  model_data <- data.frame(
    x = as.numeric(data$expression_1),
    y = as.numeric(data$expression_2)
  )
  model <- tryCatch(
    suppressWarnings(stats::loess(
      y ~ x,
      data = model_data,
      span = span,
      degree = 2L,
      control = stats::loess.control(surface = "direct")
    )),
    error = function(error) NULL
  )
  if (is.null(model)) {
    return(empty_gene_pair_loess())
  }
  x <- seq(min(model_data$x), max(model_data$x), length.out = grid_n)
  predicted <- tryCatch(
    suppressWarnings(stats::predict(
      model,
      newdata = data.frame(x = x),
      se = TRUE
    )),
    error = function(error) NULL
  )
  if (is.null(predicted) || !is.list(predicted)) {
    return(empty_gene_pair_loess())
  }
  degrees_freedom <- as.numeric(predicted$df %||% NA_real_)
  critical <- if (
    length(degrees_freedom) == 1L &&
      is.finite(degrees_freedom) &&
      degrees_freedom > 0
  ) {
    stats::qt((1 + level) / 2, df = degrees_freedom)
  } else {
    stats::qnorm((1 + level) / 2)
  }
  fit <- as.numeric(predicted$fit)
  standard_error <- as.numeric(predicted$se.fit)
  keep <- is.finite(x) & is.finite(fit) & is.finite(standard_error)
  if (!any(keep)) {
    return(empty_gene_pair_loess())
  }
  data.frame(
    x = x[keep],
    fit = fit[keep],
    lower = fit[keep] - critical * standard_error[keep],
    upper = fit[keep] + critical * standard_error[keep],
    group_label = grouped$group_label,
    cell_count = as.integer(nrow(data)),
    stringsAsFactors = FALSE
  )
}

prepare_gene_pair_density <- function(
  gene_data,
  group = "all",
  grid_n = 80L
) {
  grouped <- gene_pair_group_data(gene_data, group)
  data <- grouped$data
  if (
    nrow(data) < 3L ||
      length(unique(data$expression_1)) < 2L ||
      length(unique(data$expression_2)) < 2L
  ) {
    return(NULL)
  }
  grid_n <- suppressWarnings(as.integer(grid_n[[1L]]))
  if (is.na(grid_n) || grid_n < 20L) {
    grid_n <- 80L
  }
  density <- tryCatch(
    MASS::kde2d(
      x = as.numeric(data$expression_1),
      y = as.numeric(data$expression_2),
      n = grid_n
    ),
    error = function(error) NULL
  )
  if (is.null(density) || any(!is.finite(density$z))) {
    return(NULL)
  }
  list(
    x = density$x,
    y = density$y,
    z = density$z,
    group = grouped$group,
    group_label = grouped$group_label,
    cell_count = as.integer(nrow(data))
  )
}

plotly_alpha_color <- function(color, alpha) {
  rgb <- grDevices::col2rgb(color)
  sprintf(
    "rgba(%d,%d,%d,%.3f)",
    rgb[[1L]],
    rgb[[2L]],
    rgb[[3L]],
    alpha
  )
}

gene_pair_plotly_layout <- function(plot, genes, legend_title = NULL) {
  plotly::layout(
    plot,
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
    legend = if (is.null(legend_title)) NULL else {
      list(title = list(text = legend_title))
    },
    hovermode = "closest",
    margin = list(l = 62, r = 20, t = 18, b = 58)
  )
}

make_gene_pair_density_plotly <- function(
  gene_data,
  source,
  group = "all"
) {
  genes <- attr(gene_data, "genes")
  density <- prepare_gene_pair_density(gene_data, group = group)
  if (is.null(density)) {
    empty <- plotly::plot_ly(
      x = numeric(),
      y = numeric(),
      type = "scatter",
      mode = "markers",
      source = source,
      hoverinfo = "skip",
      showlegend = FALSE
    )
    empty <- gene_pair_plotly_layout(empty, genes)
    empty <- plotly::layout(
      empty,
      annotations = list(list(
        text = "Density unavailable for this group.",
        x = 0.5,
        y = 0.5,
        xref = "paper",
        yref = "paper",
        showarrow = FALSE
      ))
    )
    return(plotly::config(empty, displaylogo = FALSE))
  }
  plot <- plotly::plot_ly(
    x = density$x,
    y = density$y,
    z = t(density$z),
    type = "contour",
    source = source,
    name = paste("Density —", density$group_label),
    colors = c("#f5f2ec", "#d69ab2", "#9d2857", "#17232b"),
    contours = list(coloring = "heatmap", showlabels = FALSE),
    colorbar = list(title = list(text = "Cell density")),
    hovertemplate = paste0(
      genes[[1L]],
      ": %{x:.3f}<br>",
      genes[[2L]],
      ": %{y:.3f}<br>Density: %{z:.3g}<extra></extra>"
    )
  )
  gene_pair_plotly_layout(plot, genes) |>
    plotly::config(
      displaylogo = FALSE,
      toImageButtonOptions = list(filename = "espiviz-gene-pair-density")
    )
}

make_gene_pair_plotly <- function(
  gene_data,
  bundle,
  source,
  display = c("scatter", "density"),
  loess_group = "all",
  density_group = "all"
) {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    return(NULL)
  }
  display <- match.arg(display)
  if (identical(display, "density")) {
    return(make_gene_pair_density_plotly(
      gene_data,
      source = source,
      group = density_group
    ))
  }
  scope <- gene_pair_scope(gene_data)
  data <- gene_data[scope$included, , drop = FALSE]
  if (nrow(data) == 0L) {
    return(NULL)
  }
  clusters <- sort(unique(as.character(data$cluster)))
  palette <- discrete_palette(bundle, "cluster", clusters)
  plot <- plotly::plot_ly(source = source)
  trend <- prepare_gene_pair_loess(gene_data, group = loess_group)
  trend_color <- if (identical(loess_group, "all")) {
    "#17232b"
  } else {
    candidate <- unname(palette[as.character(loess_group)])
    if (length(candidate) == 1L && !is.na(candidate) && nzchar(candidate)) {
      candidate
    } else {
      "#17232b"
    }
  }
  if (nrow(trend) > 0L) {
    plot <- plotly::add_ribbons(
      plot,
      data = trend,
      x = ~x,
      ymin = ~lower,
      ymax = ~upper,
      name = "95% confidence ribbon",
      fillcolor = plotly_alpha_color(trend_color, 0.18),
      line = list(color = "transparent"),
      hoverinfo = "skip",
      showlegend = FALSE,
      inherit = FALSE
    )
  }
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
  if (nrow(trend) > 0L) {
    plot <- plotly::add_lines(
      plot,
      data = trend,
      x = ~x,
      y = ~fit,
      name = paste("Loess trend —", trend$group_label[[1L]]),
      line = list(color = trend_color, width = 3),
      hovertemplate = paste0(
        "Loess trend — ",
        trend$group_label[[1L]],
        "<br>",
        genes[[1L]],
        ": %{x:.3f}<br>",
        genes[[2L]],
        ": %{y:.3f}<extra></extra>"
      ),
      inherit = FALSE,
      showlegend = FALSE
    )
  }
  plot <- gene_pair_plotly_layout(
    plot,
    genes,
    legend_title = "Final cluster"
  )
  if (nrow(trend) == 0L) {
    trend_group <- gene_pair_group_data(gene_data, group = loess_group)
    plot <- plotly::layout(
      plot,
      annotations = list(list(
        text = paste(
          "Loess trend unavailable for",
          paste0(trend_group$group_label, ".")
        ),
        x = 0.5,
        y = 1,
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "top",
        showarrow = FALSE
      ))
    )
  }
  plotly::config(
    plot,
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
