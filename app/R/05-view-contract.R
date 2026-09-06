# Public view state, result scope, and portable figure metadata.
count_label <- function(n, noun = "gene") {
  paste(
    format(n, big.mark = ",", trim = TRUE),
    paste0(noun, if (n == 1L) "" else "s")
  )
}

analysis_gene_scope <- function(
  current,
  second = character(),
  saved = character()
) {
  genes <- compact_character(c(current, second, saved))
  genes[!duplicated(casefold_key(genes))]
}

blend_explanation <- function() {
  paste(
    "Low/high is relative log normalized expression among cells detecting each gene.",
    "Genes are scaled independently; this is not an expression ratio.",
    "Undetected genes contribute no color and double-negative cells stay neutral.",
    "Constant detected values use mid intensity."
  )
}

wrap_caption <- function(text, width = 100L) {
  paste(
    unlist(lapply(
      strsplit(paste(text, collapse = " "), "\n")[[1L]],
      strwrap,
      width = width
    )),
    collapse = "\n"
  )
}

view_option_defaults <- function() {
  list(
    secondary_gene = "",
    color_by = "expression",
    explore_results = "By cluster",
    gene_page = 1L,
    cluster_plot_type = "auto",
    comparison_plot_type = "auto",
    sample_plot_type = "auto",
    condition_plot_type = "auto",
    color_scale = "auto",
    gene_pair_display = "scatter",
    gene_pair_loess_group = "all",
    gene_pair_density_group = "all",
    marker_cluster = "1"
  )
}

resolved_plot_type <- function(value, gene_count) {
  if (is.null(value) || identical(value, "auto")) {
    if (gene_count > 6L) "dot" else "violin"
  } else {
    value
  }
}

result_uses_gene_page <- function(tab) {
  tab %in% c("By cluster", "Comparison", "By sample", "Pooled condition")
}

bundle_identity <- function(bundle) {
  list(
    data_version = as.character(bundle$data_version),
    sha256 = manifest_sha256(attr(bundle, "data_manifest") %||% list())
  )
}

selection_description <- function(bundle, selected, context = FALSE) {
  n <- length(intersect(bundle$cells$cell_id, selected))
  if (context) {
    return(paste0(
      "All ",
      count_label(nrow(bundle$cells), "cell"),
      if (n > 0L) {
        paste0(" · ", count_label(n, "cell"), " selected in Explore")
      } else {
        ""
      }
    ))
  }
  if (n == 0L) {
    "All cells — no selection applied"
  } else {
    paste(
      count_label(n, "selected cell"),
      "of",
      format(nrow(bundle$cells), big.mark = ",")
    )
  }
}

capture_view_state <- function(bundle, state, view = "explore") {
  selected <- as.character(state$selected_cells())
  clusters <- unique(as.character(bundle$cells$cluster))
  matched <- clusters[vapply(
    clusters,
    function(cluster) {
      length(selected) > 0L &&
        setequal(selected, cell_ids_for_cluster(bundle, cluster))
    },
    logical(1L)
  )]
  selection <- if (length(selected) == 0L) {
    list(type = "all")
  } else if (length(matched)) {
    list(type = "cluster", cluster = matched[[1L]])
  } else {
    list(type = "cells", cells = selected)
  }
  list(
    state_version = 1L,
    bundle = bundle_identity(bundle),
    view = view,
    gene = state$active_gene(),
    genes = state$gene_set(),
    gene_set_name = state$gene_set_name(),
    pathway = state$active_pathway(),
    selection = selection,
    options = state$view_options()
  )
}

validate_view_state <- function(value, bundle) {
  if (!is.list(value) || !identical(as.integer(value$state_version), 1L)) {
    stop("This is not a supported ESPIviz view file.", call. = FALSE)
  }
  identity <- bundle_identity(bundle)
  if (
    !identical(
      as.character(value$bundle$data_version),
      identity$data_version
    ) ||
      (!is.na(identity$sha256) &&
        !identical(value$bundle$sha256, identity$sha256))
  ) {
    stop(
      "This view uses a different data bundle. Open it with its original ESPIviz data version.",
      call. = FALSE
    )
  }
  parsed <- parse_gene_input(
    c(value$gene, value$genes, value$options$secondary_gene),
    bundle_gene_names(bundle)
  )
  if (length(parsed$missing)) {
    stop("The saved view contains unrecognized genes.", call. = FALSE)
  }
  if (length(value$gene) != 1L || !value$gene %in% bundle_gene_names(bundle)) {
    stop("The saved view has no valid current gene.", call. = FALSE)
  }
  if (
    length(value$view) != 1L ||
      !value$view %in% c("explore", "de", "pathways", "about")
  ) {
    stop("The saved view has an unknown page.", call. = FALSE)
  }
  defaults <- view_option_defaults()
  options <- utils::modifyList(defaults, value$options %||% list())
  if (length(options$secondary_gene) != 1L) {
    stop("The saved view must have at most one second gene.", call. = FALSE)
  }
  choices <- list(
    color_by = c("expression", "detection", "condition", "cluster"),
    explore_results = c(
      "By cluster",
      "Comparison",
      "Cell-level",
      "By sample",
      "Pooled condition",
      "Cluster markers",
      "Gene set"
    ),
    color_scale = c("auto", "shared"),
    gene_pair_display = c("scatter", "density"),
    gene_pair_loess_group = c(
      "none",
      "all",
      unique(as.character(bundle$cells$cluster))
    ),
    gene_pair_density_group = c(
      "all",
      unique(as.character(bundle$cells$cluster))
    ),
    marker_cluster = unique(as.character(bundle$cells$cluster))
  )
  for (key in c(
    "cluster_plot_type",
    "comparison_plot_type",
    "sample_plot_type",
    "condition_plot_type"
  )) {
    choices[[key]] <- c("auto", "dot", "violin")
  }
  for (key in names(choices)) {
    if (length(options[[key]]) != 1L || !options[[key]] %in% choices[[key]]) {
      stop(
        "The saved view contains an unsupported plot setting.",
        call. = FALSE
      )
    }
  }
  page <- suppressWarnings(as.integer(options$gene_page))
  if (length(page) != 1L || is.na(page) || page < 1L) {
    page <- 1L
  }
  options$gene_page <- page
  selection <- value$selection
  if (identical(selection$type, "cluster")) {
    if (!selection$cluster %in% unique(as.character(bundle$cells$cluster))) {
      stop("Unknown saved cluster.", call. = FALSE)
    }
    selected <- cell_ids_for_cluster(bundle, selection$cluster)
  } else if (identical(selection$type, "cells")) {
    selected <- compact_character(selection$cells)
    if (!all(selected %in% bundle$cells$cell_id)) {
      stop("Unknown saved cells.", call. = FALSE)
    }
  } else if (identical(selection$type, "all")) {
    selected <- character()
  } else {
    stop("The saved view has an invalid cell selection.", call. = FALSE)
  }
  value$options <- options[names(defaults)]
  value$selected_cells <- selected
  value$genes <- compact_character(value$genes)
  value
}

restore_view_state <- function(value, bundle, state) {
  value <- validate_view_state(value, bundle)
  state$gene_set(value$genes)
  state$gene_set_name(value$gene_set_name)
  state$active_gene(value$gene)
  state$selected_cells(value$selected_cells)
  if (
    length(value$pathway) == 1L && value$pathway %in% bundle$pathways$pathway_id
  ) {
    state$active_pathway(value$pathway)
  }
  state$view_options(value$options)
  state$restore(value$options)
  invisible(value)
}

view_link_query <- function(value, max_bytes = 6000L) {
  if (identical(value$selection$type, "cells")) {
    return(NULL)
  }
  json <- jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", na = "null")
  query <- paste0("?view_state=", utils::URLencode(json, reserved = TRUE))
  if (nchar(query, type = "bytes") > max_bytes) NULL else query
}

figure_context <- function(
  plot,
  bundle,
  title,
  genes = character(),
  selected = character(),
  context = FALSE,
  note = "",
  width = 9,
  scope_text = NULL
) {
  identity <- bundle_identity(bundle)
  gene_note <- if (length(genes) <= 8L) {
    paste(genes, collapse = ", ")
  } else {
    count_label(length(genes))
  }
  subtitle <- paste(
    compact_character(c(
      gene_note,
      scope_text %||% selection_description(bundle, selected, context)
    )),
    collapse = " · "
  )
  caption <- paste(
    compact_character(c(
      plot$labels$caption,
      note,
      paste0("ESPIviz · Data v", identity$data_version)
    )),
    collapse = " "
  )
  plot +
    ggplot2::labs(
      title = wrap_caption(title, 75L),
      subtitle = wrap_caption(subtitle, 95L),
      caption = wrap_caption(caption, as.integer(width * 12))
    ) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15,
        margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10,
        color = "#58646a",
        margin = ggplot2::margin(b = 14)
      ),
      plot.caption = ggplot2::element_text(
        size = 9,
        color = "#58646a",
        hjust = 0,
        margin = ggplot2::margin(t = 14)
      ),
      plot.margin = ggplot2::margin(16, 18, 16, 18)
    )
}

write_figure <- function(plot, file, format = "png", width = 9, height = 7) {
  if (is.null(plot)) {
    stop("Choose a result before downloading its figure.", call. = FALSE)
  }
  ggplot2::ggsave(
    file,
    plot = plot,
    device = if (format == "pdf") grDevices::cairo_pdf else "png",
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white",
    limitsize = FALSE
  )
}

make_gene_pair_ggplot <- function(
  gene_data,
  bundle,
  display = "scatter",
  trend_group = "all",
  density_group = "all"
) {
  genes <- attr(gene_data, "genes")
  if (length(genes) != 2L) {
    return(NULL)
  }
  scope <- gene_pair_scope(
    gene_data,
    if (display == "density") density_group else "all"
  )
  data <- gene_data[scope$included, , drop = FALSE]
  if (nrow(data) == 0L) {
    return(NULL)
  }
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = expression_1, y = expression_2)
  )
  if (display == "density") {
    density <- prepare_gene_pair_density(gene_data, density_group)
    if (is.null(density)) {
      return(NULL)
    }
    grid <- expand.grid(expression_1 = density$x, expression_2 = density$y)
    grid$density <- as.vector(density$z)
    plot <- plot +
      ggplot2::geom_contour_filled(data = grid, ggplot2::aes(z = density)) +
      ggplot2::labs(fill = "Probability density")
  } else {
    plot <- plot +
      ggplot2::geom_point(
        ggplot2::aes(color = factor(cluster)),
        size = 1.4,
        alpha = 0.55
      ) +
      ggplot2::scale_color_manual(
        values = discrete_palette(
          bundle,
          "cluster",
          sort(unique(as.character(data$cluster)))
        )
      ) +
      ggplot2::labs(color = "Final cluster")
    selected <- data[data$selected, , drop = FALSE]
    if (nrow(selected)) {
      plot <- plot +
        ggplot2::geom_point(data = selected, shape = 23, fill = NA, size = 2.5)
    }
    if (trend_group != "none") {
      trend <- prepare_gene_pair_loess(gene_data, trend_group)
      if (nrow(trend)) {
        plot <- plot +
          ggplot2::geom_ribbon(
            data = trend,
            ggplot2::aes(x = x, ymin = lower, ymax = upper),
            inherit.aes = FALSE,
            fill = "#14764f",
            alpha = 0.15
          ) +
          ggplot2::geom_line(
            data = trend,
            ggplot2::aes(x = x, y = fit),
            inherit.aes = FALSE,
            color = "#14764f",
            linewidth = 0.8
          )
      }
    }
  }
  plot +
    ggplot2::labs(
      x = paste0(genes[[1L]], "\nlog normalized expression"),
      y = paste0(genes[[2L]], "\nlog normalized expression"),
      caption = paste(
        "Double-negative cells are excluded.",
        if (display == "density") {
          paste(
            "Smoothed probability density for",
            density$group_label,
            "using",
            density$cell_count,
            "cells."
          )
        } else if (trend_group != "none") {
          "LOESS and 95% confidence interval describe cells; they do not estimate biological-replicate uncertainty."
        } else {
          "No trend fitted."
        }
      )
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}
