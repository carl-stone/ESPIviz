prepare_pathway_plot_data <- function(pathways) {
  required <- c(
    "pathway_id",
    "label",
    "source",
    "direction",
    "p_adjust",
    "score",
    "gene_count"
  )
  if (!is.data.frame(pathways) || !all(required %in% names(pathways))) {
    stop("Pathway plot data are missing required columns.", call. = FALSE)
  }
  if (nrow(pathways) == 0L) {
    stop("At least one enrichment result is required.", call. = FALSE)
  }

  data <- pathways
  data$label <- as.character(data$label)
  data$pathway_id <- as.character(data$pathway_id)
  data$p_adjust <- as.numeric(data$p_adjust)
  data$score <- as.numeric(data$score)
  data$gene_count <- as.numeric(data$gene_count)

  if (anyNA(data$pathway_id) || any(!nzchar(data$pathway_id))) {
    stop("Pathway IDs must be complete.", call. = FALSE)
  }
  if (anyDuplicated(data$pathway_id)) {
    stop("Pathway IDs must be unique.", call. = FALSE)
  }
  if (anyNA(data$label) || any(!nzchar(data$label))) {
    stop("Pathway labels must be complete.", call. = FALSE)
  }
  if (
    anyNA(data$score) ||
      any(!is.finite(data$score)) ||
      anyNA(data$p_adjust) ||
      any(!is.finite(data$p_adjust)) ||
      any(data$p_adjust < 0) ||
      any(data$p_adjust > 1) ||
      anyNA(data$gene_count) ||
      any(!is.finite(data$gene_count)) ||
      any(data$gene_count < 0)
  ) {
    stop(
      "Pathway scores, adjusted P values, and gene counts are invalid.",
      call. = FALSE
    )
  }

  source_levels <- c("GSEA", "ORA")
  source_values <- as.character(data$source)
  if (anyNA(source_values) || any(!source_values %in% source_levels)) {
    stop("Pathway methods must be GSEA or ORA.", call. = FALSE)
  }

  direction_levels <- c("Control", "E-Stim")
  direction_values <- as.character(data$direction)
  if (
    anyNA(direction_values) ||
      any(!direction_values %in% direction_levels)
  ) {
    stop("Pathway directions must be Control or E-Stim.", call. = FALSE)
  }

  data$source <- factor(source_values, levels = source_levels, ordered = TRUE)
  data$direction <- factor(
    direction_values,
    levels = direction_levels,
    ordered = TRUE
  )
  repeated_labels <- duplicated(data$label) |
    duplicated(data$label, fromLast = TRUE)
  data$plot_label <- ifelse(
    repeated_labels,
    paste(data$label, source_values, direction_values, sep = " · "),
    data$label
  )
  repeated_plot_labels <- duplicated(data$plot_label) |
    duplicated(data$plot_label, fromLast = TRUE)
  data$plot_label[repeated_plot_labels] <- paste(
    data$plot_label[repeated_plot_labels],
    data$pathway_id[repeated_plot_labels],
    sep = " · "
  )
  data$score_label <- ifelse(source_values == "GSEA", "NES", "Fold enrichment")
  data$baseline <- ifelse(source_values == "GSEA", 0, 1)
  data$count_label <- ifelse(
    source_values == "GSEA",
    "Source gene-set size",
    "Overlapping genes"
  )
  data$evidence_strength <- -log10(pmax(
    data$p_adjust,
    .Machine$double.xmin
  ))
  data$source_order <- match(source_values, source_levels)
  data$panel_label <- factor(
    paste(source_values, data$score_label, sep = " · "),
    levels = c("GSEA · NES", "ORA · Fold enrichment"),
    ordered = TRUE
  )
  data$pathway_key <- data$pathway_id
  data$hover_text <- paste0(
    htmltools::htmlEscape(data$label),
    "<br>Direction: ",
    direction_values,
    "<br>Method: ",
    source_values,
    "<br>",
    data$score_label,
    ": ",
    trimws(formatC(data$score, digits = 4L, format = "g")),
    "<br>Adjusted P: ",
    trimws(formatC(data$p_adjust, digits = 4L, format = "g")),
    "<br>",
    data$count_label,
    ": ",
    format(data$gene_count, big.mark = ",", scientific = FALSE, trim = TRUE)
  )
  data
}

top_pathway_results <- function(pathways, n_per_direction = 10L) {
  if (
    length(n_per_direction) != 1L ||
      is.na(n_per_direction) ||
      n_per_direction < 1 ||
      n_per_direction != as.integer(n_per_direction)
  ) {
    stop("n_per_direction must be one positive integer.", call. = FALSE)
  }

  data <- prepare_pathway_plot_data(pathways)
  data <- data[
    order(
      data$source,
      data$direction,
      data$p_adjust,
      -abs(data$score),
      data$label
    ),
    ,
    drop = FALSE
  ]
  groups <- interaction(data$source, data$direction, drop = TRUE)
  group_rank <- ave(seq_len(nrow(data)), groups, FUN = seq_along)
  selected_ids <- data$pathway_id[group_rank <= as.integer(n_per_direction)]

  pathways[
    match(selected_ids, as.character(pathways$pathway_id)),
    ,
    drop = FALSE
  ]
}

displayed_pathway_results <- function(
  pathways,
  active_pathway,
  n_per_direction = 10L
) {
  displayed <- top_pathway_results(pathways, n_per_direction)
  active_pathway <- as.character(active_pathway %||% character())
  active_pathway <- active_pathway[
    !is.na(active_pathway) & nzchar(active_pathway)
  ]
  if (
    length(active_pathway) == 0L ||
      active_pathway[[1L]] %in% as.character(displayed$pathway_id)
  ) {
    return(displayed)
  }

  active_index <- match(
    active_pathway[[1L]],
    as.character(pathways$pathway_id)
  )
  if (is.na(active_index)) {
    return(displayed)
  }

  rbind(displayed, pathways[active_index, , drop = FALSE])
}

pathway_direction_palette <- function() {
  c("Control" = "#2166AC", "E-Stim" = "#B52865")
}

pathway_direction_symbols <- function() {
  c("Control" = "circle", "E-Stim" = "diamond")
}

pathway_marker_sizes <- function(evidence_strength, size_range = c(10, 25)) {
  diameter_scale <- sqrt(pmax(as.numeric(evidence_strength), 0))
  limits <- base::range(diameter_scale, finite = TRUE)
  if (diff(limits) == 0) {
    return(rep(mean(size_range), length(diameter_scale)))
  }
  size_range[[1L]] +
    diff(size_range) * (diameter_scale - limits[[1L]]) / diff(limits)
}

make_pathway_plotly <- function(pathways, active_pathway, source) {
  data <- prepare_pathway_plot_data(pathways)
  data$active <- data$pathway_id == active_pathway
  data$marker_size <- pathway_marker_sizes(data$evidence_strength)
  methods <- levels(data$source)
  methods <- methods[methods %in% as.character(data$source)]
  direction_palette <- pathway_direction_palette()
  direction_symbols <- pathway_direction_symbols()
  direction_levels <- names(direction_palette)
  legend_seen <- character()
  panel_data <- vector("list", length(methods))
  panels <- vector("list", length(methods))

  for (method_index in seq_along(methods)) {
    method <- methods[[method_index]]
    method_data <- data[
      as.character(data$source) == method,
      ,
      drop = FALSE
    ]
    method_data <- method_data[
      order(method_data$score, method_data$label),
      ,
      drop = FALSE
    ]
    method_data$plot_label <- factor(
      method_data$plot_label,
      levels = unique(method_data$plot_label),
      ordered = TRUE
    )
    panel_data[[method_index]] <- method_data

    panel <- plotly::plot_ly(source = source)
    for (direction in direction_levels) {
      direction_data <- method_data[
        as.character(method_data$direction) == direction,
        ,
        drop = FALSE
      ]
      if (nrow(direction_data) == 0L) next

      show_direction_legend <- !direction %in% legend_seen
      legend_seen <- unique(c(legend_seen, direction))
      panel <- plotly::add_markers(
        panel,
        data = direction_data,
        x = ~score,
        y = ~plot_label,
        key = ~pathway_key,
        customdata = ~pathway_id,
        text = ~hover_text,
        hoverinfo = "text",
        name = direction,
        legendgroup = direction,
        showlegend = show_direction_legend,
        marker = list(
          size = direction_data$marker_size,
          sizemode = "diameter",
          color = unname(direction_palette[[direction]]),
          symbol = unname(direction_symbols[[direction]]),
          opacity = 0.9,
          line = list(
            color = ifelse(direction_data$active, "#111827", "#FFFFFF"),
            width = ifelse(direction_data$active, 3, 0.8)
          )
        )
      )
    }

    panels[[method_index]] <- panel |>
      plotly::layout(
        xaxis = list(
          title = list(text = method_data$score_label[[1L]]),
          automargin = TRUE,
          zeroline = FALSE
        ),
        yaxis = list(
          title = list(text = method),
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = as.character(method_data$plot_label)
        )
      )
  }

  if (length(panels) == 1L) {
    plot <- panels[[1L]]
  } else {
    heights <- pmax(vapply(panel_data, nrow, integer(1L)), 2L)
    plot <- do.call(
      plotly::subplot,
      c(
        panels,
        list(
          nrows = length(panels),
          shareX = FALSE,
          shareY = FALSE,
          titleX = TRUE,
          titleY = TRUE,
          heights = heights / sum(heights),
          margin = 0.09
        )
      )
    )
  }

  baseline_shapes <- lapply(seq_along(panel_data), function(index) {
    axis_suffix <- if (index == 1L) "" else as.character(index)
    list(
      type = "line",
      xref = paste0("x", axis_suffix),
      yref = paste0("y", axis_suffix),
      x0 = panel_data[[index]]$baseline[[1L]],
      x1 = panel_data[[index]]$baseline[[1L]],
      y0 = -0.5,
      y1 = nrow(panel_data[[index]]) - 0.5,
      line = list(color = "#6B7280", width = 1, dash = "dot"),
      layer = "below"
    )
  })

  plot <- plot |>
    plotly::layout(
      legend = list(
        title = list(text = "Direction"),
        orientation = "h",
        x = 0,
        y = 1.08
      ),
      margin = list(l = 30, r = 20, t = 55, b = 45),
      hovermode = "closest"
    )
  # subplot() initializes an empty shape list that its lazy layout merge keeps.
  # Assigning after the merge preserves one method-specific neutral baseline per
  # panel without adding non-pathway click targets.
  plot$x$layout$shapes <- baseline_shapes

  plot |>
    plotly::config(displaylogo = FALSE) |>
    plotly::event_register("plotly_click")
}

make_pathway_ggplot <- function(pathways, active_pathway) {
  data <- prepare_pathway_plot_data(pathways)
  data$active <- data$pathway_id == active_pathway
  data <- data[
    order(data$source_order, data$score, data$label),
    ,
    drop = FALSE
  ]
  data$plot_label <- factor(
    data$plot_label,
    levels = unique(data$plot_label),
    ordered = TRUE
  )
  baselines <- unique(data[c("panel_label", "baseline")])

  ggplot2::ggplot(data, ggplot2::aes(x = score, y = plot_label)) +
    ggplot2::geom_vline(
      data = baselines,
      ggplot2::aes(xintercept = baseline),
      inherit.aes = FALSE,
      color = "#6B7280",
      linewidth = 0.45,
      linetype = "dotted"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = direction,
        shape = direction,
        size = evidence_strength
      ),
      color = "white",
      stroke = 0.6,
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = data[data$active, , drop = FALSE],
      ggplot2::aes(
        fill = direction,
        shape = direction,
        size = evidence_strength
      ),
      color = "#111827",
      stroke = 1.2,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(panel_label),
      ncol = 1,
      scales = "free"
    ) +
    ggplot2::scale_fill_manual(
      name = "Direction",
      values = pathway_direction_palette(),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = "Direction",
      values = c("Control" = 21, "E-Stim" = 22),
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      name = "-log10 adjusted P",
      range = c(3, 8),
      trans = "sqrt"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(size = 4)),
      size = ggplot2::guide_legend(
        override.aes = list(fill = "#7A8793", color = "white")
      )
    ) +
    ggplot2::labs(
      x = "Method-specific enrichment score",
      y = NULL,
      caption = paste0(
        "Independent scales; dotted baselines: GSEA NES = 0, ",
        "ORA fold enrichment = 1.\n",
        "Color and shape = direction; point size = adjusted-P strength; ",
        "outline = selected pathway."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      plot.caption = ggplot2::element_text(
        color = "#52606D",
        hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

pathway_detail_ui <- function(row, genes) {
  if (nrow(row) == 0L) {
    return(NULL)
  }
  data <- prepare_pathway_plot_data(row)
  result_count <- paste(
    format(data$gene_count[[1L]], big.mark = ",", scientific = FALSE),
    if (identical(as.character(data$source[[1L]]), "GSEA")) {
      "genes in source set"
    } else {
      "overlapping genes in enrichment result"
    }
  )
  exported_label <- if (
    identical(as.character(data$source[[1L]]), "GSEA")
  ) {
    "Exported leading-edge genes"
  } else {
    "Exported overlapping genes"
  }
  displayed_genes <- paste(
    format(length(genes), big.mark = ","),
    "genes in exported list"
  )
  htmltools::div(
    class = "pathway-detail",
    htmltools::h2(data$label[[1L]]),
    htmltools::p(data$description[[1L]], class = "pathway-description"),
    htmltools::tags$dl(
      htmltools::tags$dt("Direction"),
      htmltools::tags$dd(as.character(data$direction[[1L]])),
      htmltools::tags$dt("Method"),
      htmltools::tags$dd(as.character(data$source[[1L]])),
      htmltools::tags$dt("Ontology"),
      htmltools::tags$dd("Gene Ontology (GO), Biological Process (BP)"),
      htmltools::tags$dt(data$score_label[[1L]]),
      htmltools::tags$dd(trimws(formatC(
        data$score[[1L]],
        digits = 4L,
        format = "g"
      ))),
      htmltools::tags$dt("Adjusted P"),
      htmltools::tags$dd(formatC(
        data$p_adjust[[1L]],
        digits = 4L,
        format = "g"
      )),
      htmltools::tags$dt(data$count_label[[1L]]),
      htmltools::tags$dd(result_count),
      htmltools::tags$dt(exported_label),
      htmltools::tags$dd(displayed_genes)
    )
  )
}

pathways_ui <- function(id) {
  ns <- shiny::NS(id)
  htmltools::div(
    class = "page-shell",
    htmltools::div(
      class = "page-heading",
      htmltools::div(
        htmltools::p("Complete enrichment analysis", class = "eyebrow"),
        htmltools::h1("Pathways"),
        htmltools::p(
          paste(
            "Top directional GSEA and ORA results from the primary condition",
            "analysis are shown on method-specific score scales. The complete",
            "enrichment results, including non-significant terms, are",
            "available below."
          ),
          class = "lede"
        )
      )
    ),
    bslib::layout_columns(
      col_widths = c(7, 5, 12),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Top pathway results"),
        htmltools::div(
          role = "region",
          `aria-label` = "Interactive pathway results",
          `aria-describedby` = ns("pathway_plot_note"),
          plotly::plotlyOutput(ns("pathway_plot"), height = "950px")
        ),
        htmltools::p(
          id = ns("pathway_plot_note"),
          class = "small text-body-secondary px-3 pb-3 mb-0",
          paste(
            "Each method shows the 10 terms with the lowest adjusted P",
            "values in each direction; a selected term outside that set is",
            "added to the plot.",
            "GSEA uses NES (neutral = 0); ORA uses fold enrichment",
            "(neutral = 1). Color and shape mark direction; point size",
            "marks adjusted-P strength, and an outline marks the selected",
            "pathway.",
            "Click any point to inspect it."
          )
        )
      ),
      bslib::card(
        bslib::card_header("Pathway details"),
        shiny::selectizeInput(
          ns("pathway"),
          "Pathway",
          choices = NULL,
          options = list(placeholder = "Search enrichment results")
        ),
        shiny::uiOutput(ns("pathway_detail")),
        bslib::layout_columns(
          col_widths = c(6, 6),
          shiny::actionButton(
            ns("explore_pathway"),
            "Explore gene set",
            class = "btn-primary w-100"
          ),
          shiny::actionButton(
            ns("add_pathway"),
            "Add to gene set",
            class = "btn-outline-primary w-100"
          )
        ),
        htmltools::div(
          class = "download-stack mt-3",
          shiny::downloadButton(
            ns("download_pathway"),
            "Pathway genes",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_pathway_png"),
            "Plot PNG",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_pathway_pdf"),
            "Plot PDF",
            class = "btn-outline-primary btn-sm"
          )
        )
      ),
      bslib::layout_columns(
        col_widths = c(9, 3),
        bslib::card(
          bslib::card_header("Enrichment results"),
          DT::DTOutput(ns("pathway_table"))
        ),
        bslib::card(
          bslib::card_header("Genes"),
          DT::DTOutput(ns("pathway_genes"))
        )
      )
    )
  )
}

pathways_server <- function(id, bundle, state, navigate_explore) {
  shiny::moduleServer(id, function(input, output, session) {
    source_id <- session$ns("pathway_source")
    plotted_pathways <- shiny::reactive({
      displayed_pathway_results(
        bundle$pathways,
        state$active_pathway()
      )
    })
    pathway_choices <- stats::setNames(
      as.character(bundle$pathways$pathway_id),
      paste0(
        as.character(bundle$pathways$label),
        " — ",
        as.character(bundle$pathways$source),
        ", ",
        as.character(bundle$pathways$direction)
      )
    )

    shiny::observe({
      shiny::updateSelectizeInput(
        session,
        "pathway",
        choices = pathway_choices,
        selected = state$active_pathway(),
        server = TRUE
      )
    })

    shiny::observeEvent(
      input$pathway,
      {
        if (
          !is.null(input$pathway) &&
            input$pathway %in% bundle$pathways$pathway_id
        ) {
          state$active_pathway(input$pathway)
        }
      },
      ignoreInit = TRUE
    )

    output$pathway_plot <- plotly::renderPlotly({
      make_pathway_plotly(
        plotted_pathways(),
        state$active_pathway(),
        source_id
      )
    })
    shiny::outputOptions(output, "pathway_plot", suspendWhenHidden = FALSE)

    shiny::observeEvent(
      plotly::event_data(
        "plotly_click",
        source = source_id,
        priority = "event"
      ),
      {
        event <- plotly::event_data("plotly_click", source = source_id)
        keys <- compact_character(event$key)
        if (length(keys) > 0L) state$active_pathway(keys[[1L]])
      },
      ignoreNULL = TRUE
    )

    active_row <- shiny::reactive({
      bundle$pathways[
        as.character(bundle$pathways$pathway_id) == state$active_pathway(),
        ,
        drop = FALSE
      ]
    })

    active_genes <- shiny::reactive({
      genes <- bundle$pathway_genes$gene[
        as.character(bundle$pathway_genes$pathway_id) == state$active_pathway()
      ]
      parse_gene_input(genes, bundle_gene_names(bundle))$genes
    })

    output$pathway_detail <- shiny::renderUI({
      pathway_detail_ui(active_row(), active_genes())
    })

    output$pathway_genes <- DT::renderDT({
      data <- data.frame(gene = active_genes(), stringsAsFactors = FALSE)
      DT::datatable(
        data,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        class = "compact stripe",
        options = list(pageLength = 25L, lengthChange = FALSE, dom = "ftip")
      )
    })

    output$pathway_table <- DT::renderDT({
      DT::datatable(
        bundle$pathways,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        class = "compact stripe",
        options = list(pageLength = 25L, lengthChange = FALSE, scrollX = TRUE)
      )
    })

    shiny::observeEvent(input$pathway_table_rows_selected, {
      row <- input$pathway_table_rows_selected
      if (length(row) > 0L && row[[1L]] <= nrow(bundle$pathways)) {
        state$active_pathway(as.character(bundle$pathways$pathway_id[[row[[
          1L
        ]]]]))
      }
    })

    shiny::observeEvent(input$pathway_genes_rows_selected, {
      row <- input$pathway_genes_rows_selected
      genes <- active_genes()
      if (length(row) > 0L && row[[1L]] <= length(genes)) {
        set_state_gene(state, bundle, genes[[row[[1L]]]])
      }
    })

    shiny::observeEvent(input$explore_pathway, {
      genes <- active_genes()
      if (length(genes) > 0L) {
        replace_state_gene_set(state, bundle, genes)
        navigate_explore()
      }
    })

    shiny::observeEvent(input$add_pathway, {
      genes <- active_genes()
      if (length(genes) > 0L) append_state_gene_set(state, bundle, genes)
    })

    output$download_pathway <- shiny::downloadHandler(
      filename = function() {
        row <- active_row()
        label <- if (nrow(row) > 0L) row$label[[1L]] else "pathway"
        paste0("espiviz-pathway-", safe_filename(label), ".txt")
      },
      content = function(file) writeLines(active_genes(), file, useBytes = TRUE)
    )

    output$download_pathway_png <- shiny::downloadHandler(
      filename = function() "espiviz-top-pathways.png",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_pathway_ggplot(
            plotted_pathways(),
            state$active_pathway()
          ),
          width = 9,
          height = 11,
          dpi = 320,
          bg = "white"
        )
      }
    )

    output$download_pathway_pdf <- shiny::downloadHandler(
      filename = function() "espiviz-top-pathways.pdf",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_pathway_ggplot(
            plotted_pathways(),
            state$active_pathway()
          ),
          device = grDevices::cairo_pdf,
          width = 9,
          height = 11,
          bg = "white"
        )
      }
    )
  })
}
