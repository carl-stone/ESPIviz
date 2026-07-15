de_column <- function(data, candidates, required = TRUE) {
  match <- candidates[candidates %in% names(data)]
  if (length(match) > 0L) {
    return(match[[1L]])
  }
  if (isTRUE(required)) {
    stop(
      "The differential-expression table is missing a required column.",
      call. = FALSE
    )
  }
  NULL
}

de_significance_levels <- function() {
  c(
    "Higher in Control",
    "Not significant",
    "Higher with E-Stim",
    "Adjusted P unavailable"
  )
}

de_significance_palette <- function() {
  c(
    "Higher in Control" = "#3569a8",
    "Not significant" = "#7f8b93",
    "Higher with E-Stim" = "#b64b35",
    "Adjusted P unavailable" = "#6b4f8f"
  )
}

de_significance_symbols <- function() {
  c(
    "Higher in Control" = "square",
    "Not significant" = "circle",
    "Higher with E-Stim" = "triangle-up",
    "Adjusted P unavailable" = "x"
  )
}

prepare_de_plot_data <- function(primary_de, fdr_threshold = 0.05) {
  if (
    length(fdr_threshold) != 1L ||
      is.na(fdr_threshold) ||
      !is.finite(fdr_threshold) ||
      fdr_threshold < 0 ||
      fdr_threshold > 1
  ) {
    stop("The FDR threshold must be one number from zero to one.", call. = FALSE)
  }
  gene_column <- de_column(primary_de, c("gene"))
  base_mean_column <- de_column(primary_de, c("baseMean", "base_mean"))
  fold_change_column <- de_column(
    primary_de,
    c("log2FoldChange", "log2_fold_change", "logFC")
  )
  probability_column <- de_column(
    primary_de,
    c("padj", "p_adjust", "FDR")
  )
  data <- data.frame(
    gene = as.character(primary_de[[gene_column]]),
    base_mean = as.numeric(primary_de[[base_mean_column]]),
    log2_fold_change = as.numeric(primary_de[[fold_change_column]]),
    probability = as.numeric(primary_de[[probability_column]]),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(casefold_key(data$gene))) {
    stop("Differential-expression gene keys must be unique.", call. = FALSE)
  }
  data$minus_log10_probability <- rep(NA_real_, nrow(data))
  probability_available <- !is.na(data$probability)
  data$minus_log10_probability[probability_available] <- -log10(pmax(
    data$probability[probability_available],
    .Machine$double.xmin
  ))

  significance <- rep("Not significant", nrow(data))
  significance[!probability_available] <- "Adjusted P unavailable"
  significance[
    probability_available &
      data$probability <= fdr_threshold &
      data$log2_fold_change < 0
  ] <- "Higher in Control"
  significance[
    probability_available &
      data$probability <= fdr_threshold &
      data$log2_fold_change > 0
  ] <- "Higher with E-Stim"
  data$significance <- factor(
    significance,
    levels = de_significance_levels()
  )
  data[
    c(
      "gene",
      "log2_fold_change",
      "probability",
      "minus_log10_probability",
      "base_mean",
      "significance"
    )
  ]
}

format_de_value <- function(value, digits = 4L) {
  formatted <- trimws(formatC(
    value,
    digits = digits,
    format = "fg",
    big.mark = ","
  ))
  formatted[is.na(value)] <- "Unavailable"
  formatted
}

format_de_probability <- function(value) {
  formatted <- trimws(formatC(value, digits = 4L, format = "g"))
  formatted[is.na(value)] <- "Unavailable"
  formatted
}

de_direction <- function(log2_fold_change) {
  if (length(log2_fold_change) == 0L || is.na(log2_fold_change[[1L]])) {
    return("Unavailable")
  }
  if (log2_fold_change[[1L]] > 0) {
    return("Higher with E-Stim")
  }
  if (log2_fold_change[[1L]] < 0) {
    return("Higher in Control")
  }
  "No directional change"
}

de_hover_text <- function(data) {
  paste0(
    htmltools::htmlEscape(data$gene),
    "<br>baseMean: ",
    format_de_value(data$base_mean),
    "<br>Shrunken log2 fold change: ",
    format_de_value(data$log2_fold_change),
    "<br>Adjusted P: ",
    format_de_probability(data$probability),
    "<br>FDR status: ",
    as.character(data$significance)
  )
}

de_plot_coordinates <- function(data, view) {
  view <- match.arg(view, c("ma", "volcano"))
  if (identical(view, "ma")) {
    data$plot_x <- pmax(data$base_mean, 0) + 1
    data$plot_y <- data$log2_fold_change
  } else {
    data$plot_x <- data$log2_fold_change
    data$plot_y <- data$minus_log10_probability
    data$plot_y[is.na(data$plot_y)] <- 0
  }
  data
}

make_de_plotly <- function(
  primary_de,
  active_gene,
  source,
  view = c("ma", "volcano")
) {
  view <- match.arg(view)
  data <- prepare_de_plot_data(primary_de)
  data <- de_plot_coordinates(data, view)
  data$hover <- de_hover_text(data)
  active <- casefold_key(data$gene) == casefold_key(active_gene)
  highlighted <- data[active, , drop = FALSE]
  background <- data[!active, , drop = FALSE]
  palette <- de_significance_palette()
  symbols <- de_significance_symbols()

  plot <- plotly::plot_ly(source = source)
  for (status in de_significance_levels()) {
    trace_data <- background[
      as.character(background$significance) == status,
      ,
      drop = FALSE
    ]
    if (nrow(trace_data) == 0L) next
    plot <- plotly::add_trace(
      plot,
      data = trace_data,
      x = ~plot_x,
      y = ~plot_y,
      key = ~gene,
      text = ~hover,
      hoverinfo = "text",
      type = "scattergl",
      mode = "markers",
      marker = list(
        size = if (identical(status, "Adjusted P unavailable")) 6.5 else 4.5,
        color = unname(palette[[status]]),
        symbol = unname(symbols[[status]]),
        opacity = if (identical(status, "Not significant")) 0.42 else 0.72,
        line = list(width = 0)
      ),
      name = status,
      legendgroup = status,
      inherit = FALSE,
      showlegend = TRUE
    )
  }
  if (nrow(highlighted) > 0L) {
    highlighted_status <- as.character(highlighted$significance[[1L]])
    plot <- plotly::add_trace(
      plot,
      data = highlighted,
      x = ~plot_x,
      y = ~plot_y,
      key = ~gene,
      text = ~hover,
      hoverinfo = "text",
      type = "scatter",
      mode = "markers",
      marker = list(
        size = 11,
        color = unname(palette[[highlighted_status]]),
        symbol = unname(symbols[[highlighted_status]]),
        opacity = 1,
        line = list(color = "#17232b", width = 2)
      ),
      inherit = FALSE,
      showlegend = FALSE
    )
  }

  if (identical(view, "ma")) {
    positive_x <- data$plot_x[is.finite(data$plot_x) & data$plot_x > 0]
    decade_ticks <- 10^seq(
      floor(log10(min(positive_x))),
      ceiling(log10(max(positive_x)))
    )
    xaxis <- list(
      title = "Mean normalized count (baseMean + 1; log scale)",
      type = "log",
      tickmode = "array",
      tickvals = decade_ticks,
      ticktext = format(
        decade_ticks,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      )
    )
    yaxis <- list(
      title = "Shrunken log2 fold change (E-Stim / Control)",
      zeroline = FALSE
    )
    shapes <- list(list(
      type = "line",
      x0 = 0,
      x1 = 1,
      xref = "paper",
      y0 = 0,
      y1 = 0,
      line = list(color = "#68757d", width = 1, dash = "dash")
    ))
    annotations <- list()
  } else {
    fdr_reference <- -log10(0.05)
    xaxis <- list(
      title = "Shrunken log2 fold change (E-Stim / Control)",
      zeroline = TRUE,
      zerolinecolor = "#b8c1c8"
    )
    yaxis <- list(title = "−log10 adjusted P", rangemode = "tozero")
    shapes <- list(list(
      type = "line",
      x0 = 0,
      x1 = 1,
      xref = "paper",
      y0 = fdr_reference,
      y1 = fdr_reference,
      line = list(color = "#68757d", width = 1, dash = "dash")
    ))
    annotations <- list(list(
      x = 1,
      xref = "paper",
      xanchor = "right",
      y = fdr_reference,
      yanchor = "bottom",
      text = "FDR = 0.05",
      showarrow = FALSE,
      font = list(color = "#526b7b", size = 11)
    ))
  }
  plot <- plotly::layout(
    plot,
    xaxis = xaxis,
    yaxis = yaxis,
    shapes = shapes,
    annotations = annotations,
    legend = list(
      title = list(text = "FDR status"),
      orientation = "h",
      x = 0,
      y = -0.2,
      xanchor = "left",
      yanchor = "top"
    ),
    margin = list(l = 76, r = 20, t = 18, b = 108),
    hovermode = "closest"
  )
  plot <- plotly::config(
    plot,
    displaylogo = FALSE,
    toImageButtonOptions = list(
      filename = paste0("espiviz-primary-de-", view)
    )
  )
  plotly::event_register(plot, "plotly_click")
}

make_de_ggplot <- function(
  primary_de,
  active_gene,
  view = c("ma", "volcano")
) {
  view <- match.arg(view)
  data <- prepare_de_plot_data(primary_de)
  data <- de_plot_coordinates(data, view)
  data$active <- casefold_key(data$gene) == casefold_key(active_gene)
  palette <- de_significance_palette()
  shapes <- c(
    "Higher in Control" = 15,
    "Not significant" = 16,
    "Higher with E-Stim" = 17,
    "Adjusted P unavailable" = 4
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = plot_x,
      y = plot_y,
      color = significance,
      shape = significance
    )
  )
  if (identical(view, "ma")) {
    plot <- plot +
      ggplot2::geom_hline(
        yintercept = 0,
        color = "#68757d",
        linewidth = 0.35,
        linetype = "dashed"
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        x = "Mean normalized count (baseMean + 1; log scale)",
        y = "Shrunken log2 fold change (E-Stim / Control)"
      )
  } else {
    plot <- plot +
      ggplot2::geom_hline(
        yintercept = -log10(0.05),
        color = "#68757d",
        linewidth = 0.35,
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "text",
        x = Inf,
        y = -log10(0.05),
        label = "FDR = 0.05",
        hjust = 1.05,
        vjust = -0.5,
        color = "#526b7b",
        size = 3.2
      ) +
      ggplot2::labs(
        x = "Shrunken log2 fold change (E-Stim / Control)",
        y = "−log10 adjusted P"
      )
  }
  plot +
    ggplot2::geom_point(alpha = 0.62, size = 1.15, stroke = 0.75) +
    ggplot2::geom_point(
      data = data[data$active, , drop = FALSE],
      alpha = 1,
      size = 2.8,
      stroke = 1
    ) +
    ggplot2::geom_point(
      data = data[data$active, , drop = FALSE],
      ggplot2::aes(x = plot_x, y = plot_y),
      inherit.aes = FALSE,
      color = "#17232b",
      shape = 1,
      size = 4,
      stroke = 0.85
    ) +
    ggplot2::scale_color_manual(
      values = palette,
      breaks = de_significance_levels(),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = shapes,
      breaks = de_significance_levels(),
      drop = FALSE
    ) +
    ggplot2::labs(
      color = "FDR status",
      shape = "FDR status",
      caption = paste(
        "Positive fold changes indicate higher expression with E-Stim.",
        "Crosses mark genes whose adjusted P value is unavailable."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.caption = ggplot2::element_text(
        color = "#526b7b",
        hjust = 0,
        margin = ggplot2::margin(t = 10)
      )
    )
}

current_de_gene_ui <- function(primary_de, active_gene) {
  data <- prepare_de_plot_data(primary_de)
  current <- data[
    casefold_key(data$gene) == casefold_key(active_gene),
    ,
    drop = FALSE
  ]
  if (nrow(current) == 0L) {
    return(htmltools::div(
      class = "current-gene-display",
      htmltools::strong(active_gene),
      htmltools::span("No result is available in the primary condition model.")
    ))
  }
  current <- current[1L, , drop = FALSE]
  htmltools::div(
    class = "current-gene-display de-gene-detail",
    htmltools::strong(current$gene[[1L]]),
    htmltools::span("Selected across the app"),
    htmltools::tags$dl(
      htmltools::tags$dt("baseMean"),
      htmltools::tags$dd(format_de_value(current$base_mean[[1L]])),
      htmltools::tags$dt("Shrunken log2FC"),
      htmltools::tags$dd(format_de_value(current$log2_fold_change[[1L]])),
      htmltools::tags$dt("Adjusted P"),
      htmltools::tags$dd(format_de_probability(current$probability[[1L]])),
      htmltools::tags$dt("Direction"),
      htmltools::tags$dd(de_direction(current$log2_fold_change[[1L]]))
    )
  )
}

differential_expression_ui <- function(id) {
  ns <- shiny::NS(id)
  htmltools::div(
    class = "page-shell",
    htmltools::div(
      class = "page-heading",
      htmltools::div(
        htmltools::p("Condition model", class = "eyebrow"),
        htmltools::h1("Differential expression"),
        htmltools::p(
          "Complete gene-level results from the primary Mouse × Condition pseudobulk analysis.",
          class = "lede"
        )
      ),
      htmltools::div(
        class = "heading-actions",
        shiny::downloadButton(
          ns("download_de"),
          "Download complete table",
          class = "btn-primary"
        )
      )
    ),
    bslib::layout_columns(
      col_widths = c(8, 4, 12),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header(
          htmltools::div(
            class = "de-plot-header",
            htmltools::span("Gene-level effect overview"),
            shiny::radioButtons(
              ns("de_view"),
              "Plot view",
              choices = c("MA plot" = "ma", "Volcano plot" = "volcano"),
              selected = "ma",
              inline = TRUE
            )
          )
        ),
        htmltools::p(
          paste(
            "The MA plot puts effect size in the context of expression strength;",
            "the volcano plot emphasizes statistical evidence.",
            "Positive fold changes indicate higher expression with E-Stim.",
            "Crosses mark genes whose adjusted P value is unavailable."
          ),
          id = ns("de_plot_help"),
          class = "supporting-copy de-plot-copy"
        ),
        htmltools::div(
          role = "region",
          `aria-label` = "Interactive differential-expression plot",
          `aria-describedby` = ns("de_plot_help"),
          plotly::plotlyOutput(ns("de_plot"), height = "560px")
        )
      ),
      bslib::card(
        bslib::card_header("Current gene"),
        shiny::uiOutput(ns("current_gene")),
        shiny::actionButton(
          ns("explore_current"),
          "Open in Explore",
          class = "btn-primary w-100"
        ),
        htmltools::hr(),
        htmltools::p(
          "Choose a point or table row to change the current gene.",
          class = "supporting-copy"
        ),
        htmltools::p("Plot downloads", class = "sidebar-section-title"),
        htmltools::div(
          class = "download-stack",
          shiny::downloadButton(
            ns("download_ma_png"),
            "MA PNG",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_ma_pdf"),
            "MA PDF",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_volcano_png"),
            "Volcano PNG",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_volcano_pdf"),
            "Volcano PDF",
            class = "btn-outline-primary btn-sm"
          )
        )
      ),
      bslib::card(
        class = "span-12",
        bslib::card_header("Complete primary result"),
        DT::DTOutput(ns("de_table"))
      )
    )
  )
}

differential_expression_server <- function(
  id,
  bundle,
  state,
  navigate_explore
) {
  shiny::moduleServer(id, function(input, output, session) {
    source_id <- session$ns("de_source")
    de_view <- shiny::reactive({
      view <- compact_character(input$de_view)
      if (length(view) > 0L && view[[1L]] %in% c("ma", "volcano")) {
        view[[1L]]
      } else {
        "ma"
      }
    })

    output$de_plot <- plotly::renderPlotly({
      make_de_plotly(
        bundle$primary_de,
        state$active_gene(),
        source_id,
        view = de_view()
      )
    })
    shiny::outputOptions(output, "de_plot", suspendWhenHidden = FALSE)

    output$current_gene <- shiny::renderUI({
      current_de_gene_ui(bundle$primary_de, state$active_gene())
    })

    shiny::observeEvent(
      plotly::event_data(
        "plotly_click",
        source = source_id,
        priority = "event"
      ),
      {
        event <- plotly::event_data("plotly_click", source = source_id)
        genes <- compact_character(event$key)
        if (length(genes) > 0L) set_state_gene(state, bundle, genes[[1L]])
      },
      ignoreNULL = TRUE
    )

    output$de_table <- DT::renderDT(
      {
        DT::datatable(
          bundle$primary_de,
          rownames = FALSE,
          filter = "top",
          selection = "single",
          class = "compact stripe",
          extensions = "Scroller",
          options = list(
            deferRender = TRUE,
            scrollX = TRUE,
            scrollY = 610,
            scroller = TRUE,
            pageLength = 50L,
            lengthChange = FALSE,
            search = list(caseInsensitive = TRUE)
          )
        )
      },
      server = TRUE
    )

    shiny::observeEvent(input$de_table_rows_selected, {
      row <- input$de_table_rows_selected
      if (length(row) > 0L && row[[1L]] <= nrow(bundle$primary_de)) {
        set_state_gene(state, bundle, bundle$primary_de$gene[[row[[1L]]]])
      }
    })

    shiny::observeEvent(input$explore_current, navigate_explore())

    output$download_de <- shiny::downloadHandler(
      filename = function() {
        "espiviz-primary-condition-differential-expression.csv"
      },
      content = function(file) {
        utils::write.csv(bundle$primary_de, file, row.names = FALSE, na = "")
      }
    )

    output$download_ma_png <- shiny::downloadHandler(
      filename = function() "espiviz-primary-de-ma.png",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_de_ggplot(
            bundle$primary_de,
            state$active_gene(),
            view = "ma"
          ),
          width = 8.5,
          height = 7,
          dpi = 320,
          bg = "white"
        )
      }
    )

    output$download_ma_pdf <- shiny::downloadHandler(
      filename = function() "espiviz-primary-de-ma.pdf",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_de_ggplot(
            bundle$primary_de,
            state$active_gene(),
            view = "ma"
          ),
          device = grDevices::cairo_pdf,
          width = 8.5,
          height = 7,
          bg = "white"
        )
      }
    )

    output$download_volcano_png <- shiny::downloadHandler(
      filename = function() "espiviz-primary-de-volcano.png",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_de_ggplot(
            bundle$primary_de,
            state$active_gene(),
            view = "volcano"
          ),
          width = 8.5,
          height = 7,
          dpi = 320,
          bg = "white"
        )
      }
    )

    output$download_volcano_pdf <- shiny::downloadHandler(
      filename = function() "espiviz-primary-de-volcano.pdf",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_de_ggplot(
            bundle$primary_de,
            state$active_gene(),
            view = "volcano"
          ),
          device = grDevices::cairo_pdf,
          width = 8.5,
          height = 7,
          bg = "white"
        )
      }
    )
  })
}
