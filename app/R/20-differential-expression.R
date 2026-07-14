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

prepare_de_plot_data <- function(primary_de) {
  gene_column <- de_column(primary_de, c("gene"))
  fold_change_column <- de_column(
    primary_de,
    c("log2FoldChange", "log2_fold_change", "logFC")
  )
  probability_column <- de_column(
    primary_de,
    c("padj", "p_adjust", "FDR", "pvalue", "p_value")
  )
  data <- data.frame(
    gene = as.character(primary_de[[gene_column]]),
    log2_fold_change = as.numeric(primary_de[[fold_change_column]]),
    probability = as.numeric(primary_de[[probability_column]]),
    stringsAsFactors = FALSE
  )
  data$minus_log10_probability <- ifelse(
    is.na(data$probability),
    0,
    -log10(pmax(data$probability, .Machine$double.xmin))
  )
  data
}

make_de_plotly <- function(primary_de, active_gene, source) {
  data <- prepare_de_plot_data(primary_de)
  hover <- paste0(
    data$gene,
    "<br>log2 fold change: ",
    formatC(data$log2_fold_change, digits = 4L, format = "fg"),
    "<br>Adjusted P: ",
    formatC(data$probability, digits = 4L, format = "g")
  )
  plot <- plotly::plot_ly(
    data = data,
    x = ~log2_fold_change,
    y = ~minus_log10_probability,
    key = ~gene,
    text = hover,
    hoverinfo = "text",
    type = "scattergl",
    mode = "markers",
    marker = list(
      size = 4.5,
      color = "#526b7b",
      opacity = 0.48,
      line = list(width = 0)
    ),
    source = source,
    showlegend = FALSE
  )
  highlighted <- data[
    casefold_key(data$gene) == casefold_key(active_gene),
    ,
    drop = FALSE
  ]
  if (nrow(highlighted) > 0L) {
    plot <- plotly::add_trace(
      plot,
      data = highlighted,
      x = ~log2_fold_change,
      y = ~minus_log10_probability,
      key = ~gene,
      text = paste0(
        highlighted$gene,
        "<br>log2 fold change: ",
        formatC(highlighted$log2_fold_change, digits = 4L, format = "fg"),
        "<br>Adjusted P: ",
        formatC(highlighted$probability, digits = 4L, format = "g")
      ),
      hoverinfo = "text",
      type = "scatter",
      mode = "markers",
      marker = list(
        size = 10,
        color = "#b52865",
        line = list(color = "white", width = 1.2)
      ),
      inherit = FALSE,
      showlegend = FALSE
    )
  }
  plot <- plotly::layout(
    plot,
    xaxis = list(
      title = "log2 fold change",
      zeroline = TRUE,
      zerolinecolor = "#b8c1c8"
    ),
    yaxis = list(title = "−log10 adjusted P", rangemode = "tozero"),
    margin = list(l = 70, r = 20, t = 18, b = 58),
    hovermode = "closest"
  )
  plot <- plotly::config(
    plot,
    displaylogo = FALSE,
    toImageButtonOptions = list(filename = "espiviz-primary-de-volcano")
  )
  plotly::event_register(plot, "plotly_click")
}

make_de_ggplot <- function(primary_de, active_gene) {
  data <- prepare_de_plot_data(primary_de)
  data$active <- casefold_key(data$gene) == casefold_key(active_gene)
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = log2_fold_change, y = minus_log10_probability)
  ) +
    ggplot2::geom_point(color = "#526b7b", alpha = 0.45, size = 1) +
    ggplot2::geom_point(
      data = data[data$active, , drop = FALSE],
      color = "#b52865",
      size = 2.8
    ) +
    ggplot2::labs(x = "log2 fold change", y = "−log10 adjusted P") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
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
        bslib::card_header("Volcano plot"),
        plotly::plotlyOutput(ns("volcano"), height = "590px")
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
        htmltools::div(
          class = "download-stack",
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

    output$volcano <- plotly::renderPlotly({
      make_de_plotly(bundle$primary_de, state$active_gene(), source_id)
    })
    shiny::outputOptions(output, "volcano", suspendWhenHidden = FALSE)

    output$current_gene <- shiny::renderUI({
      htmltools::div(
        class = "current-gene-display",
        htmltools::strong(state$active_gene()),
        htmltools::span("Selected across the app")
      )
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

    output$download_volcano_png <- shiny::downloadHandler(
      filename = function() "espiviz-primary-de-volcano.png",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_de_ggplot(bundle$primary_de, state$active_gene()),
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
          plot = make_de_ggplot(bundle$primary_de, state$active_gene()),
          device = grDevices::cairo_pdf,
          width = 8.5,
          height = 7,
          bg = "white"
        )
      }
    )
  })
}
