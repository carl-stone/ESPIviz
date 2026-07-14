prepare_pathway_plot_data <- function(pathways) {
  data <- pathways
  data$evidence <- -log10(pmax(as.numeric(data$p_adjust), .Machine$double.xmin))
  data$label <- as.character(data$label)
  data$pathway_id <- as.character(data$pathway_id)
  data
}

make_pathway_plotly <- function(pathways, active_pathway, source) {
  data <- prepare_pathway_plot_data(pathways)
  data <- data[order(data$evidence), , drop = FALSE]
  data$active <- data$pathway_id == active_pathway
  plotly::plot_ly(
    data = data,
    x = ~evidence,
    y = ~ factor(label, levels = label),
    key = ~pathway_id,
    text = paste0(
      data$label,
      "<br>Adjusted P: ",
      formatC(data$p_adjust, digits = 4L, format = "g"),
      "<br>Genes: ",
      data$gene_count
    ),
    hoverinfo = "text",
    type = "bar",
    orientation = "h",
    marker = list(
      color = ifelse(data$active, "#b52865", "#587687"),
      line = list(width = 0)
    ),
    source = source,
    showlegend = FALSE
  ) |>
    plotly::layout(
      xaxis = list(title = "−log10 adjusted P", rangemode = "tozero"),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 30, r = 20, t = 20, b = 55)
    ) |>
    plotly::config(displaylogo = FALSE) |>
    plotly::event_register("plotly_click")
}

make_pathway_ggplot <- function(pathways, active_pathway) {
  data <- prepare_pathway_plot_data(pathways)
  data <- data[order(data$evidence), , drop = FALSE]
  data$label <- factor(data$label, levels = data$label)
  data$active <- data$pathway_id == active_pathway
  ggplot2::ggplot(data, ggplot2::aes(x = evidence, y = label, fill = active)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#b52865", `FALSE` = "#587687"),
      guide = "none"
    ) +
    ggplot2::labs(x = "−log10 adjusted P", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

pathway_detail_ui <- function(row, genes) {
  if (nrow(row) == 0L) {
    return(NULL)
  }
  htmltools::div(
    class = "pathway-detail",
    htmltools::h2(row$label[[1L]]),
    htmltools::p(row$description[[1L]], class = "pathway-description"),
    htmltools::tags$dl(
      htmltools::tags$dt("Method"),
      htmltools::tags$dd(as.character(row$source[[1L]])),
      htmltools::tags$dt("Adjusted P"),
      htmltools::tags$dd(formatC(
        row$p_adjust[[1L]],
        digits = 4L,
        format = "g"
      )),
      htmltools::tags$dt("Genes"),
      htmltools::tags$dd(format(length(genes), big.mark = ","))
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
        htmltools::p("Manuscript-aligned terms", class = "eyebrow"),
        htmltools::h1("Pathways"),
        htmltools::p(
          "A focused set of pathway results from the primary condition analysis.",
          class = "lede"
        )
      )
    ),
    bslib::layout_columns(
      col_widths = c(7, 5, 12),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Featured pathway results"),
        plotly::plotlyOutput(ns("pathway_plot"), height = "540px")
      ),
      bslib::card(
        bslib::card_header("Pathway details"),
        shiny::selectInput(ns("pathway"), "Pathway", choices = NULL),
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
      bslib::card(
        class = "span-12",
        bslib::navset_card_tab(
          bslib::nav_panel("Genes", DT::DTOutput(ns("pathway_genes"))),
          bslib::nav_panel(
            "All featured terms",
            DT::DTOutput(ns("pathway_table"))
          )
        )
      )
    )
  )
}

pathways_server <- function(id, bundle, state, navigate_explore) {
  shiny::moduleServer(id, function(input, output, session) {
    source_id <- session$ns("pathway_source")
    pathway_choices <- stats::setNames(
      as.character(bundle$pathways$pathway_id),
      as.character(bundle$pathways$label)
    )

    shiny::observe({
      shiny::updateSelectInput(
        session,
        "pathway",
        choices = pathway_choices,
        selected = state$active_pathway()
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
      make_pathway_plotly(bundle$pathways, state$active_pathway(), source_id)
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
      filename = function() "espiviz-featured-pathways.png",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_pathway_ggplot(bundle$pathways, state$active_pathway()),
          width = 9,
          height = 6.5,
          dpi = 320,
          bg = "white"
        )
      }
    )

    output$download_pathway_pdf <- shiny::downloadHandler(
      filename = function() "espiviz-featured-pathways.pdf",
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_pathway_ggplot(bundle$pathways, state$active_pathway()),
          device = grDevices::cairo_pdf,
          width = 9,
          height = 6.5,
          bg = "white"
        )
      }
    )
  })
}
