explore_ui <- function(id, bundle) {
  ns <- shiny::NS(id)
  preset_choices <- names(bundle$featured_gene_sets)
  cluster_choices <- sort(unique(as.character(bundle$cells$cluster)))

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      title = "Explore expression",
      shiny::selectizeInput(
        ns("active_gene"),
        "Gene",
        choices = NULL,
        options = list(placeholder = "Search 38,394 genes")
      ),
      shiny::selectInput(
        ns("color_by"),
        "Color UMAP by",
        choices = c(
          "Gene expression" = "expression",
          "Cluster" = "cluster",
          "Condition" = "condition"
        )
      ),
      shiny::actionButton(
        ns("clear_selection"),
        "Clear selection",
        class = "btn-outline-secondary w-100"
      ),
      htmltools::hr(),
      htmltools::h3("Gene set", class = "sidebar-section-title"),
      shiny::selectInput(
        ns("preset"),
        "Featured set",
        choices = c("Choose a set" = "", stats::setNames(preset_choices, preset_choices))
      ),
      shiny::actionButton(ns("use_preset"), "Use featured set", class = "btn-outline-primary w-100"),
      shiny::textAreaInput(
        ns("gene_text"),
        "Paste genes",
        rows = 4,
        placeholder = "One per line, or separated by spaces or commas"
      ),
      shiny::fileInput(
        ns("gene_file"),
        "Upload a gene list",
        accept = c(".txt", ".csv", ".tsv"),
        buttonLabel = "Choose file"
      ),
      bslib::layout_columns(
        col_widths = c(6, 6),
        shiny::actionButton(ns("add_gene_text"), "Add genes", class = "btn-primary w-100"),
        shiny::actionButton(ns("add_active_gene"), "Add current", class = "btn-outline-primary w-100")
      ),
      shiny::uiOutput(ns("gene_input_status")),
      shiny::uiOutput(ns("gene_set_status")),
      bslib::layout_columns(
        col_widths = c(6, 6),
        shiny::actionButton(ns("remove_genes"), "Remove chosen", class = "btn-outline-secondary w-100"),
        shiny::actionButton(ns("clear_gene_set"), "Clear set", class = "btn-outline-secondary w-100")
      )
    ),
    bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::card(
        class = "main-figure-card",
        full_screen = TRUE,
        bslib::card_header(
          htmltools::div(
            htmltools::span("UMAP", class = "card-kicker"),
            shiny::uiOutput(ns("umap_title"), inline = TRUE)
          )
        ),
        plotly::plotlyOutput(ns("umap"), height = "650px")
      ),
      bslib::card(
        class = "selection-card",
        bslib::card_header("Current selection"),
        shiny::uiOutput(ns("selection_overview")),
        htmltools::div(class = "selection-instruction", "Click a cell or use box or lasso from the plot toolbar."),
        htmltools::hr(),
        htmltools::h3("Downloads", class = "card-section-title"),
        shiny::downloadButton(ns("download_umap_png"), "UMAP PNG", class = "btn-outline-primary btn-sm"),
        shiny::downloadButton(ns("download_umap_pdf"), "UMAP PDF", class = "btn-outline-primary btn-sm"),
        shiny::downloadButton(ns("download_gene_set"), "Gene set", class = "btn-outline-primary btn-sm"),
        shiny::downloadButton(ns("download_summary"), "Selection summary", class = "btn-outline-primary btn-sm"),
        shiny::downloadButton(ns("download_expression"), "Cell expression", class = "btn-outline-primary btn-sm"),
        shiny::downloadButton(ns("download_metadata"), "Cell metadata", class = "btn-outline-primary btn-sm")
      ),
      bslib::card(
        class = "span-12",
        bslib::navset_card_tab(
          id = ns("explore_results"),
          bslib::nav_panel(
            "Comparison",
            htmltools::div(
              class = "result-toolbar",
              shiny::uiOutput(ns("comparison_heading")),
              shiny::numericInput(ns("gene_page"), "Page", value = 1L, min = 1L, step = 1L, width = "110px")
            ),
            shiny::plotOutput(ns("gene_comparison_plot"), height = "auto"),
            DT::DTOutput(ns("comparison_table"))
          ),
          bslib::nav_panel("By condition", DT::DTOutput(ns("condition_table"))),
          bslib::nav_panel("By cluster", DT::DTOutput(ns("cluster_table"))),
          bslib::nav_panel(
            "Gene set",
            htmltools::p("The complete set stays active across pages and downloads.", class = "supporting-copy"),
            DT::DTOutput(ns("gene_set_table"))
          ),
          bslib::nav_panel(
            "Top markers",
            htmltools::div(
              class = "result-toolbar",
              shiny::selectInput(ns("marker_cluster"), "Cluster", choices = cluster_choices, width = "150px"),
              shiny::actionButton(ns("open_marker"), "Explore chosen marker", class = "btn-primary")
            ),
            DT::DTOutput(ns("marker_table"))
          )
        )
      )
    )
  )
}

summary_datatable <- function(data, page_length = 25L) {
  numeric_columns <- intersect(
    c(
      "selected_mean", "selected_median", "selected_detected_pct",
      "remaining_mean", "remaining_median", "remaining_detected_pct",
      "mean_difference", "mean_ratio", "detection_pp_difference",
      "detection_ratio", "mean_expression", "median_expression", "detected_pct"
    ),
    names(data)
  )
  table <- DT::datatable(
    data,
    rownames = FALSE,
    escape = TRUE,
    class = "compact stripe",
    options = list(
      pageLength = page_length,
      lengthChange = FALSE,
      scrollX = TRUE,
      autoWidth = TRUE,
      language = list(emptyTable = "No rows to display")
    )
  )
  if (length(numeric_columns) > 0L) {
    table <- DT::formatRound(table, columns = numeric_columns, digits = 3L)
  }
  table
}

explore_server <- function(id, bundle, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    universe <- bundle_gene_names(bundle)
    source_id <- ns("umap_source")
    input_message <- shiny::reactiveVal(NULL)
    analysis_genes <- state_analysis_genes(state)

    shiny::observe({
      shiny::updateSelectizeInput(
        session,
        "active_gene",
        choices = universe,
        selected = state$active_gene(),
        server = TRUE
      )
    })

    shiny::observeEvent(input$active_gene, {
      if (!is.null(input$active_gene) && nzchar(input$active_gene)) {
        set_state_gene(state, bundle, input$active_gene)
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(state$active_gene(), {
      if (!identical(input$active_gene, state$active_gene())) {
        shiny::updateSelectizeInput(
          session,
          "active_gene",
          selected = state$active_gene(),
          server = TRUE
        )
      }
    }, ignoreInit = TRUE)

    add_text_genes <- function(text) {
      parsed <- append_state_gene_set(state, bundle, text)
      if (length(parsed$missing) > 0L) {
        input_message(paste("Not found:", paste(parsed$missing, collapse = ", ")))
      } else {
        input_message(NULL)
      }
      parsed
    }

    shiny::observeEvent(input$add_gene_text, {
      text <- input$gene_text %||% ""
      if (!is.null(input$gene_file) && file.exists(input$gene_file$datapath)) {
        uploaded <- paste(read_gene_upload(input$gene_file), collapse = "\n")
        text <- paste(text, uploaded, sep = "\n")
      }
      add_text_genes(text)
    })

    shiny::observeEvent(input$add_active_gene, {
      append_state_gene_set(state, bundle, state$active_gene())
      input_message(NULL)
    })

    shiny::observeEvent(input$use_preset, {
      preset <- input$preset %||% ""
      if (nzchar(preset) && preset %in% names(bundle$featured_gene_sets)) {
        replace_state_gene_set(state, bundle, bundle$featured_gene_sets[[preset]])
        input_message(NULL)
      }
    })

    shiny::observeEvent(input$clear_gene_set, {
      state$gene_set(character())
      input_message(NULL)
    })

    shiny::observeEvent(input$remove_genes, {
      rows <- input$gene_set_table_rows_selected %||% integer()
      current <- state$gene_set()
      if (length(rows) > 0L && length(current) > 0L) {
        state$gene_set(current[-rows])
      }
    })

    output$gene_input_status <- shiny::renderUI({
      message <- input_message()
      if (is.null(message)) return(NULL)
      htmltools::div(message, class = "input-status", role = "status")
    })

    output$gene_set_status <- shiny::renderUI({
      count <- length(state$gene_set())
      htmltools::div(
        if (count == 0L) "No gene set active" else paste(format(count, big.mark = ","), "genes active"),
        class = "gene-set-status"
      )
    })

    shiny::observeEvent(input$clear_selection, {
      state$selected_cells(character())
    })

    output$umap_title <- shiny::renderUI({
      label <- switch(
        input$color_by %||% "expression",
        expression = state$active_gene(),
        cluster = "Final cluster",
        condition = "Condition"
      )
      htmltools::span(label, class = "figure-subtitle")
    })

    output$umap <- plotly::renderPlotly({
      make_umap_plotly(
        bundle = bundle,
        color_by = input$color_by %||% "expression",
        gene = state$active_gene(),
        selected_cell_ids = state$selected_cells(),
        source = source_id
      )
    })

    shiny::observeEvent(plotly::event_data(
      "plotly_selected",
      source = source_id,
      priority = "event"
    ), {
      event <- plotly::event_data("plotly_selected", source = source_id)
      if (!is.null(event)) {
        keys <- unique(compact_character(event$key))
        state$selected_cells(intersect(bundle$cells$cell_id, keys))
      }
    }, ignoreNULL = TRUE)

    shiny::observeEvent(plotly::event_data(
      "plotly_click",
      source = source_id,
      priority = "event"
    ), {
      event <- plotly::event_data("plotly_click", source = source_id)
      keys <- unique(compact_character(event$key))
      if (length(keys) > 0L) state$selected_cells(keys[[1L]])
    }, ignoreNULL = TRUE)

    shiny::observeEvent(plotly::event_data(
      "plotly_deselect",
      source = source_id,
      priority = "event"
    ), {
      state$selected_cells(character())
    }, ignoreNULL = TRUE)

    selection_summary <- shiny::reactive({
      visible_genes <- paginate_genes(
        analysis_genes(),
        input$gene_page %||% 1L,
        25L
      )$genes
      summarize_selection(bundle, state$selected_cells(), visible_genes)
    })

    output$selection_overview <- shiny::renderUI({
      summary <- selection_summary()
      selected_n <- length(summary$selected_cell_ids)
      remaining_n <- length(summary$remaining_cell_ids)
      explicit <- length(state$selected_cells()) > 0L
      htmltools::div(
        class = "selection-overview",
        htmltools::div(
          class = "selection-count",
          htmltools::strong(format(selected_n, big.mark = ",")),
          htmltools::span(if (explicit) "selected cells" else "cells")
        ),
        htmltools::div(
          class = "selection-detail",
          if (explicit) {
            paste(format(remaining_n, big.mark = ","), "remaining")
          } else {
            "No region selected"
          }
        )
      )
    })

    page_info <- shiny::reactive({
      paginate_genes(analysis_genes(), input$gene_page %||% 1L, 25L)
    })

    shiny::observe({
      page <- page_info()
      shiny::updateNumericInput(
        session,
        "gene_page",
        value = page$page,
        min = 1L,
        max = page$pages
      )
    })

    output$comparison_heading <- shiny::renderUI({
      page <- page_info()
      htmltools::div(
        htmltools::h3("Selected cells and remaining cells", class = "result-title"),
        htmltools::p(
          if (page$total == 0L) {
            "No genes"
          } else {
            paste0(
              "Genes ", page$start, "–", page$end, " of ",
              format(page$total, big.mark = ",")
            )
          },
          class = "supporting-copy"
        )
      )
    })

    page_comparison <- shiny::reactive({
      genes <- page_info()$genes
      selection_summary()$comparison[
        match(genes, selection_summary()$comparison$gene),
        ,
        drop = FALSE
      ]
    })

    output$gene_comparison_plot <- shiny::renderPlot(
      {
        plot <- make_gene_comparison_plot(page_comparison())
        if (is.null(plot)) return(invisible(NULL))
        plot
      },
      height = function() max(260L, 28L * max(1L, length(page_info()$genes)))
    )

    output$comparison_table <- DT::renderDT({
      summary_datatable(page_comparison(), page_length = 25L)
    })

    output$condition_table <- DT::renderDT({
      genes <- page_info()$genes
      data <- selection_summary()$selected_by_condition
      data <- data[data$gene %in% genes, , drop = FALSE]
      summary_datatable(data, page_length = 25L)
    })

    output$cluster_table <- DT::renderDT({
      genes <- page_info()$genes
      data <- selection_summary()$selected_by_cluster
      data <- data[data$gene %in% genes, , drop = FALSE]
      summary_datatable(data, page_length = 25L)
    })

    output$gene_set_table <- DT::renderDT({
      genes <- state$gene_set()
      data <- data.frame(gene = genes, stringsAsFactors = FALSE)
      DT::datatable(
        data,
        rownames = FALSE,
        selection = "multiple",
        class = "compact stripe",
        options = list(pageLength = 25L, lengthChange = FALSE, dom = "tip")
      )
    })

    marker_data <- shiny::reactive({
      cluster <- as.character(input$marker_cluster %||% "")
      data <- bundle$markers[as.character(bundle$markers$cluster) == cluster, , drop = FALSE]
      data[order(data$rank), , drop = FALSE]
    })

    output$marker_table <- DT::renderDT({
      DT::datatable(
        marker_data(),
        rownames = FALSE,
        selection = "single",
        class = "compact stripe",
        options = list(pageLength = 15L, lengthChange = FALSE, scrollX = TRUE)
      )
    })

    shiny::observeEvent(input$open_marker, {
      row <- input$marker_table_rows_selected %||% integer()
      data <- marker_data()
      if (length(row) > 0L && row[[1L]] <= nrow(data)) {
        set_state_gene(state, bundle, data$gene[[row[[1L]]]])
      }
    })

    output$download_umap_png <- shiny::downloadHandler(
      filename = function() paste0("espiviz-umap-", safe_filename(state$active_gene()), ".png"),
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_umap_ggplot(
            bundle,
            input$color_by %||% "expression",
            state$active_gene(),
            state$selected_cells()
          ),
          width = 8.5,
          height = 7,
          dpi = 320,
          bg = "white"
        )
      }
    )

    output$download_umap_pdf <- shiny::downloadHandler(
      filename = function() paste0("espiviz-umap-", safe_filename(state$active_gene()), ".pdf"),
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_umap_ggplot(
            bundle,
            input$color_by %||% "expression",
            state$active_gene(),
            state$selected_cells()
          ),
          device = grDevices::cairo_pdf,
          width = 8.5,
          height = 7,
          bg = "white"
        )
      }
    )

    output$download_gene_set <- shiny::downloadHandler(
      filename = function() "espiviz-gene-set.txt",
      content = function(file) writeLines(analysis_genes(), file, useBytes = TRUE)
    )

    output$download_summary <- shiny::downloadHandler(
      filename = function() "espiviz-selection-summary.csv",
      content = function(file) {
        utils::write.csv(
          selection_summary_export(bundle, state$selected_cells(), analysis_genes()),
          file,
          row.names = FALSE,
          na = ""
        )
      }
    )

    output$download_expression <- shiny::downloadHandler(
      filename = function() "espiviz-selected-cell-expression.rds",
      content = function(file) {
        write_selection_expression_export(
          bundle,
          state$selected_cells(),
          analysis_genes(),
          file
        )
      }
    )

    output$download_metadata <- shiny::downloadHandler(
      filename = function() "espiviz-selected-cell-metadata.csv",
      content = function(file) {
        selected <- state$selected_cells()
        if (length(selected) == 0L) selected <- bundle$cells$cell_id
        data <- bundle$cells[bundle$cells$cell_id %in% selected, , drop = FALSE]
        utils::write.csv(data, file, row.names = FALSE, na = "")
      }
    )
  })
}
