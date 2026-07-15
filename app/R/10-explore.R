explore_ui <- function(id, bundle) {
  ns <- shiny::NS(id)
  preset_choices <- names(bundle$featured_gene_sets)
  cluster_choices <- sort(unique(as.character(bundle$cells$cluster)))
  gene_pair_group_choices <- c(
    "All cells" = "all",
    stats::setNames(cluster_choices, paste("Cluster", cluster_choices))
  )

  bslib::layout_sidebar(
    fillable = FALSE,
    fill = FALSE,
    sidebar = bslib::sidebar(
      width = 330,
      title = "Explore expression",
      shiny::selectizeInput(
        ns("active_gene"),
        "Gene",
        choices = NULL,
        options = list(placeholder = "Search 38,394 genes")
      ),
      shiny::selectizeInput(
        ns("secondary_gene"),
        "Second plot gene (optional)",
        choices = NULL,
        selected = NULL,
        options = list(
          placeholder = "Add a gene for blend and scatter",
          allowEmptyOption = TRUE,
          create = FALSE
        )
      ),
      htmltools::p(
        paste(
          "Current gene is always the first plot gene.",
          "The optional second gene affects plots only; it does not change",
          "the current gene or gene set."
        ),
        class = "sidebar-help"
      ),
      shiny::selectInput(
        ns("color_by"),
        "Color UMAP by",
        choices = c(
          "Log normalized expression" = "expression",
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
        choices = c(
          "Choose a set" = "",
          stats::setNames(preset_choices, preset_choices)
        )
      ),
      shiny::actionButton(
        ns("use_preset"),
        "Use featured set",
        class = "btn-outline-primary w-100"
      ),
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
        shiny::actionButton(
          ns("add_gene_text"),
          "Add genes",
          class = "btn-primary w-100"
        ),
        shiny::actionButton(
          ns("add_active_gene"),
          "Add current",
          class = "btn-outline-primary w-100"
        )
      ),
      shiny::uiOutput(ns("gene_input_status")),
      shiny::uiOutput(ns("gene_set_status")),
      bslib::layout_columns(
        col_widths = c(6, 6),
        shiny::actionButton(
          ns("remove_genes"),
          "Remove chosen",
          class = "btn-outline-secondary w-100"
        ),
        shiny::actionButton(
          ns("clear_gene_set"),
          "Clear set",
          class = "btn-outline-secondary w-100"
        )
      )
    ),
    bslib::layout_columns(
      col_widths = c(8, 4, 12),
      fillable = FALSE,
      fill = FALSE,
      class = "explore-grid",
      bslib::card(
        class = "main-figure-card",
        full_screen = TRUE,
        bslib::card_header(
          htmltools::div(
            htmltools::span("UMAP", class = "card-kicker"),
            shiny::uiOutput(ns("umap_title"), inline = TRUE)
          )
        ),
        htmltools::div(
          role = "region",
          `aria-label` = "Interactive UMAP",
          `aria-describedby` = ns("umap_note"),
          plotly::plotlyOutput(ns("umap"), height = "650px")
        ),
        htmltools::p(
          id = ns("umap_note"),
          class = "figure-note",
          paste(
            "UMAP is a qualitative two-dimensional projection.",
            "Use local neighborhoods for exploration; do not interpret",
            "global distance, cluster area, or point density quantitatively."
          )
        )
      ),
      htmltools::div(
        class = "umap-side-rail",
        bslib::card(
          class = "selection-card selection-card-compact",
          bslib::card_header("Current selection"),
          shiny::uiOutput(ns("selection_overview")),
          htmltools::div(
            class = "selection-instruction",
            "Click a cell or use box or lasso on the UMAP."
          ),
          htmltools::h3(
            "Select cells by cluster",
            class = "card-section-title"
          ),
          shiny::selectInput(
            ns("select_cluster"),
            label = NULL,
            choices = c(
              "Choose cluster" = "",
              stats::setNames(
                cluster_choices,
                paste("Cluster", cluster_choices)
              )
            )
          ),
          shiny::actionButton(
            ns("select_cluster_cells"),
            "Select cluster",
            class = "btn-primary w-100"
          )
        ),
        bslib::card(
          class = "selection-summary-card",
          bslib::card_header("Selected-cell summary"),
          shiny::uiOutput(ns("selection_snapshot")),
          shiny::uiOutput(ns("umap_legend"))
        ),
        bslib::card(
          class = "selection-download-card",
          bslib::card_header("Downloads"),
          htmltools::div(
            class = "download-stack compact-downloads",
            shiny::downloadButton(
              ns("download_umap_png"),
              "UMAP PNG",
              class = "btn-outline-primary btn-sm"
            ),
            shiny::downloadButton(
              ns("download_umap_pdf"),
              "UMAP PDF",
              class = "btn-outline-primary btn-sm"
            ),
            shiny::downloadButton(
              ns("download_gene_set"),
              "Gene set",
              class = "btn-outline-primary btn-sm"
            ),
            shiny::downloadButton(
              ns("download_summary"),
              "Selection summary",
              class = "btn-outline-primary btn-sm"
            ),
            shiny::downloadButton(
              ns("download_expression"),
              "Cell expression",
              class = "btn-outline-primary btn-sm"
            ),
            shiny::downloadButton(
              ns("download_metadata"),
              "Cell metadata",
              class = "btn-outline-primary btn-sm"
            )
          )
        )
      ),
      bslib::card(
        class = "results-card",
        htmltools::div(
          class = "result-toolbar result-page-toolbar",
          htmltools::div(
            htmltools::h3("Expression summaries", class = "result-title"),
            shiny::uiOutput(ns("gene_page_status"), inline = TRUE)
          ),
          shiny::numericInput(
            ns("gene_page"),
            "Page",
            value = 1L,
            min = 1L,
            step = 1L,
            width = "110px"
          )
        ),
        bslib::navset_card_tab(
          id = ns("explore_results"),
          selected = "By cluster",
          bslib::nav_panel(
            "Comparison",
            shiny::uiOutput(ns("comparison_heading")),
            htmltools::div(
              class = "comparison-plot-shell",
              shiny::uiOutput(ns("gene_comparison_plot_ui"))
            ),
            DT::DTOutput(ns("comparison_table"))
          ),
          bslib::nav_panel(
            "Cell-level",
            htmltools::h3(
              "Log normalized expression distribution by final cluster",
              class = "result-title"
            ),
            htmltools::p(
              paste(
                "Violins show the current gene and optional second plot gene",
                "across final clusters. Log normalized expression is centered",
                "and can be negative;",
                "outlined points mark explicitly selected cells."
              ),
              class = "supporting-copy"
            ),
            shiny::uiOutput(ns("violin_plot_ui")),
            htmltools::hr(class = "result-divider"),
            htmltools::h3(
              "Two-gene log normalized expression",
              class = "result-title"
            ),
            htmltools::div(
              class = "result-toolbar gene-pair-controls",
              shiny::selectInput(
                ns("gene_pair_display"),
                "Display",
                choices = c(
                  "Scatter plot" = "scatter",
                  "Density plot" = "density"
                ),
                selected = "scatter",
                width = "220px"
              ),
              shiny::conditionalPanel(
                condition = "input.gene_pair_display === 'scatter'",
                ns = ns,
                shiny::selectInput(
                  ns("gene_pair_loess_group"),
                  "Loess trend",
                  choices = gene_pair_group_choices,
                  selected = "all",
                  width = "220px"
                )
              ),
              shiny::conditionalPanel(
                condition = "input.gene_pair_display === 'density'",
                ns = ns,
                shiny::selectInput(
                  ns("gene_pair_density_group"),
                  "Density cells",
                  choices = gene_pair_group_choices,
                  selected = "all",
                  width = "220px"
                )
              )
            ),
            shiny::conditionalPanel(
              condition = "input.gene_pair_display === 'scatter'",
              ns = ns,
              htmltools::p(
                paste(
                  "Each point is a cell detected by raw count for at least one gene;",
                  "color identifies its final cluster.",
                  "Larger outlined diamonds are explicitly selected cells."
                ),
                class = "supporting-copy"
              )
            ),
            shiny::conditionalPanel(
              condition = "input.gene_pair_display === 'density'",
              ns = ns,
              htmltools::p(
                paste(
                  "Filled contours show cell density for cells detected by raw",
                  "count for at least one gene."
                ),
                class = "supporting-copy"
              )
            ),
            shiny::uiOutput(ns("gene_pair_plot_ui"))
          ),
          bslib::nav_panel(
            "By sample",
            htmltools::h3(
              "Expression by biological sample",
              class = "result-title"
            ),
            htmltools::p(
              paste(
                "Samples are the biological replicates and are ordered by",
                "condition. Mean log normalized expression is shown by color; detected-cell",
                "percentage is shown by dot size."
              ),
              class = "supporting-copy"
            ),
            shiny::uiOutput(ns("sample_summary_plot_ui")),
            DT::DTOutput(ns("sample_table")),
            htmltools::hr(class = "result-divider"),
            htmltools::h3(
              "Observed cluster composition",
              class = "result-title"
            ),
            htmltools::p(
              paste(
                "Each tile reports recovered cells and their percentage",
                "within one sample. This is descriptive, not a",
                "differential-abundance analysis."
              ),
              class = "supporting-copy"
            ),
            shiny::plotOutput(ns("composition_plot"), height = "520px"),
            DT::DTOutput(ns("composition_table"))
          ),
          bslib::nav_panel(
            "Pooled condition",
            htmltools::h3(
              "Expression by condition — pooled cells",
              class = "result-title"
            ),
            htmltools::p(
              paste(
                "This descriptive view pools cells within each condition.",
                "Use the sample view above to inspect replicate consistency.",
                "Mean log normalized expression is shown by color and detected-cell",
                "percentage by dot size."
              ),
              class = "supporting-copy"
            ),
            shiny::uiOutput(ns("condition_summary_plot_ui")),
            DT::DTOutput(ns("condition_table"))
          ),
          bslib::nav_panel(
            "By cluster",
            htmltools::h3(
              "Expression by final cluster",
              class = "result-title"
            ),
            htmltools::p(
              "Mean log normalized expression is shown by color; detected-cell percentage is shown by dot size.",
              class = "supporting-copy"
            ),
            shiny::uiOutput(ns("cluster_summary_plot_ui")),
            DT::DTOutput(ns("cluster_table"))
          ),
          bslib::nav_panel(
            "Cluster markers",
            htmltools::div(
              class = "result-toolbar",
              shiny::selectInput(
                ns("marker_cluster"),
                "Marker genes from cluster",
                choices = cluster_choices,
                width = "220px"
              )
            ),
            htmltools::p(
              paste(
                "The top eight marker genes for the chosen cluster are",
                "shown across every cluster. Color is scaled within each",
                "gene, so compare patterns across clusters rather than",
                "absolute color between genes."
              ),
              class = "supporting-copy"
            ),
            shiny::uiOutput(ns("marker_overview_plot_ui")),
            htmltools::p(
              "Select a marker row below to make that gene current across the app.",
              class = "supporting-copy"
            ),
            DT::DTOutput(ns("marker_table"))
          ),
          bslib::nav_panel(
            "Gene set",
            htmltools::p(
              "The complete set stays active across pages and downloads.",
              class = "supporting-copy"
            ),
            DT::DTOutput(ns("gene_set_table"))
          )
        )
      )
    )
  )
}

summary_datatable <- function(data, page_length = 25L) {
  labels <- c(
    gene = "Gene",
    condition = "Condition",
    cluster = "Cluster",
    sample = "Sample",
    cell_count = "Cells",
    sample_total = "Sample cells",
    selected_n = "Selected cells",
    remaining_n = "Remaining cells",
    selected_mean = "Selected mean log normalized expression",
    selected_median = "Selected median log normalized expression",
    selected_detected_n = "Selected detected cells",
    selected_detected_pct = "Selected detected (%)",
    remaining_mean = "Other cells mean log normalized expression",
    remaining_median = "Remaining median log normalized expression",
    remaining_detected_n = "Remaining detected cells",
    remaining_detected_pct = "Other cells detected (%)",
    mean_difference = "Mean difference",
    detection_pp_difference = "Detection difference (pp)",
    detection_ratio = "Detection ratio",
    mean_expression = "Mean log normalized expression",
    median_expression = "Median log normalized expression",
    detected_n = "Detected cells",
    detected_pct = "Detected (%)",
    cell_pct = "Within sample (%)"
  )
  numeric_names <- c(
    "selected_mean",
    "selected_median",
    "selected_detected_pct",
    "remaining_mean",
    "remaining_median",
    "remaining_detected_pct",
    "mean_difference",
    "detection_pp_difference",
    "detection_ratio",
    "mean_expression",
    "median_expression",
    "detected_pct",
    "cell_pct"
  )
  original_names <- names(data)
  display_names <- original_names
  matched <- match(original_names, names(labels))
  display_names[!is.na(matched)] <- unname(labels[matched[!is.na(matched)]])
  names(data) <- display_names
  numeric_columns <- display_names[original_names %in% numeric_names]
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
    pair_source_id <- ns("gene_pair_source")
    input_message <- shiny::reactiveVal(NULL)
    analysis_genes <- state_analysis_genes(state)
    composition_data <- prepare_cluster_composition(bundle)
    show_selection_comparison <- function() {
      bslib::nav_select(
        "explore_results",
        selected = "Comparison",
        session = session
      )
    }

    shiny::observe({
      shiny::updateSelectizeInput(
        session,
        "active_gene",
        choices = universe,
        selected = state$active_gene(),
        server = TRUE
      )
    })

    shiny::observeEvent(
      input$active_gene,
      {
        if (!is.null(input$active_gene) && nzchar(input$active_gene)) {
          set_state_gene(state, bundle, input$active_gene)
        }
      },
      ignoreInit = TRUE
    )

    shiny::observe({
      shiny::updateSelectizeInput(
        session,
        "secondary_gene",
        choices = c(
          "No second gene" = "",
          stats::setNames(universe, universe)
        ),
        selected = "",
        server = TRUE
      )
    })

    shiny::observeEvent(
      state$active_gene(),
      {
        current <- compact_character(input$secondary_gene %||% "")
        if (
          length(current) > 0L &&
            casefold_key(current[[1L]]) == casefold_key(state$active_gene())
        ) {
          shiny::updateSelectizeInput(
            session,
            "secondary_gene",
            choices = c(
              "No second gene" = "",
              stats::setNames(universe, universe)
            ),
            selected = "",
            server = TRUE
          )
        }
      },
      ignoreInit = TRUE
    )

    plot_genes <- shiny::reactive({
      second <- parse_gene_input(
        input$secondary_gene %||% "",
        universe,
        max_genes = 1L
      )$genes
      normalize_plot_genes(
        bundle,
        c(state$active_gene(), second),
        state$active_gene()
      )
    })

    expression_summary_genes <- shiny::reactive({
      genes <- state$gene_set()
      if (length(genes) == 0L) plot_genes() else genes
    })

    plot_gene_data <- shiny::reactive({
      prepare_plot_gene_data(
        bundle,
        plot_genes(),
        state$selected_cells()
      )
    })

    add_text_genes <- function(text) {
      parsed <- append_state_gene_set(state, bundle, text)
      if (length(parsed$missing) > 0L) {
        input_message(paste(
          "Not found:",
          paste(parsed$missing, collapse = ", ")
        ))
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
        replace_state_gene_set(
          state,
          bundle,
          bundle$featured_gene_sets[[preset]]
        )
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
      if (is.null(message)) {
        return(NULL)
      }
      htmltools::div(message, class = "input-status", role = "status")
    })

    output$gene_set_status <- shiny::renderUI({
      count <- length(state$gene_set())
      htmltools::div(
        if (count == 0L) {
          "No gene set active"
        } else {
          paste(format(count, big.mark = ","), "genes active")
        },
        class = "gene-set-status"
      )
    })

    shiny::observeEvent(input$clear_selection, {
      state$selected_cells(character())
    })

    shiny::observeEvent(
      input$select_cluster_cells,
      {
        selected <- cell_ids_for_cluster(bundle, input$select_cluster %||% "")
        if (length(selected) == 0L) {
          return()
        }
        state$selected_cells(selected)
        shiny::updateSelectInput(session, "color_by", selected = "cluster")
        show_selection_comparison()
      },
      ignoreInit = TRUE
    )

    output$umap_title <- shiny::renderUI({
      label <- switch(
        input$color_by %||% "expression",
        expression = if (length(plot_genes()) == 2L) {
          paste(plot_genes(), collapse = " + ") |>
            paste("blend")
        } else {
          plot_genes()[[1L]]
        },
        cluster = "Final cluster",
        condition = "Condition"
      )
      htmltools::span(label, class = "figure-subtitle")
    })

    output$umap_legend <- shiny::renderUI({
      if (
        !identical(input$color_by %||% "expression", "expression") ||
          length(plot_genes()) != 2L
      ) {
        return(NULL)
      }
      blend_legend_ui(plot_genes())
    })

    output$umap <- plotly::renderPlotly({
      make_umap_plotly(
        bundle = bundle,
        color_by = input$color_by %||% "expression",
        gene = plot_genes(),
        selected_cell_ids = state$selected_cells(),
        source = source_id,
        gene_data = plot_gene_data()
      )
    })

    shiny::observeEvent(
      plotly::event_data(
        "plotly_selected",
        source = source_id,
        priority = "event"
      ),
      {
        event <- plotly::event_data("plotly_selected", source = source_id)
        if (!is.null(event)) {
          keys <- unique(compact_character(event$key))
          state$selected_cells(intersect(bundle$cells$cell_id, keys))
          show_selection_comparison()
        }
      },
      ignoreNULL = TRUE
    )

    shiny::observeEvent(
      plotly::event_data(
        "plotly_click",
        source = source_id,
        priority = "event"
      ),
      {
        event <- plotly::event_data("plotly_click", source = source_id)
        keys <- unique(compact_character(event$key))
        if (length(keys) > 0L) {
          state$selected_cells(keys[[1L]])
          show_selection_comparison()
        }
      },
      ignoreNULL = TRUE
    )

    selection_summary <- shiny::reactive({
      visible_genes <- paginate_genes(
        expression_summary_genes(),
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

    output$selection_snapshot <- shiny::renderUI({
      selection_snapshot_ui(prepare_selection_snapshot(plot_gene_data()))
    })

    page_info <- shiny::reactive({
      paginate_genes(
        expression_summary_genes(),
        input$gene_page %||% 1L,
        25L
      )
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
      remaining_n <- length(selection_summary()$remaining_cell_ids)
      explicit <- length(state$selected_cells()) > 0L
      heading <- if (!explicit) {
        "Expression summary — all cells"
      } else if (remaining_n == 0L) {
        "Expression summary — all cells selected"
      } else {
        "Selected cells and remaining cells"
      }
      htmltools::h3(heading, class = "result-title")
    })

    output$gene_page_status <- shiny::renderUI({
      page <- page_info()
      htmltools::span(
        if (page$total == 0L) {
          "No genes"
        } else {
          paste0(
            "Genes ",
            page$start,
            "–",
            page$end,
            " of ",
            format(page$total, big.mark = ",")
          )
        },
        class = "supporting-copy"
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

    output$gene_comparison_plot_ui <- shiny::renderUI({
      shiny::plotOutput(
        ns("gene_comparison_plot"),
        height = paste0(comparison_plot_height(length(page_info()$genes)), "px")
      )
    })

    output$gene_comparison_plot <- shiny::renderPlot(
      {
        plot <- make_gene_comparison_plot(page_comparison())
        if (is.null(plot)) {
          return(invisible(NULL))
        }
        plot
      },
      height = function() comparison_plot_height(length(page_info()$genes)),
      alt = function() {
        paste(
          "Dot plot comparing log normalized expression for",
          length(page_info()$genes),
          "genes across all cells, selected cells, and remaining cells.",
          "Color shows mean log normalized expression; dot size shows the percentage of",
          "cells with detected expression."
        )
      }
    )

    output$comparison_table <- DT::renderDT({
      concise <- prepare_comparison_table(
        page_comparison(),
        explicit_selection = length(state$selected_cells()) > 0L
      )
      summary_datatable(concise, page_length = 25L)
    })

    violin_height <- shiny::reactive({
      if (length(plot_genes()) == 2L) 650L else 430L
    })

    output$violin_plot_ui <- shiny::renderUI({
      htmltools::div(
        role = "region",
        `aria-label` = paste(
          "Log normalized expression distribution violin plot for",
          paste(plot_genes(), collapse = " and "),
          "by final cluster"
        ),
        shiny::plotOutput(
          ns("violin_plot"),
          height = paste0(violin_height(), "px")
        )
      )
    })

    output$violin_plot <- shiny::renderPlot(
      make_violin_plot(plot_gene_data(), bundle),
      height = function() violin_height(),
      alt = function() {
        paste(
          "Violin plot of log normalized expression for",
          paste(plot_genes(), collapse = " and "),
          "across final clusters. Log normalized expression is centered and",
          "can be negative.",
          "Outlined points mark explicitly selected cells; detection is based",
          "on raw counts and is reported in the selected-cell summary."
        )
      }
    )

    output$gene_pair_plot_ui <- shiny::renderUI({
      if (length(plot_genes()) != 2L) {
        return(htmltools::div(
          class = "plot-empty",
          role = "status",
          "Choose a second plot gene in the sidebar to show the two-gene plot."
        ))
      }
      gene_data <- plot_gene_data()
      scope <- gene_pair_scope(gene_data)
      if (scope$included_n == 0L) {
        return(htmltools::div(
          class = "plot-empty",
          role = "status",
          "Neither gene is detected by raw count in any cell."
        ))
      }
      display <- input$gene_pair_display %||% "scatter"
      display_label <- if (identical(display, "density")) {
        "density plot"
      } else {
        "scatter"
      }
      display_scope <- if (identical(display, "density")) {
        gene_pair_scope(
          gene_data,
          group = input$gene_pair_density_group %||% "all"
        )
      } else {
        scope
      }
      htmltools::div(
        role = "region",
        `aria-label` = paste(
          "Two-gene log normalized expression",
          display_label,
          "comparing",
          plot_genes()[[1L]],
          "and",
          plot_genes()[[2L]],
          "for cells detected for at least one gene"
        ),
        gene_pair_scope_ui(display_scope),
        plotly::plotlyOutput(ns("gene_pair_plot"), height = "540px")
      )
    })

    output$gene_pair_plot <- plotly::renderPlotly({
      shiny::req(length(plot_genes()) == 2L)
      make_gene_pair_plotly(
        plot_gene_data(),
        bundle,
        source = pair_source_id,
        display = input$gene_pair_display %||% "scatter",
        loess_group = input$gene_pair_loess_group %||% "all",
        density_group = input$gene_pair_density_group %||% "all"
      )
    })

    condition_page_summary <- shiny::reactive({
      genes <- page_info()$genes
      data <- selection_summary()$selected_by_condition
      data[data$gene %in% genes, , drop = FALSE]
    })

    sample_page_summary <- shiny::reactive({
      genes <- page_info()$genes
      data <- selection_summary()$selected_by_sample
      data[data$gene %in% genes, , drop = FALSE]
    })

    cluster_page_summary <- shiny::reactive({
      genes <- page_info()$genes
      data <- selection_summary()$selected_by_cluster
      data[data$gene %in% genes, , drop = FALSE]
    })

    output$condition_summary_plot_ui <- shiny::renderUI({
      shiny::plotOutput(
        ns("condition_summary_plot"),
        height = paste0(comparison_plot_height(length(page_info()$genes)), "px")
      )
    })

    output$condition_summary_plot <- shiny::renderPlot(
      make_group_summary_plot(
        condition_page_summary(),
        group_column = "condition",
        group_label = "Condition"
      ),
      height = function() comparison_plot_height(length(page_info()$genes)),
      alt = function() {
        paste(
          "Dot plot of log normalized expression for",
          length(page_info()$genes),
          "genes across pooled conditions. Color shows mean log normalized expression;",
          "dot size shows the percentage of selected cells with detected",
          "expression."
        )
      }
    )

    output$sample_summary_plot_ui <- shiny::renderUI({
      shiny::plotOutput(
        ns("sample_summary_plot"),
        height = paste0(comparison_plot_height(length(page_info()$genes)), "px")
      )
    })

    output$sample_summary_plot <- shiny::renderPlot(
      make_group_summary_plot(
        sample_page_summary(),
        group_column = "sample",
        group_label = "Biological sample",
        rotate_x = TRUE
      ),
      height = function() comparison_plot_height(length(page_info()$genes)),
      alt = function() {
        paste(
          "Dot plot of log normalized expression for",
          length(page_info()$genes),
          "genes across",
          length(unique(sample_page_summary()$sample)),
          "biological samples ordered by condition. Color shows mean log",
          "normalized expression; dot size shows the percentage of selected cells",
          "with detected expression."
        )
      }
    )

    output$composition_plot <- shiny::renderPlot(
      {
        make_cluster_composition_plot(composition_data, bundle)
      },
      alt = function() {
        paste(
          "Heatmap of recovered cell counts and within-sample percentages",
          "for",
          length(unique(composition_data$sample)),
          "biological samples across",
          length(unique(composition_data$cluster)),
          "final clusters, faceted by condition. This is descriptive and",
          "does not report differential abundance."
        )
      }
    )

    output$cluster_summary_plot_ui <- shiny::renderUI({
      shiny::plotOutput(
        ns("cluster_summary_plot"),
        height = paste0(comparison_plot_height(length(page_info()$genes)), "px")
      )
    })

    output$cluster_summary_plot <- shiny::renderPlot(
      make_group_summary_plot(
        cluster_page_summary(),
        group_column = "cluster",
        group_label = "Final cluster"
      ),
      height = function() comparison_plot_height(length(page_info()$genes)),
      alt = function() {
        paste(
          "Dot plot of log normalized expression for",
          length(page_info()$genes),
          "genes across final clusters. Color shows mean log normalized expression; dot",
          "size shows the percentage of selected cells with detected",
          "expression."
        )
      }
    )

    output$condition_table <- DT::renderDT({
      summary_datatable(condition_page_summary(), page_length = 25L)
    })

    output$sample_table <- DT::renderDT({
      summary_datatable(sample_page_summary(), page_length = 25L)
    })

    output$composition_table <- DT::renderDT({
      data <- composition_data[
        ,
        c(
          "condition",
          "sample",
          "cluster",
          "cell_count",
          "sample_total",
          "cell_pct"
        ),
        drop = FALSE
      ]
      summary_datatable(data, page_length = 48L)
    })

    output$cluster_table <- DT::renderDT({
      summary_datatable(cluster_page_summary(), page_length = 25L)
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
      data <- bundle$markers[
        as.character(bundle$markers$cluster) == cluster,
        ,
        drop = FALSE
      ]
      data[order(data$rank), , drop = FALSE]
    })

    marker_overview <- shiny::reactive({
      prepare_marker_overview(
        bundle,
        marker_cluster = input$marker_cluster %||% "",
        top_n = 8L
      )
    })

    output$marker_overview_plot_ui <- shiny::renderUI({
      gene_count <- length(unique(marker_overview()$gene))
      shiny::plotOutput(
        ns("marker_overview_plot"),
        height = paste0(marker_overview_height(gene_count), "px")
      )
    })

    output$marker_overview_plot <- shiny::renderPlot(
      {
        plot <- make_marker_overview_plot(marker_overview())
        if (is.null(plot)) {
          return(invisible(NULL))
        }
        plot
      },
      height = function() {
        marker_overview_height(length(unique(marker_overview()$gene)))
      },
      alt = function() {
        genes <- unique(as.character(marker_overview()$gene))
        paste0(
          "Dot plot of the top marker genes for final cluster ",
          input$marker_cluster %||% "",
          " across all final clusters. Genes shown: ",
          paste(genes, collapse = ", "),
          ". Color is relative mean log normalized expression scaled within each gene; ",
          "dot size is the percentage of cells with detected expression."
        )
      }
    )

    output$marker_table <- DT::renderDT({
      DT::datatable(
        marker_data(),
        rownames = FALSE,
        selection = "single",
        class = "compact stripe",
        options = list(pageLength = 15L, lengthChange = FALSE, scrollX = TRUE)
      )
    })

    shiny::observeEvent(input$marker_table_rows_selected, {
      row <- input$marker_table_rows_selected %||% integer()
      data <- marker_data()
      if (length(row) > 0L && row[[1L]] <= nrow(data)) {
        set_state_gene(state, bundle, data$gene[[row[[1L]]]])
      }
    }, ignoreInit = TRUE)

    output$download_umap_png <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "espiviz-umap-",
          safe_filename(paste(plot_genes(), collapse = "-")),
          ".png"
        )
      },
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_umap_ggplot(
            bundle,
            input$color_by %||% "expression",
            plot_genes(),
            state$selected_cells(),
            gene_data = plot_gene_data()
          ),
          width = 8.5,
          height = 7,
          dpi = 320,
          bg = "white"
        )
      }
    )

    output$download_umap_pdf <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "espiviz-umap-",
          safe_filename(paste(plot_genes(), collapse = "-")),
          ".pdf"
        )
      },
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = make_umap_ggplot(
            bundle,
            input$color_by %||% "expression",
            plot_genes(),
            state$selected_cells(),
            gene_data = plot_gene_data()
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
      content = function(file) {
        writeLines(analysis_genes(), file, useBytes = TRUE)
      }
    )

    output$download_summary <- shiny::downloadHandler(
      filename = function() "espiviz-selection-summary.csv",
      content = function(file) {
        utils::write.csv(
          selection_summary_export(
            bundle,
            state$selected_cells(),
            analysis_genes()
          ),
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
        if (length(selected) == 0L) {
          selected <- bundle$cells$cell_id
        }
        data <- bundle$cells[bundle$cells$cell_id %in% selected, , drop = FALSE]
        utils::write.csv(data, file, row.names = FALSE, na = "")
      }
    )
  })
}
