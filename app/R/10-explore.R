explore_ui <- function(id, bundle) {
  ns <- shiny::NS(id)
  preset_choices <- names(bundle$featured_gene_sets)
  cluster_choices <- sort(unique(as.character(bundle$cells$cluster)))
  gene_pair_group_choices <- c(
    "All cells" = "all",
    stats::setNames(cluster_choices, paste("Cluster", cluster_choices))
  )
  summary_plot_type_input <- function(id) {
    shiny::selectInput(
      ns(id),
      "Plot type",
      choices = c(
        "Automatic" = "auto",
        "Violin plot" = "violin",
        "Dot plot" = "dot"
      ),
      selected = "auto",
      width = "190px"
    )
  }
  summary_plot_note <- function(id) {
    htmltools::tags$details(
      class = "plot-help",
      htmltools::tags$summary("How to read this figure"),
      htmltools::p(
        "Log normalized expression is centered and may be negative. Detection means raw count > 0. Dot color shows mean log normalized expression; dot area shows detection percentage. Zero detection has zero area; a cross means no cells. Shared color limits span all analysis genes and summary groups, including other pages."
      ),
      htmltools::p(
        "Violin width shows density, scaled within each group. Boxes show median and interquartile range. Groups below 10 cells show individual points and a median; n counts cells. Automatic uses dot plots above six genes."
      )
    )
  }
  tool_menu <- function(label, ..., class = "") {
    htmltools::tags$details(
      class = paste("tool-menu", class),
      htmltools::tags$summary(label),
      htmltools::div(class = "tool-content", ...)
    )
  }

  htmltools::div(
    class = "page-shell explore-shell",
    htmltools::div(
      class = "page-heading",
      htmltools::h1("Explore"),
      htmltools::p(
        "Mouse retina · p27CKO (p27 inactivation), with and without electrical stimulation (E-Stim)",
        class = "page-context"
      )
    ),
    htmltools::div(
      class = "explore-toolbar",
      shiny::selectizeInput(
        ns("active_gene"),
        "Current gene",
        choices = NULL,
        options = list(placeholder = "Search genes")
      ),
      shiny::selectizeInput(
        ns("secondary_gene"),
        "Second gene (optional)",
        choices = NULL,
        selected = NULL,
        options = list(
          placeholder = "Compare another gene",
          allowEmptyOption = TRUE,
          create = FALSE
        )
      ),
      shiny::selectInput(
        ns("color_by"),
        "Color UMAP by",
        choices = c(
          "Log normalized expression" = "expression",
          "Detection" = "detection",
          "Cluster" = "cluster",
          "Condition" = "condition"
        )
      ),
      shiny::uiOutput(ns("selection_overview"))
    ),
    htmltools::div(
      class = "explore-actions",
      tool_menu(
        "Gene set",
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
          "Replace with featured set",
          class = "btn-outline-primary"
        ),
        htmltools::hr(),
        shiny::textAreaInput(
          ns("gene_text"),
          "Paste genes",
          rows = 3,
          placeholder = "One per line, or separated by spaces or commas"
        ),
        shiny::fileInput(
          ns("gene_file"),
          "Upload a gene list",
          accept = c(".txt", ".csv", ".tsv"),
          buttonLabel = "Choose file"
        ),
        htmltools::div(
          class = "button-row",
          shiny::actionButton(
            ns("add_gene_text"),
            "Add genes",
            class = "btn-primary"
          ),
          shiny::actionButton(
            ns("add_active_gene"),
            "Add current gene",
            class = "btn-outline-primary"
          )
        ),
        shiny::uiOutput(ns("gene_input_status")),
        shiny::uiOutput(ns("gene_set_status")),
        htmltools::p(
          "Choose rows in the Gene set tab to remove individual genes.",
          class = "supporting-copy"
        ),
        shiny::conditionalPanel(
          condition = "output.saved_gene_count > 0",
          ns = ns,
          htmltools::div(
            class = "button-row",
            shiny::actionButton(
              ns("remove_genes"),
              "Remove selected genes",
              class = "btn-outline-secondary"
            ),
            shiny::actionButton(
              ns("clear_gene_set"),
              "Clear set",
              class = "btn-outline-secondary"
            )
          )
        ),
        class = "gene-set-menu"
      ),
      tool_menu(
        "Select cells",
        htmltools::p(
          "Click a cell, drag a lasso, or choose a cluster. Use the plot toolbar to switch to box selection.",
          class = "supporting-copy"
        ),
        shiny::selectizeInput(
          ns("select_cell_ids"),
          "Select individual cells",
          choices = NULL,
          multiple = TRUE,
          options = list(placeholder = "Search public cell IDs")
        ),
        shiny::actionButton(
          ns("apply_cell_ids"),
          "Apply cell selection",
          class = "btn-outline-primary"
        ),
        shiny::selectInput(
          ns("select_cluster"),
          "Select cells by cluster",
          choices = c(
            "Choose cluster" = "",
            stats::setNames(cluster_choices, paste("Cluster", cluster_choices))
          )
        ),
        htmltools::div(
          class = "button-row",
          shiny::actionButton(
            ns("select_cluster_cells"),
            "Select cluster",
            class = "btn-primary"
          ),
          shiny::actionButton(
            ns("clear_selection"),
            "Clear selection",
            class = "btn-outline-secondary"
          )
        )
      ),
      shiny::uiOutput(ns("active_set_scope")),
      tool_menu(
        "Export",
        shiny::uiOutput(ns("export_scope")),
        htmltools::div(
          class = "download-stack",
          shiny::downloadButton(
            ns("download_analysis_genes"),
            "Analysis genes (TXT)"
          ),
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
            "Saved gene set (TXT)",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_summary"),
            "Selection summary (CSV)",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_expression"),
            "Cell expression (RDS)",
            class = "btn-outline-primary btn-sm"
          ),
          shiny::downloadButton(
            ns("download_metadata"),
            "Cell metadata (CSV)",
            class = "btn-outline-primary btn-sm"
          )
        ),
        htmltools::tags$details(
          class = "plot-help",
          htmltools::tags$summary("Use the expression RDS"),
          htmltools::p(
            "A genes × cells shifted sparse matrix and a separate per-cell centering vector. Reconstruct a small selection with:"
          ),
          htmltools::tags$pre(
            "x <- readRDS('espiviz-selected-cell-expression.rds')\nPFlog <- sweep(as.matrix(x$normalized_expression$sparse), 2L, x$normalized_expression$center, '-')"
          )
        ),
        shiny::downloadButton(
          ns("download_view_metadata"),
          "View metadata (JSON)"
        ),
        class = "export-menu"
      ),
      tool_menu(
        "Save view",
        htmltools::p(
          "Links include genes, a whole-cluster selection, plot settings, result tab, and data version. Save a view file for arbitrary cell selections or long gene lists.",
          class = "supporting-copy"
        ),
        shiny::uiOutput(ns("view_link_status")),
        shiny::actionButton(
          ns("copy_view"),
          "Copy view link",
          class = "btn-primary"
        ),
        shiny::downloadButton(ns("download_view"), "Save view (JSON)"),
        shiny::fileInput(
          ns("restore_view"),
          "Restore a saved view",
          accept = ".json"
        )
      )
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      fillable = FALSE,
      fill = FALSE,
      class = "explore-grid",
      bslib::card(
        class = "main-figure-card",
        full_screen = TRUE,
        bslib::card_header(
          htmltools::span("Cell map", class = "figure-heading"),
          shiny::uiOutput(ns("umap_title"), inline = TRUE)
        ),
        htmltools::div(
          role = "region",
          `aria-label` = "Interactive UMAP",
          `aria-describedby` = ns("umap_note"),
          plotly::plotlyOutput(ns("umap"), height = "460px")
        ),
        htmltools::p(
          id = ns("umap_note"),
          class = "figure-note",
          "UMAP is a qualitative projection. Explore local neighborhoods; global distance, cluster area, and point density are not quantitative."
        ),
        htmltools::div(
          class = "selection-summary",
          shiny::uiOutput(ns("selection_snapshot")),
          shiny::uiOutput(ns("umap_legend"))
        )
      ),
      bslib::card(
        class = "results-card",
        full_screen = TRUE,
        htmltools::div(
          class = "result-toolbar result-page-toolbar",
          htmltools::div(
            htmltools::h3("Expression", class = "result-title"),
            shiny::uiOutput(ns("gene_page_status"), inline = TRUE)
          ),
          shiny::conditionalPanel(
            condition = "output.gene_page_count > 1",
            ns = ns,
            htmltools::div(
              class = "gene-pager",
              shiny::actionButton(ns("previous_gene_page"), "Previous"),
              shiny::numericInput(
                ns("gene_page"),
                "Gene page",
                value = 1L,
                min = 1L,
                step = 1L,
                width = "90px"
              ),
              shiny::actionButton(ns("next_gene_page"), "Next")
            )
          )
        ),
        shiny::uiOutput(ns("result_scope")),
        shiny::conditionalPanel(
          condition = "output.summary_is_dot === 'true'",
          ns = ns,
          shiny::selectInput(
            ns("color_scale"),
            "Dot color scale",
            c(
              "Automatic per figure" = "auto",
              "Shared across analysis genes" = "shared"
            ),
            width = "270px"
          )
        ),
        bslib::navset_tab(
          id = ns("explore_results"),
          selected = "By cluster",
          bslib::nav_panel(
            "By cluster",
            htmltools::div(
              class = "result-toolbar",
              htmltools::h3(
                "Expression by final cluster",
                class = "result-title"
              ),
              summary_plot_type_input("cluster_plot_type")
            ),
            shiny::uiOutput(ns("cluster_summary_plot_ui")),
            summary_plot_note("cluster_plot_type"),
            htmltools::tags$details(
              class = "data-disclosure",
              htmltools::tags$summary("View data table"),
              DT::DTOutput(ns("cluster_table"))
            )
          ),
          bslib::nav_panel(
            "Comparison",
            htmltools::div(
              class = "result-toolbar",
              shiny::uiOutput(ns("comparison_heading")),
              summary_plot_type_input("comparison_plot_type")
            ),
            htmltools::div(
              class = "comparison-plot-shell",
              shiny::uiOutput(ns("gene_comparison_plot_ui"))
            ),
            summary_plot_note("comparison_plot_type"),
            htmltools::tags$details(
              class = "data-disclosure",
              htmltools::tags$summary("View data table"),
              DT::DTOutput(ns("comparison_table"))
            )
          ),
          bslib::nav_panel(
            "Cell-level",
            htmltools::h3(
              "log normalized expression distribution by final cluster",
              class = "result-title"
            ),
            htmltools::p(
              paste(
                "Violins show the current gene and optional second plot gene",
                "across final clusters. log normalized expression is centered",
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
                  choices = c("No trend" = "none", gene_pair_group_choices),
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
                  "Larger outlined diamonds are explicitly selected cells. LOESS and its 95% confidence interval describe cells; the interval does not estimate uncertainty among biological replicates."
                ),
                class = "supporting-copy"
              )
            ),
            shiny::conditionalPanel(
              condition = "input.gene_pair_display === 'density'",
              ns = ns,
              htmltools::p(
                paste(
                  "Filled contours show a smoothed probability density, normalized for the chosen group, for cells detected by raw",
                  "count for at least one gene."
                ),
                class = "supporting-copy"
              )
            ),
            shiny::uiOutput(ns("gene_pair_plot_ui"))
          ),
          bslib::nav_panel(
            "By sample",
            htmltools::div(
              class = "result-toolbar",
              htmltools::h3(
                "Expression by biological sample",
                class = "result-title"
              ),
              summary_plot_type_input("sample_plot_type")
            ),
            shiny::uiOutput(ns("sample_summary_plot_ui")),
            htmltools::p(
              "Samples are the biological replicates, ordered by condition.",
              class = "supporting-copy"
            ),
            summary_plot_note("sample_plot_type"),
            htmltools::tags$details(
              class = "data-disclosure",
              htmltools::tags$summary("View data table"),
              DT::DTOutput(ns("sample_table"))
            ),
            htmltools::hr(class = "result-divider"),
            htmltools::h3(
              "Observed cluster composition — all cells",
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
            htmltools::tags$details(
              class = "data-disclosure",
              htmltools::tags$summary("View data table"),
              DT::DTOutput(ns("composition_table"))
            )
          ),
          bslib::nav_panel(
            "Pooled condition",
            htmltools::div(
              class = "result-toolbar",
              htmltools::h3(
                "Expression by condition — pooled cells",
                class = "result-title"
              ),
              summary_plot_type_input("condition_plot_type")
            ),
            htmltools::p(
              paste(
                "This descriptive view pools cells within each condition.",
                "Use By sample to inspect replicate consistency."
              ),
              class = "supporting-copy"
            ),
            shiny::uiOutput(ns("condition_summary_plot_ui")),
            summary_plot_note("condition_plot_type"),
            htmltools::tags$details(
              class = "data-disclosure",
              htmltools::tags$summary("View data table"),
              DT::DTOutput(ns("condition_table"))
            )
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
            htmltools::p(
              "Marker statistics compare the chosen cluster with all remaining cells in the final analysis, independent of the current selection. Detection columns are percentages (source fractions multiplied by 100).",
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
            shiny::uiOutput(ns("empty_gene_set")),
            DT::DTOutput(ns("gene_set_table"))
          )
        ),
        htmltools::tags$details(
          class = "result-downloads data-disclosure",
          htmltools::tags$summary("Download this result"),
          shiny::selectInput(
            ns("export_result"),
            "Result",
            c("Expression by cluster" = "cluster")
          ),
          shiny::uiOutput(ns("result_export_scope")),
          htmltools::div(
            class = "button-row",
            shiny::conditionalPanel(
              condition = "output.can_export_result_figure === 'true'",
              ns = ns,
              shiny::downloadButton(ns("download_result_png"), "Figure (PNG)"),
              shiny::downloadButton(ns("download_result_pdf"), "Figure (PDF)")
            ),
            shiny::downloadButton(
              ns("download_result_data"),
              "Underlying data (CSV)"
            )
          )
        )
      )
    )
  )
}

summary_datatable <- function(data, page_length = 25L) {
  data <- label_result_conditions(data)
  if ("sample" %in% names(data)) {
    data$sample <- sample_label(data$sample)
  }
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
    remaining_mean = "Remaining cells mean log normalized expression",
    remaining_median = "Remaining median log normalized expression",
    remaining_detected_n = "Remaining detected cells",
    remaining_detected_pct = "Remaining cells detected (%)",
    mean_difference = "Selected − remaining mean",
    detection_pp_difference = "Detection difference (percentage points)",
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
    analysis_genes <- shiny::reactive(expression_summary_genes())
    composition_data <- prepare_cluster_composition(bundle)
    current_summary_type <- shiny::reactive({
      key <- switch(
        input$explore_results %||% "By cluster",
        "Comparison" = "comparison",
        "By sample" = "sample",
        "Pooled condition" = "condition",
        "cluster"
      )
      resolved_plot_type(
        input[[paste0(key, "_plot_type")]] %||% "auto",
        length(expression_summary_genes())
      )
    })
    full_scope_summary <- shiny::reactive(summarize_selection(
      bundle,
      state$selected_cells(),
      expression_summary_genes()
    ))
    shared_color_limit <- shiny::reactive({
      if (!identical(input$color_scale, "shared")) {
        return(NULL)
      }
      data <- full_scope_summary()
      expression_color_limit(c(
        data$comparison$selected_mean,
        data$comparison$remaining_mean,
        data$selected_by_cluster$mean_expression,
        data$selected_by_sample$mean_expression,
        data$selected_by_condition$mean_expression
      ))
    })
    output$summary_is_dot <- shiny::renderText(
      if (
        result_uses_gene_page(input$explore_results %||% "By cluster") &&
          current_summary_type() == "dot"
      ) {
        "true"
      } else {
        "false"
      }
    )
    shiny::outputOptions(output, "summary_is_dot", suspendWhenHidden = FALSE)
    output$saved_gene_count <- shiny::renderText(length(state$gene_set()))
    shiny::outputOptions(output, "saved_gene_count", suspendWhenHidden = FALSE)
    output$active_set_scope <- shiny::renderUI({
      n <- length(state$gene_set())
      extra <- setdiff(plot_genes(), state$gene_set())
      htmltools::div(
        class = "active-set-scope",
        role = "status",
        htmltools::strong(
          if (n) {
            paste(
              count_label(n, "saved gene"),
              state$gene_set_name() %||% "",
              sep = " · "
            )
          } else {
            "Current gene workspace"
          }
        ),
        htmltools::span(
          if (n) {
            paste(
              count_label(length(extra), "additional gene"),
              "from the current pair"
            )
          } else {
            "Summaries include the current gene and optional second gene."
          }
        )
      )
    })
    output$result_scope <- shiny::renderUI({
      tab <- input$explore_results %||% "By cluster"
      context <- tab %in% c("Cell-level", "Cluster markers")
      htmltools::p(
        class = "result-scope",
        role = "status",
        if (tab == "Gene set") {
          "Saved gene set · complete list"
        } else {
          selection_description(bundle, state$selected_cells(), context)
        },
        if (tab == "By sample") {
          " · Composition below always shows all cells."
        } else if (tab == "Comparison" && length(state$selected_cells())) {
          " · Compared with all remaining cells."
        } else if (context) {
          " · Independent of gene pagination."
        } else {
          ""
        }
      )
    })
    output$empty_gene_set <- shiny::renderUI({
      if (!length(state$gene_set())) {
        htmltools::div(
          class = "empty-state",
          htmltools::h3("Build a gene set"),
          htmltools::p(
            "No gene set yet. Choose a featured set, paste gene symbols, or add the current gene in the Gene set menu."
          )
        )
      }
    })
    shiny::observe({
      session$sendCustomMessage(
        "control-disabled",
        list(
          id = ns("remove_genes"),
          disabled = !length(input$gene_set_table_rows_selected)
        )
      )
      session$sendCustomMessage(
        "control-disabled",
        list(id = ns("download_gene_set"), disabled = !length(state$gene_set()))
      )
      session$sendCustomMessage(
        "control-disabled",
        list(id = ns("previous_gene_page"), disabled = page_info()$page <= 1L)
      )
      session$sendCustomMessage(
        "control-disabled",
        list(
          id = ns("next_gene_page"),
          disabled = page_info()$page >= page_info()$pages
        )
      )
    })
    shiny::observeEvent(
      input$previous_gene_page,
      shiny::updateNumericInput(
        session,
        "gene_page",
        value = max(1L, page_info()$page - 1L)
      )
    )
    shiny::observeEvent(
      input$next_gene_page,
      shiny::updateNumericInput(
        session,
        "gene_page",
        value = min(page_info()$pages, page_info()$page + 1L)
      )
    )
    shiny::updateSelectizeInput(
      session,
      "select_cell_ids",
      choices = as.character(bundle$cells$cell_id),
      server = TRUE
    )
    shiny::observeEvent(input$apply_cell_ids, {
      state$selected_cells(intersect(
        bundle$cells$cell_id,
        input$select_cell_ids %||% character()
      ))
      show_selection_comparison()
    })
    shiny::observe({
      defaults <- view_option_defaults()
      values <- lapply(names(defaults), function(key) {
        input[[key]] %||% defaults[[key]]
      })
      names(values) <- names(defaults)
      if (!identical(shiny::isolate(state$view_options()), values)) {
        state$view_options(values)
      }
    })
    shiny::observeEvent(
      state$restore(),
      {
        values <- state$restore()
        for (key in setdiff(
          names(values),
          c("secondary_gene", "gene_page", "explore_results")
        )) {
          shiny::updateSelectInput(session, key, selected = values[[key]])
        }
        shiny::updateSelectizeInput(
          session,
          "secondary_gene",
          choices = c(
            "No second gene" = "",
            stats::setNames(universe, universe)
          ),
          selected = values$secondary_gene,
          server = TRUE
        )
        shiny::updateNumericInput(
          session,
          "gene_page",
          value = values$gene_page
        )
        bslib::nav_select(
          "explore_results",
          values$explore_results,
          session = session
        )
      },
      ignoreNULL = TRUE
    )
    shiny::observeEvent(
      state$requested_tab(),
      {
        bslib::nav_select(
          "explore_results",
          state$requested_tab(),
          session = session
        )
        shiny::updateNumericInput(session, "gene_page", value = 1L)
        state$requested_tab(NULL)
      },
      ignoreNULL = TRUE
    )
    saved_view <- shiny::reactive(capture_view_state(bundle, state))
    link_query <- shiny::reactive(view_link_query(saved_view()))
    output$view_link_status <- shiny::renderUI({
      htmltools::p(
        role = "status",
        class = "supporting-copy",
        if (is.null(link_query())) {
          "This selection or gene list needs a saved view file. Save view (JSON) preserves it completely."
        } else {
          "Ready to copy a reproducible view link."
        }
      )
    })
    shiny::observe({
      session$sendCustomMessage(
        "control-disabled",
        list(id = ns("copy_view"), disabled = is.null(link_query()))
      )
    })
    shiny::observeEvent(input$copy_view, {
      shiny::req(link_query())
      session$sendCustomMessage(
        "copy-view",
        list(query = link_query(), statusId = ns("view_link_status"))
      )
    })
    output$download_view <- shiny::downloadHandler(
      filename = function() "espiviz-view.json",
      content = function(file) {
        jsonlite::write_json(
          saved_view(),
          file,
          auto_unbox = TRUE,
          pretty = TRUE,
          null = "null",
          na = "null"
        )
      }
    )
    shiny::observeEvent(input$restore_view, {
      tryCatch(
        {
          if (input$restore_view$size > 2e6) {
            stop("View files must be smaller than 2 MB.")
          }
          value <- jsonlite::read_json(
            input$restore_view$datapath,
            simplifyVector = TRUE
          )
          restore_view_state(value, bundle, state)
          shiny::showNotification("Saved view restored.", type = "message")
        },
        error = function(error) {
          shiny::showNotification(
            conditionMessage(error),
            type = "error",
            duration = NULL
          )
        }
      )
    })
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
      analysis_gene_scope(
        state$active_gene(),
        setdiff(plot_genes(), state$active_gene()),
        state$gene_set()
      )
    })

    plot_gene_data <- shiny::reactive({
      prepare_plot_gene_data(
        bundle,
        plot_genes(),
        state$selected_cells()
      )
    })

    add_text_genes <- function(text) {
      before <- state$gene_set()
      parsed <- append_state_gene_set(state, bundle, text)
      added <- length(setdiff(state$gene_set(), before))
      present <- length(intersect(parsed$genes, before))
      input_message(paste0(
        "Added ",
        count_label(added),
        ". ",
        present,
        " already present.",
        if (length(parsed$missing)) {
          paste0(
            " ",
            length(parsed$missing),
            if (length(parsed$missing) == 1L) {
              " unrecognized symbol"
            } else {
              " unrecognized symbols"
            },
            ": ",
            paste(parsed$missing, collapse = ", "),
            "."
          )
        } else {
          ""
        }
      ))
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
      add_text_genes(state$active_gene())
    })

    shiny::observeEvent(input$use_preset, {
      preset <- input$preset %||% ""
      if (nzchar(preset) && preset %in% names(bundle$featured_gene_sets)) {
        replace_state_gene_set(
          state,
          bundle,
          bundle$featured_gene_sets[[preset]]
        )
        state$gene_set_name(preset)
        input_message(paste(
          "Replaced gene set with",
          preset,
          "·",
          count_label(length(state$gene_set()))
        ))
      }
    })

    shiny::observeEvent(input$clear_gene_set, {
      state$gene_set(character())
      state$gene_set_name(NULL)
      input_message(NULL)
    })

    shiny::observeEvent(input$remove_genes, {
      rows <- input$gene_set_table_rows_selected %||% integer()
      current <- state$gene_set()
      if (length(rows) > 0L && length(current) > 0L) {
        state$gene_set(current[-rows])
        state$gene_set_name(NULL)
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
            paste("expression blend")
        } else {
          plot_genes()[[1L]]
        },
        detection = paste(
          paste(plot_genes(), collapse = " + "),
          "detection"
        ),
        cluster = "Final cluster",
        condition = "Condition"
      )
      htmltools::span(label, class = "figure-subtitle")
    })

    output$umap_legend <- shiny::renderUI({
      if (
        !(input$color_by %||% "expression") %in%
          c(
            "expression",
            "detection"
          ) ||
          length(plot_genes()) != 2L
      ) {
        return(NULL)
      }
      blend_legend_ui(
        plot_genes(),
        mode = input$color_by %||% "expression"
      )
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
      visible_genes <- page_info()$genes
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
            "All cells — no selection applied"
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
        if (current_summary_type() == "dot") 25L else 6L
      )
    })

    page_violin_expression_data <- shiny::reactive({
      prepare_summary_violin_data(
        bundle,
        page_info()$genes
      )
    })

    page_violin_data <- shiny::reactive({
      data <- page_violin_expression_data()
      data$selected <- data$cell_id %in% as.character(state$selected_cells())
      data
    })

    is_dot_summary_plot <- function(plot_type) {
      identical(
        resolved_plot_type(plot_type, length(expression_summary_genes())),
        "dot"
      )
    }

    summary_plot_height <- function(plot_type) {
      if (is_dot_summary_plot(plot_type)) {
        comparison_plot_height(length(page_info()$genes))
      } else {
        summary_violin_plot_height(length(page_info()$genes))
      }
    }

    summary_plot_output_ui <- function(id, plot_type) {
      shiny::plotOutput(
        ns(id),
        height = paste0(summary_plot_height(plot_type), "px")
      )
    }

    make_explore_summary_plot <- function(
      plot_type,
      dot_plot,
      group_by,
      group_label,
      rotate_x = FALSE
    ) {
      if (is_dot_summary_plot(plot_type)) {
        return(dot_plot)
      }
      make_summary_violin_plot(
        page_violin_data(),
        bundle = bundle,
        group_by = group_by,
        group_label = group_label,
        rotate_x = rotate_x
      )
    }

    explore_summary_plot_alt <- function(
      plot_type,
      scope,
      detection_scope = "selected cells"
    ) {
      if (is_dot_summary_plot(plot_type)) {
        return(paste(
          "Dot plot of log normalized expression for",
          count_label(length(page_info()$genes)),
          scope,
          "Color shows mean log normalized expression; dot size shows the",
          "percentage of",
          detection_scope,
          "with detected expression."
        ))
      }
      paste(
        "Violin plots of cell-level log normalized expression for",
        count_label(length(page_info()$genes)),
        scope
      )
    }

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

    output$gene_page_count <- shiny::renderText(
      if (result_uses_gene_page(input$explore_results %||% "By cluster")) {
        page_info()$pages
      } else {
        0L
      }
    )
    shiny::outputOptions(output, "gene_page_count", suspendWhenHidden = FALSE)

    output$gene_page_status <- shiny::renderUI({
      page <- page_info()
      tab <- input$explore_results %||% "By cluster"
      text <- if (result_uses_gene_page(tab)) {
        paste0(
          "Gene page ",
          page$page,
          " of ",
          page$pages,
          " · Showing ",
          length(page$genes),
          " of ",
          count_label(page$total, "analysis gene")
        )
      } else {
        switch(
          tab,
          "Cell-level" = paste(plot_genes(), collapse = " + "),
          "Cluster markers" = "Fixed cluster markers · all cells",
          "Gene set" = count_label(length(state$gene_set()), "saved gene")
        )
      }
      htmltools::span(text, class = "supporting-copy")
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
      summary_plot_output_ui(
        "gene_comparison_plot",
        input$comparison_plot_type
      )
    })

    output$gene_comparison_plot <- shiny::renderPlot(
      {
        plot <- make_explore_summary_plot(
          input$comparison_plot_type,
          make_gene_comparison_plot(
            page_comparison(),
            color_limit = shared_color_limit()
          ),
          group_by = "comparison",
          group_label = NULL
        )
        if (is.null(plot)) {
          return(invisible(NULL))
        }
        plot
      },
      height = function() summary_plot_height(input$comparison_plot_type),
      alt = function() {
        explore_summary_plot_alt(
          input$comparison_plot_type,
          "across the current comparison groups.",
          detection_scope = "cells"
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
          "log normalized expression distribution violin plot for",
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
          "across final clusters. log normalized expression is centered and",
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
          "Choose a second gene above to show the two-gene plot."
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
      summary_plot_output_ui(
        "condition_summary_plot",
        input$condition_plot_type
      )
    })

    output$condition_summary_plot <- shiny::renderPlot(
      {
        make_explore_summary_plot(
          input$condition_plot_type,
          make_group_summary_plot(
            condition_page_summary(),
            group_column = "condition",
            group_label = "Condition",
            color_limit = shared_color_limit()
          ),
          group_by = "condition",
          group_label = "Condition"
        )
      },
      height = function() summary_plot_height(input$condition_plot_type),
      alt = function() {
        explore_summary_plot_alt(
          input$condition_plot_type,
          "across pooled conditions."
        )
      }
    )

    output$sample_summary_plot_ui <- shiny::renderUI({
      summary_plot_output_ui(
        "sample_summary_plot",
        input$sample_plot_type
      )
    })

    output$sample_summary_plot <- shiny::renderPlot(
      {
        make_explore_summary_plot(
          input$sample_plot_type,
          make_group_summary_plot(
            sample_page_summary(),
            group_column = "sample",
            group_label = "Biological sample",
            rotate_x = TRUE,
            color_limit = shared_color_limit()
          ),
          group_by = "sample",
          group_label = "Biological sample",
          rotate_x = TRUE
        )
      },
      height = function() summary_plot_height(input$sample_plot_type),
      alt = function() {
        explore_summary_plot_alt(
          input$sample_plot_type,
          paste(
            "across",
            length(unique(sample_page_summary()$sample)),
            "biological samples ordered by condition."
          )
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
      summary_plot_output_ui(
        "cluster_summary_plot",
        input$cluster_plot_type
      )
    })

    output$cluster_summary_plot <- shiny::renderPlot(
      {
        make_explore_summary_plot(
          input$cluster_plot_type,
          make_group_summary_plot(
            cluster_page_summary(),
            group_column = "cluster",
            group_label = "Final cluster",
            color_limit = shared_color_limit()
          ),
          group_by = "cluster",
          group_label = "Final cluster"
        )
      },
      height = function() summary_plot_height(input$cluster_plot_type),
      alt = function() {
        explore_summary_plot_alt(
          input$cluster_plot_type,
          "across final clusters."
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
      data <- composition_data[,
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
      if (!length(genes)) {
        return(NULL)
      }
      data <- data.frame(Gene = genes, stringsAsFactors = FALSE)
      DT::datatable(
        data,
        rownames = FALSE,
        selection = "multiple",
        class = "compact stripe",
        callback = DT::JS(
          "table.table().node().dataset.rowSelectable = 'true';"
        ),
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
      scientific_datatable(
        marker_table_data(marker_data()),
        rownames = FALSE,
        selection = "single",
        options = list(pageLength = 15L, lengthChange = FALSE, scrollX = TRUE)
      )
    })

    shiny::observeEvent(
      input$marker_table_rows_selected,
      {
        row <- input$marker_table_rows_selected %||% integer()
        data <- marker_data()
        if (length(row) > 0L && row[[1L]] <= nrow(data)) {
          set_state_gene(state, bundle, data$gene[[row[[1L]]]])
        }
      },
      ignoreInit = TRUE
    )

    register_explore_exports(
      input,
      output,
      session,
      bundle,
      state,
      list(
        genes = expression_summary_genes,
        plot_genes = plot_genes,
        gene_data = plot_gene_data,
        page = page_info,
        summary = selection_summary,
        violin_data = page_violin_data,
        marker_overview = marker_overview,
        color_limit = shared_color_limit,
        height = function() summary_plot_height(current_summary_type())
      )
    )

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
          plot = figure_context(
            make_umap_ggplot(
              bundle,
              input$color_by %||% "expression",
              plot_genes(),
              state$selected_cells(),
              gene_data = plot_gene_data()
            ),
            bundle,
            "Final cell map",
            plot_genes(),
            state$selected_cells(),
            context = TRUE
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
          plot = figure_context(
            make_umap_ggplot(
              bundle,
              input$color_by %||% "expression",
              plot_genes(),
              state$selected_cells(),
              gene_data = plot_gene_data()
            ),
            bundle,
            "Final cell map",
            plot_genes(),
            state$selected_cells(),
            context = TRUE
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
        writeLines(state$gene_set(), file, useBytes = TRUE)
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
        utils::write.csv(
          label_result_conditions(data),
          file,
          row.names = FALSE,
          na = ""
        )
      }
    )
  })
}
