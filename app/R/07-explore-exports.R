register_explore_exports <- function(
  input,
  output,
  session,
  bundle,
  state,
  api
) {
  ns <- session$ns
  choices_for_tab <- function(tab) {
    switch(
      tab %||% "By cluster",
      "Comparison" = c("Selected and remaining cells" = "comparison"),
      "Cell-level" = c(
        "Cell distributions" = "cells",
        "Gene-pair plot" = "pair"
      ),
      "By sample" = c(
        "Expression by sample" = "sample",
        "Cluster composition" = "composition"
      ),
      "Pooled condition" = c("Expression by condition" = "condition"),
      "Cluster markers" = c("Cluster markers" = "markers"),
      "Gene set" = c("Saved gene set" = "set"),
      c("Expression by cluster" = "cluster")
    )
  }
  shiny::observeEvent(
    input$explore_results,
    {
      choices <- choices_for_tab(input$explore_results)
      shiny::updateSelectInput(
        session,
        "export_result",
        choices = choices,
        selected = unname(choices[[1L]])
      )
    },
    ignoreNULL = FALSE
  )
  result <- shiny::reactive({
    valid <- unname(choices_for_tab(input$explore_results))
    if ((input$export_result %||% "") %in% valid) {
      input$export_result
    } else {
      valid[[1L]]
    }
  })
  full_summary <- shiny::reactive(summarize_selection(
    bundle,
    state$selected_cells(),
    api$genes()
  ))
  result_data <- shiny::reactive({
    type <- result()
    if (type == "set") {
      return(data.frame(gene = state$gene_set()))
    }
    if (type == "composition") {
      return(prepare_cluster_composition(bundle))
    }
    if (type == "markers") {
      return(api$marker_overview())
    }
    if (type == "cells") {
      return(plot_gene_data_long(api$gene_data()))
    }
    if (type == "pair") {
      data <- api$gene_data()
      if (length(attr(data, "genes")) != 2L) {
        return(data.frame())
      }
      group <- if (identical(input$gene_pair_display, "density")) {
        input$gene_pair_density_group %||% "all"
      } else {
        "all"
      }
      data <- data[gene_pair_scope(data, group)$included, , drop = FALSE]
      names(data)[names(data) == "expression_1"] <- paste0(
        attr(api$gene_data(), "genes")[[1L]],
        "_PFlog"
      )
      names(data)[names(data) == "expression_2"] <- paste0(
        attr(api$gene_data(), "genes")[[2L]],
        "_PFlog"
      )
      return(data)
    }
    switch(
      type,
      comparison = full_summary()$comparison,
      cluster = full_summary()$selected_by_cluster,
      sample = full_summary()$selected_by_sample,
      condition = full_summary()$selected_by_condition
    )
  })
  result_plot <- shiny::reactive({
    type <- result()
    if (type == "set") {
      return(NULL)
    }
    plot <- switch(
      type,
      cells = make_violin_plot(api$gene_data(), bundle),
      pair = make_gene_pair_ggplot(
        api$gene_data(),
        bundle,
        input$gene_pair_display %||% "scatter",
        input$gene_pair_loess_group %||% "all",
        input$gene_pair_density_group %||% "all"
      ),
      composition = make_cluster_composition_plot(
        prepare_cluster_composition(bundle),
        bundle
      ),
      markers = make_marker_overview_plot(api$marker_overview()),
      {
        summary <- api$summary()
        grouping <- if (type == "comparison") "comparison" else type
        data <- switch(
          type,
          comparison = comparison_plot_data(summary$comparison),
          cluster = group_summary_plot_data(
            summary$selected_by_cluster,
            "cluster"
          ),
          sample = group_summary_plot_data(
            summary$selected_by_sample,
            "sample"
          ),
          condition = group_summary_plot_data(
            summary$selected_by_condition,
            "condition"
          )
        )
        plot_type <- input[[paste0(type, "_plot_type")]] %||% "auto"
        if (resolved_plot_type(plot_type, length(api$genes())) == "dot") {
          make_expression_dot_plot(
            data,
            group_label = tools::toTitleCase(type),
            rotate_x = type == "sample",
            color_limit = api$color_limit()
          ) +
            ggplot2::theme(
              legend.box = "horizontal",
              legend.spacing.x = grid::unit(1, "cm"),
              legend.title = ggplot2::element_text(size = 11),
              legend.text = ggplot2::element_text(size = 10)
            )
        } else {
          make_summary_violin_plot(
            api$violin_data(),
            bundle,
            grouping,
            tools::toTitleCase(type)
          )
        }
      }
    )
    if (is.null(plot)) {
      return(NULL)
    }
    context <- type %in% c("cells", "pair", "composition", "markers")
    genes <- if (type %in% c("cells", "pair")) {
      api$plot_genes()
    } else if (type == "markers") {
      unique(api$marker_overview()$gene)
    } else if (type == "composition") {
      character()
    } else {
      api$page()$genes
    }
    note <- if (type %in% c("comparison", "cluster", "sample", "condition")) {
      plot_type <- resolved_plot_type(
        input[[paste0(type, "_plot_type")]] %||% "auto",
        length(api$genes())
      )
      paste(
        "Descriptive summary. Figure:",
        length(genes),
        "of",
        count_label(length(api$genes())),
        "(page",
        paste0(api$page()$page, "/", api$page()$pages, ")."),
        "CSV contains every analysis gene.",
        if (plot_type == "dot") {
          paste(
            "Color: mean log normalized expression;",
            input$color_scale %||% "auto",
            "limits. Dot area: detected %. Zero area: 0%; crosses: no cells."
          )
        } else {
          "Groups below 10 cells: individual points and median. n counts cells; absent samples have 0 selected cells."
        }
      )
    } else {
      ""
    }
    scope_text <- if (type == "pair") {
      scope <- gene_pair_scope(
        api$gene_data(),
        if (identical(input$gene_pair_display, "density")) {
          input$gene_pair_density_group %||% "all"
        } else {
          "all"
        }
      )
      paste(
        count_label(scope$included_n, "included cell"),
        "of",
        count_label(scope$total_n, "cell"),
        "·",
        scope$group_label
      )
    } else if (type %in% c("markers", "composition")) {
      paste(
        "All",
        count_label(nrow(bundle$cells), "cell"),
        "· independent of selection"
      )
    } else {
      NULL
    }
    figure_context(
      plot,
      bundle,
      names(choices_for_tab(input$explore_results))[match(
        type,
        choices_for_tab(input$explore_results)
      )],
      genes,
      state$selected_cells(),
      context,
      note,
      scope_text = scope_text
    )
  })
  output$can_export_result_figure <- shiny::renderText(
    if (is.null(result_plot())) "false" else "true"
  )
  shiny::outputOptions(
    output,
    "can_export_result_figure",
    suspendWhenHidden = FALSE
  )
  output$export_scope <- shiny::renderUI({
    selected <- state$selected_cells()
    n <- if (length(selected)) length(selected) else nrow(bundle$cells)
    htmltools::div(
      class = "scope-note",
      role = "status",
      htmltools::strong(paste(
        count_label(length(api$genes())),
        "·",
        count_label(n, "cell")
      )),
      htmltools::p(
        "Analysis downloads contain current + second + saved-set genes across every page. Selection summaries also report remaining cells."
      )
    )
  })
  output$result_export_scope <- shiny::renderUI({
    type <- result()
    htmltools::p(
      class = "supporting-copy",
      if (type == "set") {
        "Saved gene set · TXT or CSV. Switch to an expression tab to export a figure."
      } else if (type %in% c("cells", "pair", "markers", "composition")) {
        if (type == "pair") {
          paste(
            count_label(nrow(result_data()), "included cell"),
            "·",
            count_label(length(api$plot_genes())),
            "· double-negative cells excluded."
          )
        } else if (type == "cells") {
          "All-cell distributions for the current gene pair, with selected cells outlined."
        } else {
          "All-cell context, independent of the current selection. CSV contains the plotted values."
        }
      } else {
        paste(
          "Figure: gene page",
          api$page()$page,
          "of",
          api$page()$pages,
          "· CSV: every analysis gene."
        )
      }
    )
  })
  output$download_result_data <- shiny::downloadHandler(
    filename = function() paste0("espiviz-", result(), ".csv"),
    content = function(file) {
      utils::write.csv(
        label_result_conditions(result_data()),
        file,
        row.names = FALSE,
        na = ""
      )
    }
  )
  for (format in c("png", "pdf")) {
    local({
      format <- format
      output[[paste0("download_result_", format)]] <- shiny::downloadHandler(
        filename = function() paste0("espiviz-", result(), ".", format),
        content = function(file) {
          height <- if (
            result() %in% c("cluster", "comparison", "sample", "condition")
          ) {
            max(7, api$height() / 90)
          } else if (result() == "markers") {
            max(
              7,
              marker_overview_height(length(unique(
                api$marker_overview()$gene
              ))) /
                90
            )
          } else {
            7
          }
          write_figure(result_plot(), file, format, height = height)
        }
      )
    })
  }
  output$download_analysis_genes <- shiny::downloadHandler(
    filename = function() "espiviz-analysis-genes.txt",
    content = function(file) writeLines(api$genes(), file)
  )
  output$download_view_metadata <- shiny::downloadHandler(
    filename = function() "espiviz-view-metadata.json",
    content = function(file) {
      metadata <- capture_view_state(bundle, state)
      metadata$result <- result()
      metadata$analysis_genes <- api$genes()
      metadata$figure_genes <- if (result() %in% c("cells", "pair")) {
        api$plot_genes()
      } else if (result() == "markers") {
        unique(api$marker_overview()$gene)
      } else if (result() == "composition") {
        character()
      } else if (result() == "set") {
        state$gene_set()
      } else {
        api$page()$genes
      }
      metadata$figure_page <- if (
        result() %in% c("cluster", "comparison", "sample", "condition")
      ) {
        api$page()$page
      } else {
        NULL
      }
      metadata$data_rows <- nrow(result_data())
      metadata$normalization <- "scclrR PFlog: target=auto, log1p=TRUE, center=TRUE. Detection: raw count > 0."
      metadata$blend <- blend_explanation()
      metadata$model <- "Fixed primary condition inference; selections are descriptive. LOESS intervals describe cells, not replicate uncertainty."
      jsonlite::write_json(
        metadata,
        file,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null"
      )
    }
  )
}
