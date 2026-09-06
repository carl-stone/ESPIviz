result_filter_ui <- function(ns, pathways = FALSE) {
  htmltools::div(
    class = "result-filters",
    shiny::textInput(
      ns("filter_search"),
      if (pathways) "Find a term" else "Find a gene",
      placeholder = if (pathways) "Search GO terms" else "Search gene symbols"
    ),
    if (pathways) {
      shiny::selectInput(
        ns("filter_method"),
        "Method",
        c("All methods" = "all", "GSEA", "ORA")
      )
    },
    shiny::selectInput(
      ns("filter_direction"),
      "Direction",
      c(
        "All directions" = "all",
        "p27CKO" = "Control",
        "p27CKO + E-Stim" = "E-Stim"
      )
    ),
    shiny::selectInput(
      ns("filter_fdr"),
      "Adjusted P",
      c(
        "All results" = "all",
        "FDR ≤ 0.05" = "significant",
        "FDR > 0.05" = "above",
        "Unavailable" = "missing"
      )
    ),
    shiny::actionButton(
      ns("reset_filters"),
      "Reset",
      class = "btn-outline-secondary"
    )
  )
}

reset_result_filters <- function(session, pathways = FALSE) {
  shiny::updateTextInput(session, "filter_search", value = "")
  shiny::updateSelectInput(session, "filter_direction", selected = "all")
  shiny::updateSelectInput(session, "filter_fdr", selected = "all")
  if (pathways) {
    shiny::updateSelectInput(session, "filter_method", selected = "all")
  }
}

filter_result_rows <- function(
  data,
  direction = "all",
  fdr = "all",
  method = "all",
  search = ""
) {
  pathway <- "pathway_id" %in% names(data)
  p <- if (pathway) data$p_adjust else data$padj
  d <- if (pathway) {
    as.character(data$direction)
  } else {
    ifelse(
      data$log2FoldChange > 0,
      "E-Stim",
      ifelse(data$log2FoldChange < 0, "Control", "None")
    )
  }
  keep <- rep(TRUE, nrow(data))
  if (!is.null(direction) && direction != "all") {
    keep <- keep & !is.na(d) & d == direction
  }
  if (identical(fdr, "significant")) {
    keep <- keep & !is.na(p) & p <= 0.05
  }
  if (identical(fdr, "above")) {
    keep <- keep & !is.na(p) & p > 0.05
  }
  if (identical(fdr, "missing")) {
    keep <- keep & is.na(p)
  }
  if (pathway && !is.null(method) && method != "all") {
    keep <- keep & as.character(data$source) == method
  }
  text <- if (pathway) paste(data$label, data$pathway_id) else data$gene
  if (length(search) == 1L && nzchar(search)) {
    keep <- keep & grepl(tolower(search), tolower(text), fixed = TRUE)
  }
  data <- data[which(keep), , drop = FALSE]
  if (pathway) {
    data[
      order(data$p_adjust, data$label, data$pathway_id, na.last = TRUE),
      ,
      drop = FALSE
    ]
  } else {
    data[
      order(data$padj, -abs(data$log2FoldChange), data$gene, na.last = TRUE),
      ,
      drop = FALSE
    ]
  }
}

result_count_ui <- function(shown, total) {
  htmltools::p(
    class = "result-count",
    role = "status",
    paste(
      "Showing",
      format(shown, big.mark = ","),
      "of",
      format(total, big.mark = ","),
      "results"
    )
  )
}

pathway_table_data <- function(data) {
  data.frame(
    Term = data$label,
    Method = as.character(data$source),
    Direction = condition_label(data$direction),
    `Adjusted P` = data$p_adjust,
    NES = ifelse(as.character(data$source) == "GSEA", data$score, NA_real_),
    `Fold enrichment` = ifelse(
      as.character(data$source) == "ORA",
      data$score,
      NA_real_
    ),
    `Source set genes` = ifelse(
      as.character(data$source) == "GSEA",
      data$gene_count,
      NA_integer_
    ),
    `Overlapping genes` = ifelse(
      as.character(data$source) == "ORA",
      data$gene_count,
      NA_integer_
    ),
    check.names = FALSE
  )
}

marker_table_data <- function(data) {
  data.frame(
    Gene = data$gene,
    Rank = data$rank,
    `Average log2 fold change` = data$avg_log2FC,
    `Chosen cluster detected (%)` = 100 * data$pct.1,
    `Remaining cells detected (%)` = 100 * data$pct.2,
    `Adjusted P` = data$p_val_adj,
    check.names = FALSE
  )
}

scientific_datatable <- function(
  data,
  rownames = FALSE,
  selection = "single",
  options = list(),
  ...
) {
  table <- DT::datatable(
    data,
    rownames = rownames,
    escape = TRUE,
    selection = selection,
    class = "stripe scientific-table",
    callback = DT::JS(paste0(
      "table.table().node().dataset.rowSelectable = '",
      if (identical(selection, "none")) "false" else "true",
      "';"
    )),
    options = utils::modifyList(
      list(
        pageLength = 15L,
        lengthChange = FALSE,
        scrollX = TRUE,
        order = list(),
        dom = "tip",
        language = list(
          emptyTable = "No matching results. Change a filter or select Reset."
        )
      ),
      options
    ),
    ...
  )
  numeric_columns <- names(data)[vapply(data, is.numeric, logical(1L))]
  probabilities <- intersect(numeric_columns, c("Adjusted P", "P value"))
  decimals <- setdiff(
    numeric_columns,
    c(probabilities, "Rank", "Source set genes", "Overlapping genes")
  )
  if (length(probabilities)) {
    table <- DT::formatSignif(table, probabilities, digits = 3L)
  }
  if (length(decimals)) {
    table <- DT::formatRound(table, decimals, digits = 2L)
  }
  table
}

pathway_export_height <- function(pathways) {
  labels <- prepare_pathway_plot_data(pathways)$plot_label
  lines <- sum(vapply(
    labels,
    function(label) length(strwrap(label, 49L)),
    integer(1L)
  ))
  max(11, 3.2 + (lines + 0.6 * length(labels)) * 0.15)
}
