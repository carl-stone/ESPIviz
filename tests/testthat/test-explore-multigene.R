test_that("plot-gene data keeps two canonical genes, PFlog, and raw detection", {
  expect_app_helper("normalize_plot_genes")
  expect_app_helper("prepare_plot_gene_data")
  bundle <- synthetic_bundle()

  expect_identical(
    normalize_plot_genes(bundle, c("egfp", "Glul", "Mcm2"), "Other"),
    c("EGFP", "Glul")
  )
  expect_identical(
    normalize_plot_genes(bundle, character(), "Other"),
    "Other"
  )

  observed <- prepare_plot_gene_data(
    bundle,
    c("Glul", "EGFP"),
    selected_cell_ids = c("cell_1", "not-a-cell")
  )

  expect_identical(attr(observed, "genes"), c("Glul", "EGFP"))
  expect_equal(nrow(observed), nrow(bundle$cells))
  expect_identical(observed$cell_id, bundle$cells$cell_id)
  expect_equal(
    observed$expression_1,
    as.numeric(expression_matrix(bundle, "Glul")[, 1L])
  )
  expect_equal(
    observed$expression_2,
    as.numeric(expression_matrix(bundle, "EGFP")[, 1L])
  )
  expect_identical(
    observed$detected_1,
    unname(bundle$counts[, "Glul"] > 0)
  )
  expect_identical(
    observed$detected_2,
    unname(bundle$counts[, "EGFP"] > 0)
  )
  expect_identical(observed$selected, bundle$cells$cell_id == "cell_1")

  low_detected <- observed[observed$cell_id == "cell_6", , drop = FALSE]
  expect_true(low_detected$detected_2)
  expect_lt(low_detected$expression_2, 0)
})

test_that("two-gene expression blend scales each normalized-expression axis", {
  expect_app_helper("prepare_umap_expression_blend_data")
  expect_app_helper("blend_legend_ui")
  bundle <- synthetic_bundle()
  gene_data <- prepare_plot_gene_data(bundle, c("Glul", "EGFP"))

  gene_data$detected_1 <- FALSE
  gene_data$detected_2 <- FALSE
  observed <- prepare_umap_expression_blend_data(gene_data)

  expect_equal(
    observed$strength_1,
    (observed$expression_1 - min(observed$expression_1)) /
      diff(range(observed$expression_1))
  )
  expect_equal(
    observed$strength_2,
    (observed$expression_2 - min(observed$expression_2)) /
      diff(range(observed$expression_2))
  )
  expect_true(all(grepl("^#[0-9A-F]{6}$", observed$blend_color)))
  expect_identical(
    observed$blend_color[observed$cell_id == "cell_4"],
    unname(blend_palette()[["Neither detected"]])
  )
  expect_true(all(observed$blend_strength >= 0))
  expect_true(all(observed$blend_strength <= 2))

  legend <- htmltools::renderTags(
    blend_legend_ui(c("Glul", "EGFP"), mode = "expression")
  )$html
  for (label in c(
    "Low Glul + low EGFP",
    "High Glul",
    "High EGFP",
    "High Glul + high EGFP"
  )) {
    expect_match(legend, label, fixed = TRUE)
  }
  expect_no_match(legend, "detected", fixed = TRUE)
  expect_match(legend, 'role="list"', fixed = TRUE)
})

test_that("two-gene detection mode uses four raw-count detection classes", {
  expect_app_helper("prepare_umap_detection_data")
  bundle <- synthetic_bundle()
  gene_data <- prepare_plot_gene_data(bundle, c("Glul", "EGFP"))

  observed <- prepare_umap_detection_data(gene_data)
  classes <- stats::setNames(observed$blend_class, observed$cell_id)

  expect_identical(classes[["cell_1"]], "Glul detected")
  expect_identical(classes[["cell_2"]], "EGFP detected")
  expect_identical(classes[["cell_3"]], "Both detected")
  expect_identical(classes[["cell_4"]], "Neither detected")
  expect_identical(classes[["cell_6"]], "EGFP detected")
  expect_identical(
    observed$blend_color[observed$cell_id == "cell_2"],
    observed$blend_color[observed$cell_id == "cell_6"]
  )

  legend <- htmltools::renderTags(
    blend_legend_ui(c("Glul", "EGFP"), mode = "detection")
  )$html
  for (label in c(
    "Neither detected",
    "Glul detected",
    "EGFP detected",
    "Both detected"
  )) {
    expect_match(legend, label, fixed = TRUE)
  }
})

test_that("interactive and download UMAPs separate expression and detection blends", {
  bundle <- synthetic_bundle()
  genes <- c("Glul", "EGFP")

  interactive <- plotly::plotly_build(make_umap_plotly(
    bundle = bundle,
    color_by = "expression",
    gene = genes,
    selected_cell_ids = character(),
    source = "blend_umap_test"
  ))
  cell_traces <- Filter(
    function(trace) {
      identical(trace$type, "scattergl") && length(trace$key %||% character()) > 0L
    },
    interactive$x$data
  )
  expect_length(cell_traces, 1L)
  expect_setequal(as.character(cell_traces[[1L]]$key), bundle$cells$cell_id)
  expect_identical(anyDuplicated(as.character(cell_traces[[1L]]$key)), 0L)
  expect_false(isTRUE(cell_traces[[1L]]$marker$showscale))
  expected_expression <- prepare_umap_expression_blend_data(
    prepare_plot_gene_data(bundle, genes)
  )
  expect_equal(
    as.character(cell_traces[[1L]]$marker$color),
    unname(expected_expression$blend_color)
  )

  download <- make_umap_ggplot(
    bundle = bundle,
    color_by = "expression",
    gene = genes,
    selected_cell_ids = character()
  )
  expect_s3_class(download, "ggplot")
  built <- ggplot2::ggplot_build(download)
  expect_equal(nrow(built$data[[1L]]), nrow(bundle$cells))
  expect_match(download$labels$caption, "low both", fixed = TRUE)
  expect_no_match(download$labels$caption, "detected", fixed = TRUE)

  detection <- plotly::plotly_build(make_umap_plotly(
    bundle = bundle,
    color_by = "detection",
    gene = genes,
    selected_cell_ids = character(),
    source = "detection_umap_test"
  ))
  detection_trace <- Filter(
    function(trace) {
      identical(trace$type, "scattergl") && length(trace$key %||% character()) > 0L
    },
    detection$x$data
  )[[1L]]
  expected_detection <- prepare_umap_detection_data(
    prepare_plot_gene_data(bundle, genes)
  )
  expect_equal(
    as.character(detection_trace$marker$color),
    unname(expected_detection$blend_color)
  )
  expect_true(all(grepl("Raw detected", detection_trace$text, fixed = TRUE)))

  detection_download <- make_umap_ggplot(
    bundle = bundle,
    color_by = "detection",
    gene = genes,
    selected_cell_ids = character()
  )
  expect_match(
    detection_download$labels$caption,
    "Detection uses raw counts",
    fixed = TRUE
  )
})

test_that("cell-level violins retain PFlog values and selected-cell overlays", {
  expect_app_helper("plot_gene_data_long")
  expect_app_helper("make_violin_plot")
  bundle <- synthetic_bundle()
  gene_data <- prepare_plot_gene_data(
    bundle,
    c("Glul", "EGFP"),
    selected_cell_ids = c("cell_1", "cell_3")
  )

  long <- plot_gene_data_long(gene_data)
  expect_equal(nrow(long), 2L * nrow(bundle$cells))
  expect_setequal(as.character(long$gene), c("Glul", "EGFP"))
  expect_equal(sum(long$selected), 4L)
  low_detected <- long[
    long$cell_id == "cell_6" & as.character(long$gene) == "EGFP",
    ,
    drop = FALSE
  ]
  expect_true(low_detected$detected)
  expect_lt(low_detected$expression, 0)

  plot <- make_violin_plot(gene_data, bundle)
  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)
  expect_true(any(vapply(built$data, nrow, integer(1L)) == 4L))
  expect_match(plot$labels$y, "Log normalized expression", fixed = TRUE)
  expect_match(plot$labels$caption, "raw counts", fixed = TRUE)
})

test_that("two-gene PFlog scatter excludes double-negative cells", {
  expect_app_helper("gene_pair_scope")
  expect_app_helper("gene_pair_scope_ui")
  expect_app_helper("make_gene_pair_plotly")
  bundle <- synthetic_bundle()
  gene_data <- prepare_plot_gene_data(
    bundle,
    c("Glul", "EGFP"),
    selected_cell_ids = c("cell_1", "cell_3")
  )

  built <- plotly::plotly_build(make_gene_pair_plotly(
    gene_data,
    bundle,
    source = "pair_scatter_test"
  ))
  traces <- Filter(
    function(trace) identical(trace$type, "scattergl"),
    built$x$data
  )
  keys <- unlist(lapply(traces, function(trace) {
    as.character(trace$key %||% character())
  }), use.names = FALSE)
  scope <- gene_pair_scope(gene_data)
  expect_equal(scope$included_n, 5L)
  expect_equal(scope$excluded_n, 1L)
  expect_equal(scope$excluded_pct, 100 / 6)
  expect_equal(length(keys), scope$included_n)
  expect_identical(anyDuplicated(keys), 0L)
  expect_setequal(keys, bundle$cells$cell_id[scope$included])
  expect_false("cell_4" %in% keys)
  expect_setequal(
    sub("^Cluster ", "", vapply(traces, `[[`, character(1L), "name")),
    c("0", "1", "2")
  )
  expect_true(all(unlist(lapply(traces, function(trace) {
    grepl("Raw detected", trace$text, fixed = TRUE)
  }))))
  expect_true(all(unlist(lapply(traces, function(trace) {
    grepl("Explicitly selected", trace$text, fixed = TRUE)
  }))))

  cell_6_trace <- traces[vapply(traces, function(trace) {
    "cell_6" %in% as.character(trace$key %||% character())
  }, logical(1L))][[1L]]
  cell_6_index <- match("cell_6", as.character(cell_6_trace$key))
  expect_lt(as.numeric(cell_6_trace$y[[cell_6_index]]), 0)
  expect_match(
    built$x$layout$xaxis$title$text,
    "Glul log normalized expression",
    fixed = TRUE
  )
  expect_match(
    built$x$layout$yaxis$title$text,
    "EGFP log normalized expression",
    fixed = TRUE
  )
  scope_html <- htmltools::renderTags(gene_pair_scope_ui(scope))$html
  expect_match(scope_html, "Showing 5 of 6 cells", fixed = TRUE)
  expect_match(scope_html, "artificial diagonal", fixed = TRUE)

  cluster_scope <- gene_pair_scope(gene_data, group = "1")
  cluster_scope_html <- htmltools::renderTags(
    gene_pair_scope_ui(cluster_scope)
  )$html
  expect_match(
    cluster_scope_html,
    "Showing 1 of 2 cells in Cluster 1",
    fixed = TRUE
  )

  single <- prepare_plot_gene_data(bundle, "Glul")
  expect_null(make_gene_pair_plotly(single, bundle, "single_gene_test"))
})

test_that("gene-pair loess trends use all included cells or one cluster", {
  expect_app_helper("prepare_gene_pair_loess")
  x <- seq(-2, 2, length.out = 60L)
  gene_data <- data.frame(
    cell_id = paste0("cell_", seq_along(x)),
    cluster = rep(c("0", "1"), each = 30L),
    expression_1 = x,
    expression_2 = 0.5 * x + sin(x) + rep(c(-0.1, 0.1), each = 30L),
    detected_1 = TRUE,
    detected_2 = TRUE,
    selected = FALSE,
    stringsAsFactors = FALSE
  )
  attr(gene_data, "genes") <- c("GeneA", "GeneB")

  all_cells <- prepare_gene_pair_loess(gene_data, group = "all")
  cluster_one <- prepare_gene_pair_loess(gene_data, group = "1")

  expect_identical(all_cells$group_label[[1L]], "All cells")
  expect_identical(unique(all_cells$cell_count), 60L)
  expect_identical(cluster_one$group_label[[1L]], "Cluster 1")
  expect_identical(unique(cluster_one$cell_count), 30L)
  expect_true(all(c("x", "fit", "lower", "upper") %in% names(all_cells)))
  expect_true(all(is.finite(all_cells$x)))
  expect_true(all(is.finite(all_cells$fit)))
  expect_true(all(all_cells$lower <= all_cells$fit))
  expect_true(all(all_cells$fit <= all_cells$upper))
})

test_that("gene-pair density prepares one smooth grid for the requested cells", {
  expect_app_helper("prepare_gene_pair_density")
  set.seed(1L)
  x <- seq(-2, 2, length.out = 60L)
  gene_data <- data.frame(
    cell_id = paste0("cell_", seq_along(x)),
    cluster = rep(c("0", "1"), each = 30L),
    expression_1 = x,
    expression_2 = x^2 + stats::rnorm(length(x), sd = 0.1),
    detected_1 = TRUE,
    detected_2 = TRUE,
    selected = FALSE,
    stringsAsFactors = FALSE
  )
  attr(gene_data, "genes") <- c("GeneA", "GeneB")

  observed <- prepare_gene_pair_density(gene_data, group = "0", grid_n = 40L)

  expect_identical(observed$group_label, "Cluster 0")
  expect_identical(observed$cell_count, 30L)
  expect_length(observed$x, 40L)
  expect_length(observed$y, 40L)
  expect_identical(dim(observed$z), c(40L, 40L))
  expect_true(all(is.finite(observed$z)))
  expect_true(all(observed$z >= 0))
})

test_that("gene-pair display switches between scatter and density traces", {
  x <- seq(-2, 2, length.out = 60L)
  gene_data <- data.frame(
    cell_id = paste0("cell_", seq_along(x)),
    cluster = rep(c("0", "1"), each = 30L),
    expression_1 = x,
    expression_2 = 0.5 * x + sin(x),
    detected_1 = TRUE,
    detected_2 = TRUE,
    selected = FALSE,
    stringsAsFactors = FALSE
  )
  attr(gene_data, "genes") <- c("GeneA", "GeneB")
  bundle <- synthetic_bundle()

  scatter <- plotly::plotly_build(make_gene_pair_plotly(
    gene_data,
    bundle,
    source = "pair_scatter_options_test",
    display = "scatter",
    loess_group = "all"
  ))
  density <- plotly::plotly_build(make_gene_pair_plotly(
    gene_data,
    bundle,
    source = "pair_density_options_test",
    display = "density",
    density_group = "all"
  ))
  scatter_names <- vapply(
    scatter$x$data,
    function(trace) trace$name %||% "",
    character(1L)
  )

  expect_true(any(vapply(scatter$x$data, function(trace) {
    identical(trace$type, "scattergl")
  }, logical(1L))))
  expect_true(any(grepl("Loess trend", scatter_names, fixed = TRUE)))
  expect_true(any(grepl("95% confidence ribbon", scatter_names, fixed = TRUE)))
  expect_false(any(vapply(density$x$data, function(trace) {
    identical(trace$type, "scattergl")
  }, logical(1L))))
  expect_true(any(vapply(density$x$data, function(trace) {
    identical(trace$type, "contour")
  }, logical(1L))))
  expect_s3_class(
    make_gene_pair_plotly(
      gene_data,
      bundle,
      source = "pair_missing_cluster_test",
      display = "scatter",
      loess_group = "missing"
    ),
    "plotly"
  )
  expect_s3_class(
    make_gene_pair_plotly(
      gene_data,
      bundle,
      source = "density_missing_cluster_test",
      display = "density",
      density_group = "missing"
    ),
    "plotly"
  )
  missing_trend <- plotly::plotly_build(make_gene_pair_plotly(
    gene_data,
    bundle,
    source = "pair_missing_cluster_annotation_test",
    display = "scatter",
    loess_group = "missing"
  ))
  expect_match(
    missing_trend$x$layout$annotations[[1L]]$text,
    "Loess trend unavailable for Cluster missing.",
    fixed = TRUE
  )
})

test_that("selection snapshot is compact and uses raw-count detection", {
  expect_app_helper("prepare_selection_snapshot")
  expect_app_helper("selection_snapshot_ui")
  bundle <- synthetic_bundle()
  gene_data <- prepare_plot_gene_data(
    bundle,
    c("Glul", "EGFP"),
    selected_cell_ids = c("cell_1", "cell_3")
  )

  observed <- prepare_selection_snapshot(gene_data)
  expect_true(observed$explicit_selection)
  expect_equal(observed$selected_n, 2L)
  expect_equal(observed$total_n, 6L)
  expect_equal(observed$selected_pct, 100 * 2 / 6)
  expect_equal(observed$cluster_n, 2L)
  expect_equal(observed$sample_n, 2L)

  glul <- observed$genes[observed$genes$gene == "Glul", , drop = FALSE]
  expected_glul <- expression_matrix(bundle, "Glul")[c(1L, 3L), 1L]
  expect_equal(glul$mean_pflog, mean(expected_glul))
  expect_equal(glul$detected_n, 2L)
  expect_equal(glul$detected_pct, 100)
  egfp <- observed$genes[observed$genes$gene == "EGFP", , drop = FALSE]
  expect_equal(egfp$detected_n, 1L)
  expect_equal(egfp$detected_pct, 50)

  html <- htmltools::renderTags(selection_snapshot_ui(observed))$html
  expect_match(html, "2 of 6 cells", fixed = TRUE)
  expect_match(html, "Mean log normalized expression", fixed = TRUE)
  expect_match(html, "detected by raw count", fixed = TRUE)
  expect_match(html, 'aria-live="polite"', fixed = TRUE)

  all_cells <- prepare_selection_snapshot(prepare_plot_gene_data(bundle, "Glul"))
  expect_false(all_cells$explicit_selection)
  expect_equal(all_cells$selected_n, nrow(bundle$cells))
  expect_equal(all_cells$genes$detected_n, 3L)
})

test_that("comparison table keeps only interpretable columns with plain headers", {
  expect_app_helper("prepare_comparison_table")
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[1:2]
  comparison <- summarize_selection(
    bundle,
    selected,
    c("Glul", "EGFP")
  )$comparison

  concise <- prepare_comparison_table(comparison, explicit_selection = TRUE)
  expect_identical(
    names(concise),
    c(
      "gene",
      "selected_mean",
      "selected_detected_pct",
      "remaining_mean",
      "remaining_detected_pct",
      "mean_difference",
      "detection_pp_difference"
    )
  )
  expect_false(any(grepl("median|ratio|_n$", names(concise))))
  displayed <- summary_datatable(concise)
  expect_identical(
    names(displayed$x$data),
    c(
      "Gene",
      "Selected mean log normalized expression",
      "Selected detected (%)",
      "Other cells mean log normalized expression",
      "Other cells detected (%)",
      "Mean difference",
      "Detection difference (pp)"
    )
  )

  all_cells <- prepare_comparison_table(
    summarize_selection(bundle, character(), "Glul")$comparison,
    explicit_selection = FALSE
  )
  expect_identical(names(all_cells), c("gene", "mean_expression", "detected_pct"))
  expect_identical(
    names(summary_datatable(all_cells)$x$data),
    c("Gene", "Mean log normalized expression", "Detected (%)")
  )
})

test_that("Explore exposes a compact two-gene cell-level workflow", {
  html <- htmltools::renderTags(
    explore_ui("explore_test", synthetic_bundle())
  )$html

  for (label in c(
    "Second plot gene (optional)",
    "Current gene is always the first plot gene",
    "Selected-cell summary",
    "Cell-level",
    "Log normalized expression distribution by final cluster",
    "Two-gene log normalized expression"
  )) {
    expect_match(html, label, fixed = TRUE)
  }
  for (id in c(
    "secondary_gene",
    "selection_snapshot",
    "umap_legend",
    "violin_plot_ui",
    "gene_pair_display",
    "gene_pair_loess_group",
    "gene_pair_density_group",
    "gene_pair_plot_ui"
  )) {
    expect_match(html, paste0("explore_test-", id), fixed = TRUE)
  }
  for (label in c("Display", "Scatter plot", "Density plot", "Loess trend")) {
    expect_match(html, label, fixed = TRUE)
  }
  expect_match(html, 'class="umap-side-rail"', fixed = TRUE)
})

test_that("optional second plot gene does not mutate global gene state", {
  bundle <- synthetic_bundle()
  state <- new_app_state(bundle)

  shiny::testServer(
    explore_server,
    args = list(bundle = bundle, state = state),
    {
      session$setInputs(
        active_gene = "Glul",
        secondary_gene = "EGFP",
        color_by = "expression"
      )
      expect_identical(state$active_gene(), "Glul")
      expect_identical(state$gene_set(), character())
      title_html <- htmltools::renderTags(output$umap_title)$html
      snapshot_html <- htmltools::renderTags(output$selection_snapshot)$html
      pair_html <- htmltools::renderTags(output$gene_pair_plot_ui)$html
      expect_match(title_html, "Glul + EGFP expression blend", fixed = TRUE)
      expect_match(snapshot_html, "Glul", fixed = TRUE)
      expect_match(snapshot_html, "EGFP", fixed = TRUE)
      expect_match(
        pair_html,
        paste(
          "Two-gene log normalized expression scatter comparing Glul and EGFP",
          "for cells detected for at least one gene"
        ),
        fixed = TRUE
      )
      expect_match(pair_html, "Excluded 1 double-negative cell", fixed = TRUE)

      session$setInputs(
        gene_pair_display = "density",
        gene_pair_density_group = "1"
      )
      density_html <- htmltools::renderTags(output$gene_pair_plot_ui)$html
      expect_match(
        density_html,
        "Two-gene log normalized expression density plot comparing",
        fixed = TRUE
      )
      expect_match(
        density_html,
        "Showing 1 of 2 cells in Cluster 1",
        fixed = TRUE
      )

      session$setInputs(active_gene = "Mcm2")
      expect_identical(state$active_gene(), "Mcm2")
      expect_identical(state$gene_set(), character())
      expect_match(
        htmltools::renderTags(output$umap_title)$html,
        "Mcm2 + EGFP expression blend",
        fixed = TRUE
      )

      session$setInputs(color_by = "detection")
      expect_match(
        htmltools::renderTags(output$umap_title)$html,
        "Mcm2 + EGFP detection",
        fixed = TRUE
      )
    }
  )
})

test_that("manual plot-gene pairs populate every expression summary", {
  bundle <- synthetic_bundle()
  state <- new_app_state(bundle)

  shiny::testServer(
    explore_server,
    args = list(bundle = bundle, state = state),
    {
      expect_expression_summary_genes <- function(expected) {
        expect_identical(page_info()$genes, expected)
        expect_identical(page_comparison()$gene, expected)
        expect_setequal(unique(sample_page_summary()$gene), expected)
        expect_setequal(unique(condition_page_summary()$gene), expected)
        expect_setequal(unique(cluster_page_summary()$gene), expected)
      }

      session$setInputs(
        active_gene = "Glul",
        secondary_gene = "EGFP"
      )

      expect_expression_summary_genes(c("Glul", "EGFP"))

      state$gene_set(c("Mcm2", "Other"))
      session$flushReact()
      expect_expression_summary_genes(c("Glul", "EGFP", "Mcm2", "Other"))
    }
  )
})

test_that("cell selections reuse prepared violin expression data", {
  bundle <- synthetic_bundle()
  state <- new_app_state(bundle)
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  original_prepare <- prepare_summary_violin_data
  rlang::local_bindings(
    prepare_summary_violin_data = function(...) {
      calls$count <- calls$count + 1L
      original_prepare(...)
    },
    .env = environment(explore_server)
  )

  shiny::testServer(
    explore_server,
    args = list(bundle = bundle, state = state),
    {
      initial <- page_violin_data()
      initial_calls <- calls$count

      state$selected_cells(bundle$cells$cell_id[[1L]])
      session$flushReact()
      selected <- page_violin_data()

      expect_gt(initial_calls, 0L)
      expect_identical(calls$count, initial_calls)
      expect_false(any(initial$selected))
      expect_identical(sum(selected$selected), length(unique(selected$gene)))
    }
  )
})
