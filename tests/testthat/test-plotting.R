test_that("comparison dot-plot data exposes expression and detection by group", {
  expect_app_helper("comparison_plot_data")
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[1:2]
  comparison <- summarize_selection(
    bundle,
    selected,
    c("Glul", "EGFP")
  )$comparison

  plotted <- comparison_plot_data(comparison)

  expect_setequal(unique(plotted$group), c("Selected cells", "Remaining cells"))
  expect_true(all(c("mean_expression", "detected_pct") %in% names(plotted)))
  expect_equal(nrow(plotted), 2L * 2L)

  selected_glul <- plotted[
    plotted$gene == "Glul" &
      plotted$group == "Selected cells",
    ,
    drop = FALSE
  ]
  expected <- comparison$selected_detected_pct[comparison$gene == "Glul"]
  expect_equal(selected_glul$detected_pct, expected)

  all_cells <- summarize_selection(bundle, character(), "Glul")$comparison
  plotted_all <- comparison_plot_data(all_cells)
  expect_identical(unique(plotted_all$group), "All cells")
  expect_false(any(is.na(plotted_all$mean_expression)))
  expect_false(any(is.na(plotted_all$detected_pct)))
})

test_that("grouped dot-plot data keeps cluster and condition summaries intact", {
  expect_app_helper("group_summary_plot_data")
  bundle <- synthetic_bundle()
  summary <- summarize_selection(bundle, character(), c("Glul", "EGFP"))

  by_cluster <- group_summary_plot_data(summary$selected_by_cluster, "cluster")
  expect_setequal(by_cluster$group, c("0", "1", "2"))
  expect_identical(levels(by_cluster$group), c("0", "1", "2"))
  expect_setequal(by_cluster$gene, c("Glul", "EGFP"))
  expect_equal(nrow(by_cluster), nrow(summary$selected_by_cluster))

  target <- by_cluster[by_cluster$group == "0" & by_cluster$gene == "Glul", ]
  source <- summary$selected_by_cluster[
    as.character(summary$selected_by_cluster$cluster) == "0" &
      summary$selected_by_cluster$gene == "Glul",
    ,
    drop = FALSE
  ]
  expect_equal(target$mean_expression, source$mean_expression)
  expect_equal(target$detected_pct, source$detected_pct)

  by_condition <- group_summary_plot_data(
    summary$selected_by_condition,
    "condition"
  )
  expect_setequal(by_condition$group, levels(bundle$cells$condition))
  expect_identical(levels(by_condition$group), levels(bundle$cells$condition))
})

test_that("selection and grouped expression dot plots build across edge cases", {
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[1:2]
  comparisons <- list(
    summarize_selection(bundle, character(), "Glul")$comparison,
    summarize_selection(bundle, selected, "Glul")$comparison,
    summarize_selection(bundle, selected, "ZeroGene")$comparison
  )

  for (comparison in comparisons) {
    expected_rows <- nrow(comparison_plot_data(comparison))
    built <- ggplot2::ggplot_build(make_gene_comparison_plot(comparison))
    expect_equal(nrow(built$data[[1L]]), expected_rows)
    expect_true(all(is.finite(built$data[[1L]]$x)))
    expect_true(all(is.finite(built$data[[1L]]$y)))
  }

  summary <- summarize_selection(bundle, character(), c("Glul", "EGFP"))
  cluster_built <- ggplot2::ggplot_build(make_group_summary_plot(
    summary$selected_by_cluster,
    "cluster",
    "Final cluster"
  ))
  condition_built <- ggplot2::ggplot_build(make_group_summary_plot(
    summary$selected_by_condition,
    "condition",
    "Condition"
  ))
  expect_equal(
    nrow(cluster_built$data[[1L]]),
    nrow(summary$selected_by_cluster)
  )
  expect_equal(
    nrow(condition_built$data[[1L]]),
    nrow(summary$selected_by_condition)
  )
})

test_that("dot-plot color scales start at zero without failing on all-zero genes", {
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[1:2]
  comparison <- summarize_selection(bundle, selected, "Glul")$comparison
  plotted <- comparison_plot_data(comparison)
  plot <- make_gene_comparison_plot(comparison)
  color_scale <- plot$scales$get_scales("colour")

  expect_equal(color_scale$limits, c(0, max(plotted$mean_expression)))

  zero <- summarize_selection(bundle, selected, "ZeroGene")$comparison
  zero_plot <- make_gene_comparison_plot(zero)
  zero_scale <- zero_plot$scales$get_scales("colour")
  expect_equal(zero_scale$limits, c(0, 1))
  expect_silent(ggplot2::ggplot_build(zero_plot))
})

test_that("comparison plot height stays visible and scales through one page", {
  expect_app_helper("comparison_plot_height")

  expect_identical(comparison_plot_height(0L), 300L)
  expect_identical(comparison_plot_height(1L), 300L)
  expect_identical(comparison_plot_height(10L), 380L)
  expect_identical(comparison_plot_height(25L), 680L)
  expect_identical(comparison_plot_height(200L), 680L)
})

test_that("UMAP point styling adapts to the production cell density", {
  expect_app_helper("umap_point_style")

  small <- umap_point_style(100L)
  production <- umap_point_style(3456L)

  expect_gt(small$interactive_size, production$interactive_size)
  expect_gt(small$interactive_opacity, production$interactive_opacity)
  expect_lte(production$interactive_size, 3)
  expect_lte(production$interactive_opacity, 0.6)
  expect_lte(production$static_size, 1)
})

test_that("UMAP selection outlines stay legible for whole clusters", {
  expect_app_helper("selection_overlay_style")
  point_style <- umap_point_style(3456L)
  clicked <- selection_overlay_style(3456L, 1L)
  cluster <- selection_overlay_style(3456L, 968L)

  expect_gt(clicked$interactive_size, cluster$interactive_size)
  expect_gt(clicked$interactive_line_width, cluster$interactive_line_width)
  expect_lte(cluster$interactive_size, point_style$interactive_size + 1)
  expect_lte(cluster$interactive_line_width, 0.6)
  expect_lte(cluster$static_size, point_style$static_size + 0.35)
  expect_lte(cluster$static_stroke, 0.25)
})

test_that("gene-expression UMAP draws high-expressing cells last", {
  bundle <- synthetic_bundle()
  plotted <- umap_plot_data(bundle, "expression", "Glul")

  expect_true(all(diff(plotted$color_value) >= 0))
  expect_setequal(plotted$cell_id, bundle$cells$cell_id)

  clusters <- umap_plot_data(bundle, "cluster", "Glul")
  expect_identical(clusters$cell_id, bundle$cells$cell_id)
})

test_that("cluster labels use one on-data position per cluster", {
  expect_app_helper("cluster_label_positions")
  bundle <- synthetic_bundle()

  labels <- cluster_label_positions(bundle$cells)

  expect_setequal(labels$cluster, unique(as.character(bundle$cells$cluster)))
  expect_equal(nrow(labels), length(unique(bundle$cells$cluster)))
  expect_true(all(is.finite(labels$umap_1)))
  expect_true(all(is.finite(labels$umap_2)))
  expect_true(all(
    paste(labels$umap_1, labels$umap_2) %in%
      paste(bundle$cells$umap_1, bundle$cells$umap_2)
  ))

  singleton <- bundle$cells[1L, , drop = FALSE]
  expect_equal(
    cluster_label_positions(singleton)[c("umap_1", "umap_2")],
    singleton[c("umap_1", "umap_2")]
  )
})

test_that("cluster UMAP labels appear interactively and in downloads", {
  bundle <- synthetic_bundle()
  expected <- sort(unique(as.character(bundle$cells$cluster)))

  interactive <- plotly::plotly_build(make_umap_plotly(
    bundle,
    color_by = "cluster",
    gene = "Glul",
    selected_cell_ids = character(),
    source = "cluster_label_test"
  ))
  text_traces <- Filter(
    function(trace) grepl("text", trace$mode %||% "", fixed = TRUE),
    interactive$x$data
  )
  expect_length(text_traces, 1L)
  expect_setequal(as.character(text_traces[[1L]]$text), expected)
  point_traces <- Filter(
    function(trace) identical(trace$type, "scattergl"),
    interactive$x$data
  )
  expect_true(all(vapply(
    point_traces,
    function(trace) identical(trace$showlegend, FALSE),
    logical(1)
  )))

  download <- make_umap_ggplot(
    bundle,
    color_by = "cluster",
    gene = "Glul",
    selected_cell_ids = character()
  )
  label_layers <- vapply(
    download$layers,
    function(layer) inherits(layer$geom, "GeomLabel"),
    logical(1)
  )
  expect_equal(sum(label_layers), 1L)
})

test_that("interactive UMAP traces preserve every cell key and selected key", {
  bundle <- synthetic_bundle()

  for (color_by in c("expression", "condition", "cluster")) {
    built <- plotly::plotly_build(make_umap_plotly(
      bundle,
      color_by = color_by,
      gene = "Glul",
      selected_cell_ids = character(),
      source = paste0("key_test_", color_by)
    ))
    base_traces <- Filter(
      function(trace) identical(trace$type, "scattergl"),
      built$x$data
    )
    keys <- unlist(lapply(base_traces, `[[`, "key"), use.names = FALSE)
    expect_identical(length(keys), nrow(bundle$cells))
    expect_identical(anyDuplicated(keys), 0L)
    expect_setequal(keys, bundle$cells$cell_id)
  }

  selected <- bundle$cells$cell_id[c(2, 5)]
  built <- plotly::plotly_build(make_umap_plotly(
    bundle,
    color_by = "expression",
    gene = "Glul",
    selected_cell_ids = selected,
    source = "selected_key_test"
  ))
  selected_trace <- Filter(
    function(trace) {
      text <- as.character(trace$text %||% character())
      length(text) > 0L && all(startsWith(text, "Selected cell:"))
    },
    built$x$data
  )
  expect_length(selected_trace, 1L)
  expect_identical(as.character(selected_trace[[1L]]$key), selected)
})

test_that("whole-cluster selection returns exact cell IDs in bundle order", {
  expect_app_helper("cell_ids_for_cluster")
  bundle <- synthetic_bundle()

  expect_identical(
    cell_ids_for_cluster(bundle, "1"),
    bundle$cells$cell_id[as.character(bundle$cells$cluster) == "1"]
  )
  expect_identical(
    cell_ids_for_cluster(bundle, factor("2")),
    bundle$cells$cell_id[as.character(bundle$cells$cluster) == "2"]
  )
  expect_identical(cell_ids_for_cluster(bundle, "not-a-cluster"), character())
  expect_identical(cell_ids_for_cluster(bundle, character()), character())
})

test_that("whole-cluster selection replaces and persists until explicitly cleared", {
  bundle <- synthetic_bundle()
  state <- new_app_state(bundle)
  expected <- cell_ids_for_cluster(bundle, "1")

  shiny::testServer(
    explore_server,
    args = list(bundle = bundle, state = state),
    {
      session$setInputs(select_cluster = "", select_cluster_cells = 0L)
      session$setInputs(select_cluster = "1", select_cluster_cells = 1L)
      expect_identical(state$selected_cells(), expected)

      session$setInputs(active_gene = "EGFP", color_by = "condition")
      expect_identical(state$selected_cells(), expected)

      deselect_input <- paste0(
        "plotly_deselect-",
        session$ns("umap_source")
      )
      do.call(
        session$rootScope()$setInputs,
        stats::setNames(list("{}"), deselect_input)
      )
      expect_identical(state$selected_cells(), expected)

      session$setInputs(clear_selection = 1L)
      expect_identical(state$selected_cells(), character())
    }
  )
})
