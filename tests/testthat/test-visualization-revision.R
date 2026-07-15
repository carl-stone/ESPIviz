visualization_plotly_points <- function(plot) {
  built <- plotly::plotly_build(plot)
  rows <- lapply(built$x$data, function(trace) {
    keys <- trace$key
    if (is.null(keys) || length(keys) == 0L) {
      return(NULL)
    }
    keys <- as.character(keys)
    data.frame(
      key = keys,
      x = as.numeric(trace$x)[seq_along(keys)],
      y = as.numeric(trace$y)[seq_along(keys)],
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      key = character(),
      x = numeric(),
      y = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

test_that("selection summaries preserve sample-level replicates", {
  expect_app_helper("summarize_selection")
  bundle <- synthetic_bundle()

  observed <- summarize_selection(bundle, character(), "Glul")
  by_sample <- observed$selected_by_sample

  expect_true(all(c(
    "sample",
    "condition",
    "gene",
    "cell_count",
    "mean_expression",
    "detected_pct"
  ) %in% names(by_sample)))
  expect_equal(nrow(by_sample), 4L)
  expect_setequal(
    as.character(by_sample$sample),
    unique(as.character(bundle$cells$sample))
  )
  expect_identical(
    stats::setNames(
      as.character(by_sample$condition),
      as.character(by_sample$sample)
    )[c("10_control", "33_control", "3_estim", "30_estim")],
    c(
      `10_control` = "p27CKO",
      `33_control` = "p27CKO",
      `3_estim` = "p27CKO +EStim",
      `30_estim` = "p27CKO +EStim"
    )
  )

  control <- by_sample[
    as.character(by_sample$sample) == "10_control",
    ,
    drop = FALSE
  ]
  expect_equal(control$cell_count, 2L)
  expected_expression <- expression_matrix(bundle, "Glul")[, 1L]
  expect_equal(control$mean_expression, mean(expected_expression[1:2]))
  expect_equal(control$detected_pct, 50)

  estim <- by_sample[
    as.character(by_sample$sample) == "3_estim",
    ,
    drop = FALSE
  ]
  expect_equal(estim$cell_count, 2L)
  expect_equal(estim$mean_expression, mean(expected_expression[4:5]))
  expect_equal(estim$detected_pct, 50)
})

test_that("marker overview completes selected markers across every cluster", {
  expect_app_helper("prepare_marker_overview")
  expect_app_helper("make_marker_overview_plot")
  bundle <- synthetic_bundle()

  observed <- prepare_marker_overview(bundle, "0", top_n = 2L)

  expect_true(all(c(
    "gene",
    "cluster",
    "mean_expression",
    "detected_pct",
    "within_gene_mean"
  ) %in% names(observed)))
  expect_setequal(observed$gene, c("Glul", "EGFP"))
  expect_setequal(as.character(observed$cluster), c("0", "1", "2"))
  expect_equal(nrow(observed), 6L)
  expect_identical(
    anyDuplicated(paste(observed$gene, observed$cluster, sep = "\r")),
    0L
  )

  glul_cluster_zero <- observed[
    observed$gene == "Glul" & as.character(observed$cluster) == "0",
    ,
    drop = FALSE
  ]
  glul_expression <- expression_matrix(bundle, "Glul")[, 1L]
  expect_equal(glul_cluster_zero$mean_expression, mean(glul_expression[1:2]))
  expect_equal(glul_cluster_zero$detected_pct, 50)

  egfp_cluster_two <- observed[
    observed$gene == "EGFP" & as.character(observed$cluster) == "2",
    ,
    drop = FALSE
  ]
  egfp_expression <- expression_matrix(bundle, "EGFP")[, 1L]
  expect_equal(egfp_cluster_two$mean_expression, mean(egfp_expression[5:6]))
  expect_equal(egfp_cluster_two$detected_pct, 50)
  expect_true(all(is.finite(observed$within_gene_mean)))
  expect_true(all(observed$within_gene_mean >= 0))
  expect_true(all(observed$within_gene_mean <= 1))
  expect_equal(
    range(observed$within_gene_mean[observed$gene == "Glul"]),
    c(0, 1)
  )

  plot <- make_marker_overview_plot(observed)
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("cluster composition keeps zero-count sample-cluster combinations", {
  expect_app_helper("prepare_cluster_composition")
  expect_app_helper("make_cluster_composition_plot")
  bundle <- synthetic_bundle()

  observed <- prepare_cluster_composition(bundle)
  sample_count <- length(unique(as.character(bundle$cells$sample)))
  cluster_count <- length(unique(as.character(bundle$cells$cluster)))

  expect_true(all(c(
    "sample",
    "condition",
    "cluster",
    "cell_count",
    "sample_total",
    "cell_pct"
  ) %in% names(observed)))
  expect_equal(nrow(observed), sample_count * cluster_count)
  expect_identical(
    anyDuplicated(paste(observed$sample, observed$cluster, sep = "\r")),
    0L
  )
  expect_equal(sum(observed$cell_count), nrow(bundle$cells))
  expect_equal(
    as.numeric(tapply(observed$cell_pct, observed$sample, sum)),
    rep(100, sample_count)
  )

  split_estim <- observed[
    as.character(observed$sample) == "3_estim" &
      as.character(observed$cluster) == "1",
    ,
    drop = FALSE
  ]
  expect_identical(as.character(split_estim$condition), "p27CKO +EStim")
  expect_equal(split_estim$cell_count, 1L)
  expect_equal(split_estim$sample_total, 2L)
  expect_equal(split_estim$cell_pct, 50)

  absent_cluster <- observed[
    as.character(observed$sample) == "10_control" &
      as.character(observed$cluster) == "2",
    ,
    drop = FALSE
  ]
  expect_equal(absent_cluster$cell_count, 0L)
  expect_equal(absent_cluster$sample_total, 2L)
  expect_equal(absent_cluster$cell_pct, 0)

  plot <- make_cluster_composition_plot(observed, bundle)
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("differential-expression data groups missing FDR with non-significant genes", {
  expect_app_helper("prepare_de_plot_data")
  primary_de <- synthetic_bundle()$primary_de
  primary_de$padj[primary_de$gene == "Other"] <- NA_real_

  expect_identical(
    unname(de_significance_symbols()),
    c("square", "circle", "triangle-up")
  )
  expect_identical(
    unname(de_significance_symbols("volcano")),
    c("square", "circle", "x", "triangle-up")
  )

  observed <- prepare_de_plot_data(primary_de, fdr_threshold = 0.05)

  expect_equal(observed$base_mean, c(100, 20, 75, 500))
  expect_equal(observed$probability[1:3], c(4e-5, 0.5, 4e-4))
  expect_equal(
    observed$minus_log10_probability[1:3],
    -log10(c(4e-5, 0.5, 4e-4))
  )
  expect_true(is.na(observed$probability[[4L]]))
  expect_true(is.na(observed$minus_log10_probability[[4L]]))
  expect_identical(
    levels(observed$significance),
    c(
      "Higher in Control",
      "Not significant",
      "Higher with E-Stim"
    )
  )
  expect_identical(
    as.character(observed$significance),
    c(
      "Higher with E-Stim",
      "Not significant",
      "Higher in Control",
      "Not significant"
    )
  )
  expect_identical(
    as.character(de_plot_status(observed, "volcano")),
    c(
      "Higher with E-Stim",
      "Not significant",
      "Higher in Control",
      "Adjusted P unavailable"
    )
  )

  stricter <- prepare_de_plot_data(primary_de, fdr_threshold = 1e-4)
  expect_identical(
    as.character(stricter$significance),
    c(
      "Higher with E-Stim",
      "Not significant",
      "Not significant",
      "Not significant"
    )
  )
})

test_that("differential-expression table projects readable result columns", {
  expect_app_helper("prepare_de_table_data")
  primary_de <- synthetic_bundle()$primary_de

  observed <- prepare_de_table_data(primary_de)

  expect_identical(
    names(observed),
    c(
      "Gene",
      "Base Mean",
      "log2FC",
      "LFC SE",
      "Adjusted P",
      "Mean Count Control",
      "Mean Count E-Stim"
    )
  )
  expect_identical(observed$Gene, primary_de$gene)
  expect_identical(observed$`Base Mean`, primary_de$baseMean)
  expect_identical(observed$log2FC, primary_de$log2FoldChange)
  expect_identical(observed$`LFC SE`, primary_de$lfcSE)
  expect_identical(observed$`Adjusted P`, primary_de$padj)
  expect_identical(
    observed$`Mean Count Control`,
    primary_de$mean_count_control
  )
  expect_identical(
    observed$`Mean Count E-Stim`,
    primary_de$mean_count_estim
  )
  expect_true(all(c("stat", "pvalue", "design") %in% names(primary_de)))
  expect_false(any(c("stat", "pvalue", "design") %in% names(observed)))
})

test_that("MA and volcano plots retain one clickable key per gene", {
  expect_app_helper("make_de_plotly")
  expect_app_helper("make_de_ggplot")
  primary_de <- synthetic_bundle()$primary_de
  primary_de$padj[primary_de$gene == "Other"] <- NA_real_

  ma <- make_de_plotly(
    primary_de,
    active_gene = "Glul",
    source = "de_ma_test",
    view = "ma"
  )
  volcano <- make_de_plotly(
    primary_de,
    active_gene = "Glul",
    source = "de_volcano_test",
    view = "volcano"
  )
  ma_points <- visualization_plotly_points(ma)
  volcano_points <- visualization_plotly_points(volcano)

  for (points in list(ma_points, volcano_points)) {
    expect_equal(nrow(points), nrow(primary_de))
    expect_identical(anyDuplicated(points$key), 0L)
    expect_setequal(points$key, primary_de$gene)
  }
  ma_built <- plotly::plotly_build(ma)
  expect_identical(ma_built$x$layout$legend$title$text, "MA category")
  ma_symbols <- unlist(lapply(ma_built$x$data, function(trace) {
    trace$marker$symbol %||% character()
  }), use.names = FALSE)
  expect_true(all(c("square", "circle", "triangle-up") %in% ma_symbols))
  expect_false("x" %in% ma_symbols)
  volcano_built <- plotly::plotly_build(volcano)
  expect_identical(volcano_built$x$layout$legend$title$text, "FDR status")
  volcano_symbols <- unlist(lapply(volcano_built$x$data, function(trace) {
    trace$marker$symbol %||% character()
  }), use.names = FALSE)
  expect_true("x" %in% volcano_symbols)

  ma_mcm2 <- ma_points[ma_points$key == "Mcm2", , drop = FALSE]
  expect_equal(ma_mcm2$x, 75 + 1)
  expect_equal(ma_mcm2$y, -1.5)
  ma_xaxis <- plotly::plotly_build(ma)$x$layout$xaxis
  expect_identical(ma_xaxis$type, "log")
  expect_identical(ma_xaxis$tickmode, "array")
  expect_true(all(log10(ma_xaxis$tickvals) %% 1 == 0))
  volcano_mcm2 <- volcano_points[
    volcano_points$key == "Mcm2",
    ,
    drop = FALSE
  ]
  expect_equal(volcano_mcm2$x, -1.5)
  expect_equal(volcano_mcm2$y, -log10(4e-4))
  expect_equal(
    volcano_points$y[volcano_points$key == "Other"],
    0
  )

  for (view in c("ma", "volcano")) {
    plot <- make_de_ggplot(primary_de, "Glul", view = view)
    expect_s3_class(plot, "ggplot")
    built <- ggplot2::ggplot_build(plot)
    expect_silent(built)
    full_point_layer <- Filter(
      function(layer) nrow(layer) == nrow(primary_de) && "shape" %in% names(layer),
      built$data
    )[[1L]]
    expected_shapes <- if (identical(view, "ma")) {
      c(15, 16, 17)
    } else {
      c(4, 15, 16, 17)
    }
    expect_setequal(unique(full_point_layer$shape), expected_shapes)
  }
})

test_that("pathway plots keep GSEA and ORA on method-specific axes", {
  expect_app_helper("prepare_pathway_plot_data")
  expect_app_helper("make_pathway_plotly")
  expect_app_helper("make_pathway_ggplot")
  pathways <- synthetic_bundle()$pathways
  pathways$source[[2L]] <- "ORA"
  pathways$direction[[1L]] <- "Control"
  pathways$direction[[2L]] <- "E-Stim"
  pathways$score[[2L]] <- 2.3

  observed <- prepare_pathway_plot_data(pathways)

  expect_identical(as.character(observed$source), c("GSEA", "ORA"))
  expect_identical(levels(observed$source), c("GSEA", "ORA"))
  expect_identical(observed$score_label, c("NES", "Fold enrichment"))
  expect_equal(observed$baseline, c(0, 1))
  expect_identical(
    observed$count_label,
    c("Source gene-set size", "Overlapping genes")
  )
  expect_equal(observed$evidence_strength, -log10(pathways$p_adjust))
  expect_identical(observed$source_order, c(1L, 2L))
  expect_identical(
    levels(observed$panel_label),
    c("GSEA · NES", "ORA · Fold enrichment")
  )
  expect_identical(observed$pathway_key, pathways$pathway_id)
  expect_identical(anyDuplicated(observed$pathway_key), 0L)
  expect_true(all(grepl("Method", observed$hover_text, fixed = TRUE)))
  expect_true(all(grepl("Adjusted P", observed$hover_text, fixed = TRUE)))
  expect_true(grepl("Source gene-set size", observed$hover_text[[1L]], fixed = TRUE))
  expect_true(grepl("Overlapping genes", observed$hover_text[[2L]], fixed = TRUE))

  interactive <- make_pathway_plotly(
    pathways,
    active_pathway = "oxidative_stress",
    source = "pathway_axis_test"
  )
  interactive_built <- plotly::plotly_build(interactive)
  keys <- unlist(lapply(interactive_built$x$data, function(trace) {
    as.character(trace$key %||% character())
  }), use.names = FALSE)
  expect_equal(length(keys), nrow(pathways))
  expect_identical(anyDuplicated(keys), 0L)
  expect_setequal(keys, pathways$pathway_id)
  pathway_symbols <- unlist(lapply(interactive_built$x$data, function(trace) {
    trace$marker$symbol %||% character()
  }), use.names = FALSE)
  expect_setequal(pathway_symbols, c("circle", "diamond"))
  x_axes <- grep(
    "^xaxis[0-9]*$",
    names(interactive_built$x$layout),
    value = TRUE
  )
  expect_gte(length(x_axes), 2L)

  static <- make_pathway_ggplot(pathways, "oxidative_stress")
  expect_s3_class(static, "ggplot")
  static_built <- ggplot2::ggplot_build(static)
  expect_equal(length(static_built$layout$panel_scales_x), 2L)
  expect_identical(anyDuplicated(static$data$pathway_key), 0L)
  expect_setequal(unique(static_built$data[[2L]]$shape), c(21, 22))

  duplicate <- rbind(pathways, pathways[1L, , drop = FALSE])
  expect_error(
    prepare_pathway_plot_data(duplicate),
    "must be unique",
    fixed = TRUE
  )

  expect_error(
    prepare_pathway_plot_data(pathways[0, , drop = FALSE]),
    "At least one",
    fixed = TRUE
  )

  duplicate_label <- pathways
  duplicate_label$label[[2L]] <- duplicate_label$label[[1L]]
  duplicate_label$source[[2L]] <- duplicate_label$source[[1L]]
  observed_duplicate_label <- prepare_pathway_plot_data(duplicate_label)
  expect_identical(anyDuplicated(observed_duplicate_label$plot_label), 0L)
  expect_true(all(mapply(
    grepl,
    pattern = as.character(observed_duplicate_label$direction),
    x = observed_duplicate_label$plot_label,
    MoreArgs = list(fixed = TRUE)
  )))
})

test_that("pathway plots use the top ten results per method and direction", {
  pathways <- do.call(rbind, lapply(c("GSEA", "ORA"), function(source) {
    do.call(rbind, lapply(c("Control", "E-Stim"), function(direction) {
      index <- seq_len(15L)
      data.frame(
        pathway_id = paste(tolower(source), direction, index, sep = "_"),
        label = paste(source, direction, index),
        source = source,
        direction = direction,
        description = paste(source, direction, index),
        p_value = index / 100,
        p_adjust = index / 20,
        score = if (identical(source, "GSEA")) {
          if (identical(direction, "E-Stim")) index else -index
        } else {
          index + 1
        },
        gene_count = index,
        stringsAsFactors = FALSE
      )
    }))
  }))

  observed <- top_pathway_results(pathways)
  groups <- interaction(observed$source, observed$direction, drop = TRUE)

  expect_equal(nrow(observed), 40L)
  expect_true(all(table(groups) == 10L))
  expect_true(all(observed$p_adjust <= 0.5))
  expect_true(any(observed$p_adjust >= 0.05))
})
