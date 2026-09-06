test_that("blend strength respects detection even for negative and constant PFlog", {
  expect_equal(
    scale_expression_strength(c(-4, -3, -2), c(FALSE, TRUE, TRUE)),
    c(0, 0, 1)
  )
  expect_equal(
    scale_expression_strength(c(-4, -3, -2), rep(FALSE, 3)),
    c(0, 0, 0)
  )
  expect_equal(
    scale_expression_strength(c(-3, -3, 20), c(TRUE, TRUE, FALSE)),
    c(.5, .5, 0)
  )
  bundle <- synthetic_bundle()
  data <- prepare_plot_gene_data(bundle, c("Glul", "ZeroGene"))
  blended <- prepare_umap_expression_blend_data(data)
  expect_true(all(blended$strength_2 == 0))
  expect_true(all(
    blended$blend_color[!blended$detected_1 & !blended$detected_2] ==
      unname(blend_palette()[["Neither detected"]])
  ))
})

test_that("analysis and saved scopes stay separate and exports include both plot genes", {
  bundle <- synthetic_bundle(extra_genes = 35L)
  state <- new_app_state(bundle)
  shiny::testServer(
    explore_server,
    args = list(bundle = bundle, state = state),
    {
      session$setInputs(active_gene = "Glul", secondary_gene = "EGFP")
      expect_identical(analysis_genes(), c("Glul", "EGFP"))
      path <- output$download_summary
      expect_setequal(utils::read.csv(path)$gene, c("Glul", "EGFP"))
      path <- output$download_expression
      expect_setequal(readRDS(path)$genes$gene, c("Glul", "EGFP"))
      saved <- bundle$genes$gene[3:33]
      state$gene_set(saved)
      session$flushReact()
      expect_identical(
        analysis_genes(),
        analysis_gene_scope("Glul", "EGFP", saved)
      )
      expect_setequal(
        utils::read.csv(output$download_summary)$gene,
        analysis_genes()
      )
      expect_identical(readLines(output$download_gene_set), saved)
      expect_equal(length(page_info()$genes), 25L)
      session$setInputs(cluster_plot_type = "violin")
      expect_equal(length(page_info()$genes), 6L)
      session$setInputs(explore_results = "Cell-level")
      expect_identical(output$gene_page_count, "0")
    }
  )
})

test_that("saved view preserves every exploration setting and validates data identity", {
  bundle <- synthetic_bundle()
  state <- new_app_state(bundle)
  shiny::isolate({
    state$gene_set(c("EGFP", "Mcm2"))
    state$active_gene("Glul")
    state$gene_set_name("Example")
    state$selected_cells(bundle$cells$cell_id[c(2L, 4L)])
    options <- view_option_defaults()
    options$secondary_gene <- "EGFP"
    options$color_by <- "cluster"
    options$explore_results <- "By sample"
    options$sample_plot_type <- "dot"
    options$gene_page <- 2L
    state$view_options(options)
    captured <- capture_view_state(bundle, state)
    expect_null(view_link_query(captured))
    encoded <- jsonlite::toJSON(
      captured,
      auto_unbox = TRUE,
      null = "null",
      na = "null"
    )
    decoded <- jsonlite::fromJSON(encoded)
    target <- new_app_state(bundle)
    restore_view_state(decoded, bundle, target)
    expect_identical(target$active_gene(), "Glul")
    expect_identical(target$gene_set(), c("EGFP", "Mcm2"))
    expect_identical(target$selected_cells(), state$selected_cells())
    expect_equal(target$view_options(), options)
    captured$bundle$data_version <- "different"
    expect_error(validate_view_state(captured, bundle), "different data bundle")
    state$selected_cells(cell_ids_for_cluster(bundle, "1"))
    captured <- capture_view_state(bundle, state)
    expect_match(view_link_query(captured), "view_state=")
  })
})

test_that("replicate summaries retain absent samples without fabricated expression", {
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[1L]
  summary <- summarize_selection(
    bundle,
    selected,
    c("Glul", "EGFP")
  )$selected_by_sample
  expect_setequal(
    as.character(summary$sample),
    unique(as.character(bundle$cells$sample))
  )
  missing <- summary$cell_count == 0L
  expect_true(any(missing))
  expect_true(all(is.na(summary$mean_expression[missing])))
  expect_true(all(is.na(summary$detected_pct[missing])))
  expect_true(all(summary$detected_n[missing] == 0L))
  plot <- make_summary_violin_plot(
    prepare_summary_violin_data(bundle, "Glul", selected),
    bundle,
    "sample",
    "Sample"
  )
  built <- ggplot2::ggplot_build(plot)
  labels <- built$layout$panel_params[[1L]]$y$get_labels()
  expect_equal(length(labels), length(unique(bundle$cells$sample)))
  expect_true(any(grepl("0 selected cells", labels, fixed = TRUE)))
})

test_that("scientific result filters rank complete data and preserve missing probabilities", {
  data <- synthetic_bundle()$primary_de
  data$padj[4L] <- NA_real_
  all <- filter_result_rows(data)
  expect_equal(tail(all$padj, 1L), NA_real_)
  expect_equal(nrow(filter_result_rows(data, fdr = "missing")), 1L)
  expect_true(all(filter_result_rows(data, fdr = "significant")$padj <= .05))
  expect_equal(nrow(filter_result_rows(data, search = "no-such-gene")), 0L)
  pathways <- synthetic_bundle()$pathways
  expect_true(all(
    filter_result_rows(pathways, method = "GSEA")$source == "GSEA"
  ))
  expect_false("description" %in% names(pathway_table_data(pathways)))
})

test_that("current-gene cards report FDR independently from effect direction", {
  data <- synthetic_bundle()$primary_de
  data$log2FoldChange[1L] <- 1
  for (p in c(NA_real_, 0.2, 0.01)) {
    data$padj[1L] <- p
    card <- as.character(current_de_gene_ui(data, data$gene[1L]))
    expect_match(card, "Higher in p27CKO + E-Stim", fixed = TRUE)
    status <- if (is.na(p)) {
      "Adjusted P unavailable"
    } else if (p > 0.05) {
      "FDR &gt; 0.05"
    } else {
      "FDR ≤ 0.05"
    }
    expect_match(card, status, fixed = TRUE)
  }
})

test_that("zero detection has zero dot area and shared limits remain explicit", {
  bundle <- synthetic_bundle()
  summary <- summarize_selection(
    bundle,
    character(),
    c("Glul", "ZeroGene")
  )$selected_by_cluster
  plot <- make_group_summary_plot(
    summary,
    "cluster",
    "Cluster",
    color_limit = 8
  )
  built <- ggplot2::ggplot_build(plot)
  expect_true(all(built$data[[1L]]$size[summary$detected_pct == 0] == 0))
  expect_equal(plot$scales$get_scales("colour")$limits, c(-8, 8))
})
