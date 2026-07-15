test_that("application exposes exactly the four public navigation areas", {
  expect_app_helper("app_ui")
  html <- htmltools::renderTags(app_ui(synthetic_bundle()))$html

  for (label in c("Explore", "Differential expression", "Pathways", "About")) {
    expect_match(html, paste0(">", label, "<"), fixed = TRUE)
  }
  expect_no_match(
    html,
    ">(Descriptive|Replicate-aware|Secondary|Enriched|Depleted)<",
    perl = TRUE
  )
})

test_that("desktop pages use natural height so plots and tables do not collapse", {
  html <- htmltools::renderTags(app_ui(synthetic_bundle()))$html

  expect_no_match(html, '<body class="bslib-page-fill', fixed = TRUE)
})

test_that("Explore result panels span the full plotting grid", {
  html <- htmltools::renderTags(
    explore_ui("explore_test", synthetic_bundle())
  )$html

  expect_match(html, 'col-widths-sm="8,4,12"', fixed = TRUE)
})

test_that("secondary full-width panels span their plotting grids", {
  de_html <- htmltools::renderTags(differential_expression_ui("de_test"))$html
  pathway_html <- htmltools::renderTags(pathways_ui("pathways_test"))$html
  about_html <- htmltools::renderTags(about_ui("about_test"))$html

  expect_match(de_html, 'col-widths-sm="8,4,12"', fixed = TRUE)
  expect_match(pathway_html, 'col-widths-sm="7,5,12"', fixed = TRUE)
  expect_match(about_html, 'col-widths-sm="7,5,12"', fixed = TRUE)
})

test_that("differential-expression table has a meaningful result heading", {
  html <- htmltools::renderTags(
    differential_expression_ui("de_test")
  )$html

  expect_match(html, "Pseudobulk DE results", fixed = TRUE)
  expect_no_match(html, "Complete primary result", fixed = TRUE)
})

test_that("Explore limits UMAP colors and exposes every selection mode", {
  expect_app_helper("explore_ui")
  html <- htmltools::renderTags(explore_ui(
    "explore_test",
    synthetic_bundle()
  ))$html

  expect_match(html, "Log normalized expression", fixed = TRUE)
  expect_match(html, ">Detection<", fixed = TRUE)
  expect_match(html, ">Cluster<", fixed = TRUE)
  expect_match(html, ">Condition<", fixed = TRUE)
  expect_match(html, "Clear selection", fixed = TRUE)
  expect_match(html, "Click a cell or use box or lasso", fixed = TRUE)
  expect_no_match(html, "QC", fixed = TRUE)
})

test_that("Explore summary tables use plain-language expression labels", {
  table <- summary_datatable(data.frame(
    gene = "Glul",
    mean_expression = 1.25,
    detected_pct = 90,
    stringsAsFactors = FALSE
  ))

  expect_identical(
    names(table$x$data),
    c("Gene", "Mean log normalized expression", "Detected (%)")
  )
})

test_that("Explore exposes an explicit whole-cluster selection control", {
  html <- htmltools::renderTags(
    explore_ui("explore_test", synthetic_bundle())
  )$html

  expect_match(html, "Select cells by cluster", fixed = TRUE)
  expect_match(html, "explore_test-select_cluster", fixed = TRUE)
  expect_match(html, "explore_test-select_cluster_cells", fixed = TRUE)
})

test_that("Explore provides plotted summaries by cluster and condition", {
  html <- htmltools::renderTags(
    explore_ui("explore_test", synthetic_bundle())
  )$html

  expect_match(html, "explore_test-cluster_summary_plot_ui", fixed = TRUE)
  expect_match(html, "explore_test-condition_summary_plot_ui", fixed = TRUE)
  expect_match(html, "Mean log normalized expression", fixed = TRUE)
  expect_match(html, "dot size", fixed = TRUE)
})

test_that("expression summary plots default to violins with dot plots available", {
  html <- htmltools::renderTags(
    explore_ui("explore_test", synthetic_bundle())
  )$html

  for (id in c(
    "comparison_plot_type",
    "sample_plot_type",
    "condition_plot_type",
    "cluster_plot_type"
  )) {
    select_id <- paste0("explore_test-", id)
    pattern <- paste0(
      '(?s)<select[^>]*id="',
      select_id,
      '"[^>]*>.*?</select>'
    )
    matched <- regmatches(html, regexpr(pattern, html, perl = TRUE))

    expect_length(matched, 1L)
    expect_match(
      matched,
      '<option value="violin" selected>Violin plot</option>',
      fixed = TRUE
    )
    expect_match(
      matched,
      '<option value="dot">Dot plot</option>',
      fixed = TRUE
    )
  }
})

test_that("expression summary plot selectors switch every rendered plot", {
  bundle <- synthetic_bundle()
  state <- new_app_state(bundle)
  plot_ui_ids <- c(
    "gene_comparison_plot_ui",
    "sample_summary_plot_ui",
    "condition_summary_plot_ui",
    "cluster_summary_plot_ui"
  )
  plot_ids <- c(
    "gene_comparison_plot",
    "sample_summary_plot",
    "condition_summary_plot",
    "cluster_summary_plot"
  )

  shiny::testServer(
    explore_server,
    args = list(bundle = bundle, state = state),
    {
      session$flushReact()
      for (id in plot_ui_ids) {
        expect_match(
          htmltools::renderTags(output[[id]])$html,
          'height:360px;',
          fixed = TRUE
        )
      }
      for (id in plot_ids) {
        expect_silent(output[[id]])
      }

      session$setInputs(
        comparison_plot_type = "dot",
        sample_plot_type = "dot",
        condition_plot_type = "dot",
        cluster_plot_type = "dot"
      )
      session$flushReact()
      for (id in plot_ui_ids) {
        expect_match(
          htmltools::renderTags(output[[id]])$html,
          'height:300px;',
          fixed = TRUE
        )
      }
      for (id in plot_ids) {
        expect_silent(output[[id]])
      }
    }
  )
})

test_that("Explore exposes replicate-aware visualization panels", {
  html <- htmltools::renderTags(
    explore_ui("explore_test", synthetic_bundle())
  )$html

  for (label in c(
    "By sample",
    "Pooled condition",
    "By cluster",
    "Cluster markers"
  )) {
    expect_match(html, paste0(">", label, "<"), fixed = TRUE)
  }
  for (id in c(
    "sample_summary_plot_ui",
    "sample_table",
    "composition_plot",
    "composition_table",
    "marker_overview_plot_ui"
  )) {
    expect_match(html, paste0("explore_test-", id), fixed = TRUE)
  }
  expect_match(
    html,
    "Select a marker row below to make that gene current",
    fixed = TRUE
  )
  expect_no_match(html, "Explore chosen marker", fixed = TRUE)

  source_text <- paste(
    readLines(
      file.path(repo_root, "app", "R", "10-explore.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  for (server_rendered_id in c(
    "sample_summary_plot",
    "marker_overview_plot"
  )) {
    expect_match(
      source_text,
      paste0('"', server_rendered_id, '"'),
      fixed = TRUE
    )
  }
  plot_chunks <- strsplit(source_text, "\n    output$", fixed = TRUE)[[1L]]
  static_plot_ids <- c(
    "gene_comparison_plot",
    "violin_plot",
    "condition_summary_plot",
    "sample_summary_plot",
    "composition_plot",
    "cluster_summary_plot",
    "marker_overview_plot"
  )
  for (plot_id in static_plot_ids) {
    chunk <- plot_chunks[
      startsWith(plot_chunks, paste0(plot_id, " <- shiny::renderPlot"))
    ]
    expect_length(chunk, 1L)
    expect_match(chunk, "alt = function()", fixed = TRUE)
  }
})

test_that("all requested application downloads are wired into the UI", {
  module_files <- file.path(
    repo_root,
    "app",
    "R",
    c(
      "10-explore.R",
      "20-differential-expression.R",
      "30-pathways.R",
      "40-about.R"
    )
  )
  source_text <- paste(
    unlist(lapply(module_files, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  download_ids <- c(
    "download_umap_png",
    "download_umap_pdf",
    "download_gene_set",
    "download_summary",
    "download_expression",
    "download_metadata",
    "download_de",
    "download_volcano_png",
    "download_volcano_pdf",
    "download_pathway",
    "download_bundle"
  )

  for (id in download_ids) {
    expect_match(source_text, id, fixed = TRUE)
  }
})

test_that("the maintenance surface renders without application data", {
  expect_app_helper("maintenance_ui")
  html <- htmltools::renderTags(maintenance_ui())$html

  expect_match(html, "temporarily unavailable", fixed = TRUE)
  expect_match(html, "data could not be loaded", fixed = TRUE)
})

test_that("pathway details render with supported HTML tags", {
  expect_app_helper("pathway_detail_ui")
  bundle <- synthetic_bundle()
  gsea_row <- bundle$pathways[1L, , drop = FALSE]
  gsea_genes <- bundle$pathway_genes$gene[
    bundle$pathway_genes$pathway_id == gsea_row$pathway_id[[1L]]
  ]
  gsea_html <- htmltools::renderTags(
    pathway_detail_ui(gsea_row, gsea_genes)
  )$html
  ora_row <- bundle$pathways[2L, , drop = FALSE]
  ora_row$source <- "ORA"
  ora_genes <- bundle$pathway_genes$gene[
    bundle$pathway_genes$pathway_id == ora_row$pathway_id[[1L]]
  ]
  ora_html <- htmltools::renderTags(
    pathway_detail_ui(ora_row, ora_genes)
  )$html

  expect_match(gsea_html, "<dl>", fixed = TRUE)
  expect_match(gsea_html, "<dt>Method</dt>", fixed = TRUE)
  expect_match(gsea_html, "<dt>Ontology</dt>", fixed = TRUE)
  expect_match(
    gsea_html,
    "Gene Ontology (GO), Biological Process (BP)",
    fixed = TRUE
  )
  expect_match(gsea_html, "Source gene-set size", fixed = TRUE)
  expect_match(gsea_html, "Exported leading-edge genes", fixed = TRUE)
  expect_match(ora_html, "Overlapping genes", fixed = TRUE)
  expect_match(ora_html, "Exported overlapping genes", fixed = TRUE)
})

test_that("Pathways presents top plots and side-by-side enrichment tables", {
  html <- htmltools::renderTags(pathways_ui("pathways_test"))$html

  expect_match(html, "Top pathway results", fixed = TRUE)
  expect_match(html, ">Enrichment results<", fixed = TRUE)
  expect_match(html, ">Genes<", fixed = TRUE)
  expect_match(html, 'col-widths-sm="9,3"', fixed = TRUE)
  expect_lt(
    regexpr(">Enrichment results<", html, fixed = TRUE)[[1L]],
    regexpr(">Genes<", html, fixed = TRUE)[[1L]]
  )
  expect_no_match(html, "Featured", fixed = TRUE)
  expect_no_match(html, "featured", fixed = TRUE)
})

test_that("About data links render with supported list tags", {
  expect_app_helper("data_links_ui")
  bundle <- synthetic_bundle()
  attr(bundle, "data_manifest") <- list(
    asset_url = "https://example.org/espiviz-data-v1.0.0.rds",
    asset_sha256 = paste(rep("a", 64L), collapse = "")
  )
  html <- htmltools::renderTags(data_links_ui(bundle))$html

  expect_match(html, "<ul", fixed = TRUE)
  expect_match(html, "<li>", fixed = TRUE)
  expect_match(html, "Processed app bundle", fixed = TRUE)
  expect_match(html, "Application source code", fixed = TRUE)
})

test_that("deep-link view names map without falling back early", {
  expect_app_helper("app_view_label")
  expect_app_helper("app_view_slug")

  expect_identical(app_view_label("de"), "Differential expression")
  expect_identical(
    app_view_label("DIFFERENTIAL_EXPRESSION"),
    "Differential expression"
  )
  expect_identical(app_view_label("pathways"), "Pathways")
  expect_identical(app_view_label("about"), "About")
  expect_null(app_view_label("unknown"))
  expect_identical(app_view_slug("Differential expression"), "de")
  expect_identical(app_view_slug("Pathways"), "pathways")
  expect_identical(app_view_slug("About"), "about")
})

test_that("narrow desktop layouts collapse and overlay the Explore sidebar", {
  styles <- paste(
    readLines(file.path(repo_root, "app", "www", "styles.css"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(styles, "@media (max-width: 900px)", fixed = TRUE)
  expect_match(styles, "--bslib-sidebar-js-window-size: mobile", fixed = TRUE)
  expect_match(
    styles,
    ".bslib-sidebar-layout:not(.sidebar-collapsed) > .sidebar",
    fixed = TRUE
  )
})
