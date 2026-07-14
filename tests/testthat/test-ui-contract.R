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

test_that("Explore limits UMAP colors and exposes every selection mode", {
  expect_app_helper("explore_ui")
  html <- htmltools::renderTags(explore_ui(
    "explore_test",
    synthetic_bundle()
  ))$html

  expect_match(html, "Gene expression", fixed = TRUE)
  expect_match(html, ">Cluster<", fixed = TRUE)
  expect_match(html, ">Condition<", fixed = TRUE)
  expect_match(html, "Clear selection", fixed = TRUE)
  expect_match(html, "Click a cell or use box or lasso", fixed = TRUE)
  expect_no_match(html, "QC", fixed = TRUE)
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
  expect_match(html, "Mean expression", fixed = TRUE)
  expect_match(html, "dot size", fixed = TRUE)
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
  row <- bundle$pathways[1L, , drop = FALSE]
  genes <- bundle$pathway_genes$gene[
    bundle$pathway_genes$pathway_id == row$pathway_id[[1L]]
  ]
  html <- htmltools::renderTags(pathway_detail_ui(row, genes))$html

  expect_match(html, "<dl>", fixed = TRUE)
  expect_match(html, "<dt>Method</dt>", fixed = TRUE)
  expect_match(html, "<dd>", fixed = TRUE)
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
