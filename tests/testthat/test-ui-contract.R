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

test_that("Explore limits UMAP colors and exposes every selection mode", {
  expect_app_helper("explore_ui")
  html <- htmltools::renderTags(explore_ui("explore_test", synthetic_bundle()))$html

  expect_match(html, "Gene expression", fixed = TRUE)
  expect_match(html, ">Cluster<", fixed = TRUE)
  expect_match(html, ">Condition<", fixed = TRUE)
  expect_match(html, "Clear selection", fixed = TRUE)
  expect_match(html, "Click a cell or use box or lasso", fixed = TRUE)
  expect_no_match(html, "QC", fixed = TRUE)
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
