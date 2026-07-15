test_that("selected-cell expression export contains exact long-form values", {
  expect_app_helper("selection_expression_export")
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[c(2, 5)]
  genes <- c("Glul", "EGFP", "ZeroGene")

  exported <- selection_expression_export(bundle, selected, genes)

  expect_equal(nrow(exported), length(selected) * length(genes))
  expect_true(all(c(
    "cell_id", "gene", "count", "normalized_expression",
    "condition", "cluster", "Mouse", "sample"
  ) %in% names(exported)))
  expect_setequal(exported$cell_id, selected)
  expect_setequal(exported$gene, genes)

  target <- exported[
    exported$cell_id == selected[[1]] & exported$gene == "EGFP",
    ,
    drop = FALSE
  ]
  cell_index <- match(selected[[1]], bundle$cells$cell_id)
  expected_count <- as.numeric(bundle$counts[cell_index, "EGFP"])
  expect_equal(target$count, expected_count)
  expect_equal(
    target$normalized_expression,
    expression_matrix(bundle, "EGFP", selected[[1]])[1L, 1L]
  )
})

test_that("exports retain a complete gene set beyond the visual page", {
  expect_app_helper("selection_expression_export")
  expect_app_helper("selection_summary_export")
  bundle <- synthetic_bundle(extra_genes = 35L)
  genes <- bundle$genes$gene
  selected <- bundle$cells$cell_id[1:2]

  expression <- selection_expression_export(bundle, selected, genes)
  summaries <- selection_summary_export(bundle, selected, genes)

  expect_setequal(expression$gene, genes)
  expect_equal(nrow(expression), length(selected) * length(genes))
  expect_setequal(summaries$gene, genes)
})

test_that("large expression exports use scclrR sparse-plus-center storage", {
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[c(1, 4, 6)]
  genes <- c("Glul", "EGFP", "ZeroGene")
  path <- tempfile(fileext = ".rds")

  write_selection_expression_export(bundle, selected, genes, path)
  exported <- readRDS(path)
  normalized <- exported$normalized_expression

  expect_named(normalized, c("sparse", "center", "k", "alpha"))
  expect_s4_class(normalized$sparse, "dgCMatrix")
  expect_identical(dim(normalized$sparse), c(length(genes), length(selected)))
  expect_identical(rownames(normalized$sparse), genes)
  expect_identical(colnames(normalized$sparse), selected)
  expect_identical(names(normalized$center), selected)

  reconstructed <- sweep(
    as.matrix(normalized$sparse),
    2L,
    normalized$center,
    "-"
  )
  expected <- expression_matrix(bundle, genes, selected)
  attr(expected, "missing_genes") <- NULL
  expect_equal(t(reconstructed), expected, tolerance = 1e-12)
})

test_that("the public bundle URL comes from the immutable manifest", {
  expect_app_helper("public_bundle_url")
  expected_url <- paste0(
    "https://github.com/carl-stone/ESPIviz/releases/download/",
    "data-v1.0.0/espiviz-data-v1.0.0.rds"
  )
  manifest <- list(
    asset_url = expected_url,
    sha256 = paste(rep("a", 64), collapse = ""),
    schema_version = "1.0.0"
  )

  expect_identical(public_bundle_url(manifest), expected_url)
  expect_true(is.na(public_bundle_url(list())))
})

test_that("UMAP plots write valid PNG and PDF files", {
  expect_app_helper("make_umap_ggplot")
  bundle <- synthetic_bundle()
  plot <- make_umap_ggplot(
    bundle = bundle,
    color_by = "expression",
    gene = "Glul",
    selected_cell_ids = bundle$cells$cell_id[1:2]
  )
  png_path <- tempfile(fileext = ".png")
  pdf_path <- tempfile(fileext = ".pdf")

  ggplot2::ggsave(png_path, plot, width = 7, height = 5, units = "in", dpi = 120)
  ggplot2::ggsave(pdf_path, plot, width = 7, height = 5, units = "in")

  expect_gt(file.info(png_path)$size, 1000)
  expect_gt(file.info(pdf_path)$size, 1000)
})
