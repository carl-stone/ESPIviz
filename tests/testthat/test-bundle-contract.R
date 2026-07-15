test_that("synthetic bundle satisfies the public schema", {
  expect_app_helper("validate_bundle")
  bundle <- synthetic_bundle()

  expect_silent(validate_bundle(bundle, expected = synthetic_expectations(bundle)))
  expect_identical(rownames(bundle$counts), bundle$cells$cell_id)
  expect_identical(colnames(bundle$counts), bundle$genes$gene)
  expect_true(inherits(bundle$counts, "dgCMatrix"))
})

test_that("bundle validation rejects dimension and cell-order drift", {
  expect_app_helper("validate_bundle")
  bundle <- synthetic_bundle()
  expected <- synthetic_expectations(bundle)

  wrong_dimensions <- bundle
  wrong_dimensions$counts <- wrong_dimensions$counts[, -1, drop = FALSE]
  expect_error(
    validate_bundle(wrong_dimensions, expected = expected),
    "dimension|gene|count",
    ignore.case = TRUE
  )

  wrong_order <- bundle
  wrong_order$cells <- wrong_order$cells[rev(seq_len(nrow(wrong_order$cells))), ]
  expect_error(
    validate_bundle(wrong_order, expected = expected),
    "order|cell",
    ignore.case = TRUE
  )
})

test_that("bundle validation rejects invalid raw counts and library sizes", {
  expect_app_helper("validate_bundle")
  bundle <- synthetic_bundle()
  expected <- synthetic_expectations(bundle)

  fractional_counts <- bundle
  fractional_counts$counts@x[[1L]] <- 0.5
  expect_error(
    validate_bundle(fractional_counts, expected = expected),
    "raw counts",
    fixed = TRUE
  )

  negative_counts <- bundle
  negative_counts$counts@x[[1L]] <- -1
  expect_error(
    validate_bundle(negative_counts, expected = expected),
    "raw counts",
    fixed = TRUE
  )

  nonfinite_counts <- bundle
  nonfinite_counts$counts@x[[1L]] <- Inf
  expect_error(
    validate_bundle(nonfinite_counts, expected = expected),
    "raw counts",
    fixed = TRUE
  )

  mismatched_library_size <- bundle
  mismatched_library_size$cells$library_size[[1L]] <-
    mismatched_library_size$cells$library_size[[1L]] + 1
  expect_error(
    validate_bundle(mismatched_library_size, expected = expected),
    "raw-count row sums",
    fixed = TRUE
  )
})

test_that("bundle validation rejects schema and primary-result drift", {
  expect_app_helper("validate_bundle")
  bundle <- synthetic_bundle()
  expected <- synthetic_expectations(bundle)

  wrong_schema <- bundle
  wrong_schema$schema_version <- "999.0.0"
  expect_error(
    validate_bundle(wrong_schema, expected = expected),
    "schema",
    ignore.case = TRUE
  )

  wrong_de <- bundle
  wrong_de$primary_de <- wrong_de$primary_de[-1, , drop = FALSE]
  expect_error(
    validate_bundle(wrong_de, expected = expected),
    "24,601|24601|differential|DE|row",
    ignore.case = TRUE
  )

  missing_base_mean <- bundle
  missing_base_mean$primary_de$baseMean <- NULL
  expect_error(
    validate_bundle(missing_base_mean, expected = expected),
    "differential-expression",
    fixed = TRUE
  )

  for (column in c("lfcSE", "mean_count_control", "mean_count_estim")) {
    missing_display_column <- bundle
    missing_display_column$primary_de[[column]] <- NULL
    expect_error(
      validate_bundle(missing_display_column, expected = expected),
      "differential-expression",
      fixed = TRUE,
      info = paste("Missing required DE display column:", column)
    )
  }

  invalid_base_mean <- bundle
  invalid_base_mean$primary_de$baseMean[[1L]] <- -1
  expect_error(
    validate_bundle(invalid_base_mean, expected = expected),
    "differential-expression",
    fixed = TRUE
  )

  invalid_log2fc <- bundle
  invalid_log2fc$primary_de$log2FoldChange[[1L]] <- Inf
  expect_error(
    validate_bundle(invalid_log2fc, expected = expected),
    "differential-expression",
    fixed = TRUE
  )

  invalid_lfcse <- bundle
  invalid_lfcse$primary_de$lfcSE <-
    as.character(invalid_lfcse$primary_de$lfcSE)
  expect_error(
    validate_bundle(invalid_lfcse, expected = expected),
    "differential-expression",
    fixed = TRUE
  )

  invalid_padj <- bundle
  invalid_padj$primary_de$padj[[1L]] <- -1
  expect_error(
    validate_bundle(invalid_padj, expected = expected),
    "differential-expression",
    fixed = TRUE
  )

  invalid_pvalue <- bundle
  invalid_pvalue$primary_de$pvalue[[1L]] <- 2
  expect_error(
    validate_bundle(invalid_pvalue, expected = expected),
    "differential-expression",
    fixed = TRUE
  )

  invalid_mean_count <- bundle
  invalid_mean_count$primary_de$mean_count_control[[1L]] <- -1
  expect_error(
    validate_bundle(invalid_mean_count, expected = expected),
    "differential-expression",
    fixed = TRUE
  )
})

test_that("bundle validation accepts repeated labels but requires unique rows", {
  expect_app_helper("validate_bundle")
  bundle <- synthetic_bundle()
  expected <- synthetic_expectations(bundle)

  empty_pathways <- bundle
  empty_pathways$pathways <- empty_pathways$pathways[0, , drop = FALSE]
  empty_pathways$pathway_genes <- empty_pathways$pathway_genes[0, , drop = FALSE]
  expect_error(
    validate_bundle(empty_pathways, expected = expected),
    "enrichment result tables",
    fixed = TRUE
  )

  duplicate_pathways <- bundle
  duplicate_pathways$pathways$pathway_id[[2L]] <-
    duplicate_pathways$pathways$pathway_id[[1L]]
  expect_error(
    validate_bundle(duplicate_pathways, expected = expected),
    "enrichment result tables",
    fixed = TRUE
  )

  duplicate_labels <- bundle
  duplicate_labels$pathways$label[[2L]] <-
    duplicate_labels$pathways$label[[1L]]
  expect_silent(validate_bundle(duplicate_labels, expected = expected))
})

test_that("bundle loading verifies the asset checksum", {
  expect_app_helper("load_bundle")
  bundle <- synthetic_bundle()
  path <- tempfile(fileext = ".rds")
  saveRDS(bundle, path, version = 3)
  hash <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  manifest <- list(sha256 = hash, schema_version = "1.0.0")

  expect_silent(
    load_bundle(path, manifest, expected = synthetic_expectations(bundle))
  )

  manifest$sha256 <- paste(rep("0", 64), collapse = "")
  expect_error(
    load_bundle(path, manifest, expected = synthetic_expectations(bundle)),
    "checksum|sha",
    ignore.case = TRUE
  )
})

test_that("bundle surface excludes private and analysis-development fields", {
  bundle <- synthetic_bundle()
  public_names <- c(
    names(bundle),
    names(bundle$cells),
    names(bundle$provenance)
  )
  forbidden <- c(
    "barcode", "nFeature", "nCount", "percent.mt", "qc", "paired",
    "sensitivity", "pflog", "seurat", "cluster_abundance", "box_path"
  )

  for (term in forbidden) {
    expect_false(
      any(grepl(tolower(term), tolower(public_names), fixed = TRUE)),
      info = paste("Forbidden public field:", term)
    )
  }
})
