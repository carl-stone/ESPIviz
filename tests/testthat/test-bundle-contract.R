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
