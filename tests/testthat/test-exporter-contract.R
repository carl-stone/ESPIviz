make_source_manifest <- function(source_sha256) {
  names <- c(
    "source_object",
    "primary_de",
    "markers",
    "ora_up",
    "ora_down",
    "gsea",
    "featured_genes_config",
    "featured_pathways_config"
  )
  directory <- tempfile("espiviz-source-")
  dir.create(directory)
  entries <- lapply(names, function(name) {
    path <- file.path(directory, paste0(name, ".txt"))
    writeLines(name, path)
    list(path = normalizePath(path, mustWork = TRUE))
  })
  names(entries) <- names
  entries$source_object$sha256 <- source_sha256
  c(list(manifest_version = "1.0.0"), entries)
}

write_source_manifest <- function(manifest) {
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE)
  path
}

make_synthetic_seurat_source <- function() {
  counts <- Matrix::Matrix(
    matrix(c(10, 0, 2, 0, 5, 1, 2, 1, 0, 3, 2, 8), nrow = 3L),
    sparse = TRUE
  )
  dimnames(counts) <- list(c("Glul", "EGFP", "Other"), paste0("barcode_", 1:4))
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  normalized <- log1p(
    sweep(as.matrix(counts), 2L, Matrix::colSums(counts), "/") * 10000
  )
  normalized <- Matrix::Matrix(normalized, sparse = TRUE)
  dimnames(normalized) <- dimnames(counts)
  SeuratObject::LayerData(object[["RNA"]], layer = "data") <- normalized

  embeddings <- matrix(seq_len(8L), ncol = 2L)
  rownames(embeddings) <- colnames(counts)
  colnames(embeddings) <- c("UMAP_1", "UMAP_2")
  object[["test_umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "UMAP_",
    assay = "RNA"
  )
  object$test_cluster <- c("1", "1", "2", "2")
  object$Condition <- c("p27CKO", "p27CKO", "p27CKO +EStim", "p27CKO +EStim")
  object$Mouse <- c("10", "10", "3", "3")
  object$sample_id <- c("10_control", "10_control", "3_estim", "3_estim")
  object
}

synthetic_source_contract <- function() {
  utils::modifyList(
    espiviz_contract(),
    list(
      reduction = "test_umap",
      cluster_column = "test_cluster",
      genes = 3L,
      cells = 4L,
      clusters = c("1", "2")
    )
  )
}

test_that("exporter pins the frozen structural source contract", {
  expect_true(exists("espiviz_contract", mode = "function"))
  contract <- espiviz_contract()

  expect_identical(
    contract$reduction,
    "umap_pflog_mg_selected_no_filter_cc_dims20"
  )
  expect_identical(
    contract$cluster_column,
    "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5"
  )
  expect_identical(contract$genes, 38394L)
  expect_identical(contract$cells, 3456L)
  expect_identical(contract$clusters, as.character(1:8))
  expect_identical(contract$primary_de_rows, 24601L)
  expect_identical(contract$primary_de_design, "primary_unpaired_condition")
})

test_that("source-manifest example requires one complete source hash", {
  manifest <- jsonlite::read_json(
    file.path(repo_root, "config", "source-manifest.example.json"),
    simplifyVector = TRUE
  )

  expect_match(manifest$source_object$sha256, "^[0-9a-f]{64}$")
})

test_that("exporter rejects source content that does not match its declared hash", {
  manifest <- make_source_manifest(paste(rep("0", 64), collapse = ""))
  path <- write_source_manifest(manifest)

  expect_error(
    build_espiviz_bundle(path),
    "source_object SHA-256",
    ignore.case = TRUE
  )
})

test_that("exporter rejects incomplete source hashes", {
  manifest <- make_source_manifest("not-a-hash")
  path <- write_source_manifest(manifest)

  expect_error(
    build_espiviz_bundle(path),
    "complete SHA-256",
    ignore.case = TRUE
  )
})

test_that("public-bundle validator accepts only the narrow data surface", {
  expect_true(exists("validate_public_bundle", mode = "function"))
  bundle <- synthetic_bundle()

  expect_silent(validate_public_bundle(bundle, enforce_frozen = FALSE))

  wrong_cell_order <- bundle
  wrong_cell_order$cells <- wrong_cell_order$cells[rev(seq_len(nrow(bundle$cells))), ]
  expect_error(
    validate_public_bundle(wrong_cell_order, enforce_frozen = FALSE),
    "row order|cell IDs",
    ignore.case = TRUE
  )

  with_qc <- bundle
  with_qc$cells$nFeature_RNA <- 1000L
  expect_error(
    validate_public_bundle(with_qc, enforce_frozen = FALSE),
    "unsupported|columns",
    ignore.case = TRUE
  )

  with_barcode <- bundle
  with_barcode$cells$barcode <- paste0("original_", seq_len(nrow(bundle$cells)))
  expect_error(
    validate_public_bundle(with_barcode, enforce_frozen = FALSE),
    "unsupported|columns",
    ignore.case = TRUE
  )

  with_private_path <- bundle
  with_private_path$provenance$input_path <- "/Users/example/Library/CloudStorage/input.rds"
  expect_error(
    validate_public_bundle(with_private_path, enforce_frozen = FALSE),
    "forbidden|path|Box",
    ignore.case = TRUE
  )
})

test_that("public-bundle validator rejects alternate DE surfaces and schemas", {
  bundle <- synthetic_bundle()

  paired <- bundle
  paired$primary_de$design <- "paired_sensitivity"
  expect_error(
    validate_public_bundle(paired, enforce_frozen = FALSE),
    "non-primary|paired",
    ignore.case = TRUE
  )

  extra_table <- bundle
  extra_table$paired_sensitivity <- bundle$primary_de
  expect_error(
    validate_public_bundle(extra_table, enforce_frozen = FALSE),
    "fields must be exactly",
    ignore.case = TRUE
  )

  missing_result_column <- bundle
  missing_result_column$primary_de$padj <- NULL
  expect_error(
    validate_public_bundle(missing_result_column, enforce_frozen = FALSE),
    "missing required columns",
    ignore.case = TRUE
  )
})

test_that("exporter normalization check verifies source data-layer identity", {
  counts <- methods::as(Matrix::Matrix(
    matrix(c(10, 0, 5, 20, 15, 0), nrow = 2L),
    sparse = TRUE
  ), "dgCMatrix")
  dimnames(counts) <- list(c("Glul", "EGFP"), c("cell_a", "cell_b", "cell_c"))
  library_size <- Matrix::colSums(counts)
  normalized <- log1p(sweep(as.matrix(counts), 2L, library_size, "/") * 10000)
  normalized <- methods::as(Matrix::Matrix(normalized, sparse = TRUE), "dgCMatrix")

  observed <- espiviz_validate_normalization(counts, normalized)
  expect_lte(observed$max_abs_error, 1e-12)

  normalized[1, 1] <- normalized[1, 1] + 0.01
  expect_error(
    espiviz_validate_normalization(counts, normalized),
    "does not match",
    ignore.case = TRUE
  )
})

test_that("source extraction enforces dimensions, reduction, and cluster column", {
  object <- make_synthetic_seurat_source()
  contract <- synthetic_source_contract()

  extracted <- espiviz_extract_source(object, contract)
  expect_equal(dim(extracted$counts), c(4L, 3L))
  expect_identical(extracted$cells$cell_id, paste0("cell_", 1:4))
  expect_identical(extracted$cells$sample, c(
    "10_control", "10_control", "3_estim", "3_estim"
  ))
  expect_false(any(grepl("barcode", extracted$cells$cell_id, fixed = TRUE)))
  expect_lte(extracted$normalization_check$max_abs_error, 1e-12)

  wrong_dimensions <- contract
  wrong_dimensions$genes <- 4L
  expect_error(
    espiviz_extract_source(object, wrong_dimensions),
    "dimensions",
    ignore.case = TRUE
  )

  wrong_reduction <- contract
  wrong_reduction$reduction <- "missing_umap"
  expect_error(
    espiviz_extract_source(object, wrong_reduction),
    "missing reduction",
    ignore.case = TRUE
  )

  wrong_cluster <- contract
  wrong_cluster$cluster_column <- "missing_cluster"
  expect_error(
    espiviz_extract_source(object, wrong_cluster),
    "missing required columns",
    ignore.case = TRUE
  )
})

test_that("source extraction rejects a normalized layer that drifts from counts", {
  object <- make_synthetic_seurat_source()
  normalized <- SeuratObject::LayerData(object[["RNA"]], layer = "data")
  normalized[1, 1] <- normalized[1, 1] + 0.01
  SeuratObject::LayerData(object[["RNA"]], layer = "data") <- normalized

  expect_error(
    espiviz_extract_source(object, synthetic_source_contract()),
    "does not match",
    ignore.case = TRUE
  )
})

test_that("build-data script can be sourced without executing the CLI", {
  build_script <- file.path(repo_root, "scripts", "build-data.R")
  environment <- new.env(parent = globalenv())

  expect_silent(source(build_script, local = environment))
  expect_true(exists("espiviz_run_cli", envir = environment, mode = "function"))
})
