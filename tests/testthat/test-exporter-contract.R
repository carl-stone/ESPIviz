make_source_manifest <- function(source_sha256) {
  names <- c(
    "source_object",
    "primary_de",
    "markers",
    "ora_up",
    "ora_down",
    "gsea",
    "featured_genes_config"
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
  expect_identical(contract$data_version, "1.2.0")
})

test_that("source-manifest example requires one complete source hash", {
  manifest <- jsonlite::read_json(
    file.path(repo_root, "config", "source-manifest.example.json"),
    simplifyVector = TRUE
  )

  expect_match(manifest$source_object$sha256, "^[0-9a-f]{64}$")
  expect_setequal(
    names(manifest),
    c(
      "manifest_version",
      "source_object",
      "primary_de",
      "markers",
      "ora_up",
      "ora_down",
      "gsea",
      "featured_genes_config"
    )
  )
})

test_that("exporter includes every GSEA and ORA enrichment result", {
  directory <- tempfile("espiviz-pathways-")
  dir.create(directory)
  paths <- stats::setNames(
    file.path(directory, paste0(c("ora_up", "ora_down", "gsea"), ".tsv")),
    c("ora_up", "ora_down", "gsea")
  )
  gsea <- data.frame(
    ID = c("GO:0001", "GO:0002"),
    Description = c("Shared term", "Down term"),
    pvalue = c(0.01, 0.4),
    p.adjust = c(0.02, 0.8),
    NES = c(2, -1.5),
    setSize = c(20L, 30L),
    core_enrichment = c("1/2", "3"),
    check.names = FALSE
  )
  ora_up <- data.frame(
    ID = "GO:0001",
    Description = "Shared term",
    pvalue = 0.3,
    p.adjust = 0.6,
    FoldEnrichment = 3,
    Count = 2L,
    geneID = "Glul/Other",
    direction = "E-Stim",
    check.names = FALSE
  )
  ora_down <- data.frame(
    ID = "GO:0003",
    Description = "Control term",
    pvalue = 0.04,
    p.adjust = 0.08,
    FoldEnrichment = 2,
    Count = 1L,
    geneID = "Mcm2",
    direction = "Control",
    check.names = FALSE
  )
  utils::write.table(
    gsea,
    paths[["gsea"]],
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    ora_up,
    paths[["ora_up"]],
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    ora_down,
    paths[["ora_down"]],
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  observed <- espiviz_read_pathway_results(
    result_paths = paths,
    genes = c("Glul", "Other", "Mcm2"),
    map_entrez = function(ids) {
      unname(c("1" = "Glul", "2" = "Other", "3" = "Mcm2")[ids])
    },
    map_ontology = function(ids) {
      rep("BP", length(ids))
    }
  )

  expect_equal(nrow(observed$pathways), 4L)
  expect_identical(anyDuplicated(observed$pathways$pathway_id), 0L)
  expect_gt(anyDuplicated(observed$pathways$label), 0L)
  expect_setequal(observed$pathways$p_adjust, c(0.02, 0.8, 0.6, 0.08))
  expect_true(any(observed$pathways$p_adjust > 0.05))
  expect_setequal(
    observed$pathways$direction,
    c("E-Stim", "Control")
  )
  expect_setequal(observed$pathway_genes$gene, c("Glul", "Other", "Mcm2"))

  expect_error(
    espiviz_read_pathway_results(
      result_paths = paths,
      genes = c("Glul", "Other", "Mcm2"),
      map_entrez = function(ids) {
        mapped <- unname(c(
          "1" = "Glul",
          "2" = "Other",
          "3" = "Mcm2"
        )[ids])
        mapped[[1L]] <- NA_character_
        mapped
      },
      map_ontology = function(ids) rep("BP", length(ids))
    ),
    "every GSEA Entrez ID",
    fixed = TRUE
  )

  expect_error(
    espiviz_read_pathway_results(
      result_paths = paths,
      genes = c("Glul", "Other", "Mcm2"),
      map_entrez = function(ids) {
        unname(c("1" = "Glul", "2" = "Other", "3" = "Mcm2")[ids])
      },
      map_ontology = function(ids) {
        c("BP", rep("MF", length(ids) - 1L))
      }
    ),
    "Biological Process",
    fixed = TRUE
  )
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

  missing_base_mean <- bundle
  missing_base_mean$primary_de$baseMean <- NULL
  expect_error(
    validate_public_bundle(missing_base_mean, enforce_frozen = FALSE),
    "missing required columns",
    ignore.case = TRUE
  )

  for (column in c("lfcSE", "mean_count_control", "mean_count_estim")) {
    missing_display_column <- bundle
    missing_display_column$primary_de[[column]] <- NULL
    expect_error(
      validate_public_bundle(missing_display_column, enforce_frozen = FALSE),
      "missing required columns",
      ignore.case = TRUE,
      info = paste("Missing required DE display column:", column)
    )
  }

  invalid_base_mean <- bundle
  invalid_base_mean$primary_de$baseMean[[1L]] <- -1
  expect_error(
    validate_public_bundle(invalid_base_mean, enforce_frozen = FALSE),
    "baseMean",
    fixed = TRUE
  )

  invalid_log2fc <- bundle
  invalid_log2fc$primary_de$log2FoldChange[[1L]] <- Inf
  expect_error(
    validate_public_bundle(invalid_log2fc, enforce_frozen = FALSE),
    "log2FoldChange",
    fixed = TRUE
  )

  invalid_lfcse <- bundle
  invalid_lfcse$primary_de$lfcSE <-
    as.character(invalid_lfcse$primary_de$lfcSE)
  expect_error(
    validate_public_bundle(invalid_lfcse, enforce_frozen = FALSE),
    "lfcSE",
    fixed = TRUE
  )

  invalid_padj <- bundle
  invalid_padj$primary_de$padj[[1L]] <- -1
  expect_error(
    validate_public_bundle(invalid_padj, enforce_frozen = FALSE),
    "padj",
    fixed = TRUE
  )

  invalid_pvalue <- bundle
  invalid_pvalue$primary_de$pvalue[[1L]] <- 2
  expect_error(
    validate_public_bundle(invalid_pvalue, enforce_frozen = FALSE),
    "pvalue",
    fixed = TRUE
  )

  invalid_mean_count <- bundle
  invalid_mean_count$primary_de$mean_count_control[[1L]] <- -1
  expect_error(
    validate_public_bundle(invalid_mean_count, enforce_frozen = FALSE),
    "mean_count_control",
    fixed = TRUE
  )
})

test_that("public-bundle validator requires usable enrichment results", {
  bundle <- synthetic_bundle()

  empty_pathways <- bundle
  empty_pathways$pathways <- empty_pathways$pathways[0, , drop = FALSE]
  empty_pathways$pathway_genes <- empty_pathways$pathway_genes[0, , drop = FALSE]
  expect_error(
    validate_public_bundle(empty_pathways, enforce_frozen = FALSE),
    "unique IDs",
    fixed = TRUE
  )

  duplicate_ids <- bundle
  duplicate_ids$pathways$pathway_id[[2L]] <-
    duplicate_ids$pathways$pathway_id[[1L]]
  expect_error(
    validate_public_bundle(duplicate_ids, enforce_frozen = FALSE),
    "unique IDs",
    fixed = TRUE
  )

  duplicate_labels <- bundle
  duplicate_labels$pathways$label[[2L]] <-
    duplicate_labels$pathways$label[[1L]]
  expect_silent(
    validate_public_bundle(duplicate_labels, enforce_frozen = FALSE)
  )

  invalid_method <- bundle
  invalid_method$pathways$source[[1L]] <- "Other"
  expect_error(
    validate_public_bundle(invalid_method, enforce_frozen = FALSE),
    "pathway statistics",
    fixed = TRUE
  )

  invalid_probability <- bundle
  invalid_probability$pathways$p_adjust[[1L]] <- 2
  expect_error(
    validate_public_bundle(invalid_probability, enforce_frozen = FALSE),
    "pathway statistics",
    fixed = TRUE
  )

  invalid_score <- bundle
  invalid_score$pathways$score[[1L]] <- Inf
  expect_error(
    validate_public_bundle(invalid_score, enforce_frozen = FALSE),
    "pathway statistics",
    fixed = TRUE
  )

  invalid_count <- bundle
  invalid_count$pathways$gene_count[[1L]] <- NA_real_
  expect_error(
    validate_public_bundle(invalid_count, enforce_frozen = FALSE),
    "pathway statistics",
    fixed = TRUE
  )
})

test_that("exporter stores the exact scclrR PFlog result", {
  testthat::skip_if_not_installed("scclrR")
  counts <- methods::as(Matrix::Matrix(
    matrix(c(10, 0, 5, 20, 15, 0), nrow = 2L),
    sparse = TRUE
  ), "dgCMatrix")
  dimnames(counts) <- list(c("Glul", "EGFP"), c("cell_a", "cell_b", "cell_c"))
  reference <- scclrR::normalize_matrix(counts, target = "auto")
  observed <- espiviz_pflog_normalization(counts)

  expect_named(observed, c(
    "method", "target", "log1p", "centered", "sparse", "center", "k",
    "alpha", "package_version", "package_remote_sha"
  ))
  expect_match(observed$method, "scclrR|PFlog", ignore.case = TRUE)
  expect_identical(observed$target, "auto")
  expect_true(observed$log1p)
  expect_true(observed$centered)
  expect_equal(observed$alpha, reference$alpha, tolerance = 1e-12)
  expect_equal(observed$k, reference$k, tolerance = 1e-12)
  expect_equal(observed$sparse, reference$sparse, tolerance = 1e-12)
  expect_equal(observed$center, reference$center, tolerance = 1e-12)
  expect_identical(
    observed$package_version,
    as.character(utils::packageVersion("scclrR"))
  )
  expect_match(observed$package_remote_sha, "^[0-9a-f]{40}$")
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
  expect_equal(
    extracted$normalization$sparse,
    scclrR::normalize_matrix(
      methods::as(Matrix::t(extracted$counts), "dgCMatrix"),
      target = "auto"
    )$sparse,
    tolerance = 1e-12
  )
  expect_identical(
    colnames(extracted$normalization$sparse),
    extracted$cells$cell_id
  )

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

test_that("source extraction does not require a Seurat normalized data layer", {
  object <- make_synthetic_seurat_source()

  expect_false("data" %in% SeuratObject::Layers(object[["RNA"]]))
  expect_silent(
    espiviz_extract_source(object, synthetic_source_contract())
  )
})

test_that("build-data script can be sourced without executing the CLI", {
  build_script <- file.path(repo_root, "scripts", "build-data.R")
  environment <- new.env(parent = globalenv())

  expect_silent(source(build_script, local = environment))
  expect_true(exists("espiviz_run_cli", envir = environment, mode = "function"))
  expect_true(exists(
    "espiviz_format_pflog_status",
    envir = environment,
    mode = "function"
  ))
  status <- environment$espiviz_format_pflog_status(list(
    alpha = 0.25,
    k = 100,
    package_version = "0.1.0",
    package_remote_sha = paste(rep("a", 40L), collapse = "")
  ))
  expect_identical(
    status,
    "PFlog alpha 0.25; K 100; scclrR 0.1.0 (aaaaaaaa)."
  )
})
