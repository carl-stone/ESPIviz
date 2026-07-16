ESPIVIZ_SCHEMA_VERSION <- "1.1.0"

ESPIVIZ_PRODUCTION_EXPECTED <- list(
  schema_version = ESPIVIZ_SCHEMA_VERSION,
  genes = 38394L,
  cells = 3456L,
  clusters = 8L,
  de_rows = 24601L,
  reduction = "umap_pflog_mg_selected_no_filter_cc_dims20",
  cluster_column = "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5"
)

bundle_gene_names <- function(bundle) {
  genes <- bundle$genes
  if (is.data.frame(genes)) {
    if (!"gene" %in% names(genes)) {
      stop("The bundle gene table must contain a gene column.", call. = FALSE)
    }
    genes <- genes$gene
  }
  as.character(genes)
}

read_data_manifest <- function(path) {
  if (!file.exists(path)) {
    stop("The data manifest is missing.", call. = FALSE)
  }
  manifest <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (!is.list(manifest)) {
    stop("The data manifest is invalid.", call. = FALSE)
  }
  manifest
}

public_bundle_url <- function(manifest) {
  url <- manifest$asset_url %||% manifest$url %||% manifest$asset$url
  url <- compact_character(url)
  if (length(url) == 0L) NA_character_ else url[[1L]]
}

manifest_sha256 <- function(manifest) {
  hash <- manifest$sha256 %||% manifest$asset_sha256 %||% manifest$asset$sha256
  hash <- compact_character(hash)
  if (length(hash) == 0L) NA_character_ else tolower(hash[[1L]])
}

manifest_expected <- function(manifest) {
  dimensions <- manifest$dimensions %||% list()
  list(
    schema_version = manifest$schema_version %||% ESPIVIZ_SCHEMA_VERSION,
    genes = manifest$gene_count %||% manifest$genes %||% dimensions$genes,
    cells = manifest$cell_count %||% manifest$cells %||% dimensions$cells,
    clusters = manifest$cluster_count %||%
      manifest$clusters %||%
      dimensions$clusters,
    de_rows = manifest$de_row_count %||%
      manifest$de_rows %||%
      dimensions$primary_de_rows %||%
      dimensions$de_rows,
    source_sha256 = manifest$source_sha256,
    input_sha256 = manifest$input_sha256
  )
}

production_manifest_expected <- function(manifest) {
  declared <- manifest_expected(manifest)
  for (field in c("genes", "cells", "clusters", "de_rows")) {
    if (
      is.null(declared[[field]]) ||
        as.integer(declared[[field]]) != ESPIVIZ_PRODUCTION_EXPECTED[[field]]
    ) {
      stop(
        "The data manifest does not match the fixed production contract.",
        call. = FALSE
      )
    }
  }
  expected <- ESPIVIZ_PRODUCTION_EXPECTED
  expected$source_sha256 <- declared$source_sha256
  expected$input_sha256 <- declared$input_sha256
  expected
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    stop("The data bundle is missing.", call. = FALSE)
  }
  tolower(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

validate_bundle <- function(bundle, expected = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "The Matrix package is required to validate the data bundle.",
      call. = FALSE
    )
  }
  expected <- expected %||% ESPIVIZ_PRODUCTION_EXPECTED
  required <- c(
    "schema_version",
    "data_version",
    "provenance",
    "palette",
    "cells",
    "genes",
    "counts",
    "normalization",
    "primary_de",
    "markers",
    "pathways",
    "pathway_genes",
    "featured_gene_sets"
  )
  missing <- setdiff(required, names(bundle))
  if (length(missing) > 0L) {
    stop(
      "The bundle is missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  unexpected <- setdiff(names(bundle), required)
  if (length(unexpected) > 0L) {
    stop("The bundle contains fields outside the public schema.", call. = FALSE)
  }

  forbidden <- grep(
    "paired|sensitivity|(^|_)qc($|_)|barcode|seurat|pflog|alternative|abundance",
    names(bundle),
    ignore.case = TRUE,
    value = TRUE
  )
  if (length(forbidden) > 0L) {
    stop("The bundle contains excluded fields.", call. = FALSE)
  }

  schema <- as.character(bundle$schema_version)
  expected_schema <- as.character(
    expected$schema_version %||% ESPIVIZ_SCHEMA_VERSION
  )
  if (!identical(schema, expected_schema)) {
    stop(
      "The bundle schema version does not match the application.",
      call. = FALSE
    )
  }

  cell_columns <- c(
    "cell_id",
    "umap_1",
    "umap_2",
    "cluster",
    "condition",
    "Mouse",
    "sample",
    "library_size"
  )
  if (
    !is.data.frame(bundle$cells) ||
      !identical(names(bundle$cells), cell_columns)
  ) {
    stop(
      "The bundle cell table does not match the required schema.",
      call. = FALSE
    )
  }

  if (
    any(grepl("barcode|(^|_)qc($|_)", names(bundle$cells), ignore.case = TRUE))
  ) {
    stop("The bundle cell table contains excluded columns.", call. = FALSE)
  }
  if (
    anyDuplicated(bundle$cells$cell_id) ||
      !all(grepl("^cell_[0-9]+$", bundle$cells$cell_id))
  ) {
    stop(
      "The bundle must use unique sequential internal cell IDs.",
      call. = FALSE
    )
  }
  cell_sequence <- suppressWarnings(as.integer(sub(
    "^cell_",
    "",
    bundle$cells$cell_id
  )))
  if (
    anyNA(cell_sequence) || !identical(cell_sequence, seq_along(cell_sequence))
  ) {
    stop("The internal cell IDs are out of order.", call. = FALSE)
  }
  if (
    any(!is.finite(bundle$cells$library_size)) ||
      any(bundle$cells$library_size <= 0)
  ) {
    stop("Cell library sizes must be positive finite values.", call. = FALSE)
  }

  if (
    !is.data.frame(bundle$genes) ||
      !identical(names(bundle$genes), c("gene", "gene_index")) ||
      !identical(
        as.integer(bundle$genes$gene_index),
        seq_len(nrow(bundle$genes))
      )
  ) {
    stop(
      "The bundle gene table does not match the required schema.",
      call. = FALSE
    )
  }
  genes <- bundle_gene_names(bundle)
  if (
    length(genes) == 0L ||
      anyNA(genes) ||
      any(!nzchar(genes)) ||
      anyDuplicated(casefold_key(genes))
  ) {
    stop("The bundle gene universe is invalid or ambiguous.", call. = FALSE)
  }
  if (
    !inherits(bundle$counts, "dgCMatrix") ||
      !identical(dim(bundle$counts), c(nrow(bundle$cells), length(genes)))
  ) {
    stop(
      "The sparse count matrix dimensions do not match the bundle.",
      call. = FALSE
    )
  }
  if (
    !identical(rownames(bundle$counts), as.character(bundle$cells$cell_id)) ||
      !identical(colnames(bundle$counts), genes)
  ) {
    stop(
      "The sparse count matrix order does not match the bundle.",
      call. = FALSE
    )
  }
  if (
    any(!is.finite(bundle$counts@x)) ||
      any(bundle$counts@x < 0) ||
      any(abs(bundle$counts@x - round(bundle$counts@x)) > 1e-8)
  ) {
    stop(
      "The sparse count matrix must contain finite, non-negative raw counts.",
      call. = FALSE
    )
  }
  if (
    !isTRUE(all.equal(
      as.numeric(Matrix::rowSums(bundle$counts)),
      bundle$cells$library_size,
      tolerance = 0
    ))
  ) {
    stop(
      "Cell library sizes must equal sparse raw-count row sums.",
      call. = FALSE
    )
  }

  de_columns <- c(
    "gene",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj",
    "mean_count_control",
    "mean_count_estim",
    "design"
  )
  allowed_de_columns <- c(
    "gene",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "unshrunkenLog2FoldChange",
    "mean_count_control",
    "mean_count_estim",
    "contrast",
    "design",
    "lfc_shrink_type"
  )
  if (
    !is.data.frame(bundle$primary_de) ||
      !all(de_columns %in% names(bundle$primary_de)) ||
      any(!names(bundle$primary_de) %in% allowed_de_columns) ||
      anyDuplicated(casefold_key(bundle$primary_de$gene)) ||
      any(!bundle$primary_de$gene %in% genes) ||
      !is.numeric(bundle$primary_de$baseMean) ||
      anyNA(bundle$primary_de$baseMean) ||
      any(!is.finite(bundle$primary_de$baseMean)) ||
      any(bundle$primary_de$baseMean < 0) ||
      !is.numeric(bundle$primary_de$log2FoldChange) ||
      anyNA(bundle$primary_de$log2FoldChange) ||
      any(!is.finite(bundle$primary_de$log2FoldChange)) ||
      !is.numeric(bundle$primary_de$lfcSE) ||
      anyNA(bundle$primary_de$lfcSE) ||
      any(!is.finite(bundle$primary_de$lfcSE)) ||
      any(bundle$primary_de$lfcSE < 0) ||
      !is.numeric(bundle$primary_de$pvalue) ||
      any(
        !is.na(bundle$primary_de$pvalue) &
          (
            !is.finite(bundle$primary_de$pvalue) |
              bundle$primary_de$pvalue < 0 |
              bundle$primary_de$pvalue > 1
          )
      ) ||
      !is.numeric(bundle$primary_de$padj) ||
      any(
        !is.na(bundle$primary_de$padj) &
          (
            !is.finite(bundle$primary_de$padj) |
              bundle$primary_de$padj < 0 |
              bundle$primary_de$padj > 1
          )
      ) ||
      !is.numeric(bundle$primary_de$mean_count_control) ||
      anyNA(bundle$primary_de$mean_count_control) ||
      any(!is.finite(bundle$primary_de$mean_count_control)) ||
      any(bundle$primary_de$mean_count_control < 0) ||
      !is.numeric(bundle$primary_de$mean_count_estim) ||
      anyNA(bundle$primary_de$mean_count_estim) ||
      any(!is.finite(bundle$primary_de$mean_count_estim)) ||
      any(bundle$primary_de$mean_count_estim < 0) ||
      !identical(
        unique(as.character(bundle$primary_de$design)),
        "primary_unpaired_condition"
      )
  ) {
    stop("The primary differential-expression table is invalid.", call. = FALSE)
  }
  marker_columns <- c(
    "cluster",
    "gene",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "rank"
  )
  if (
    !is.data.frame(bundle$markers) ||
      !identical(names(bundle$markers), marker_columns) ||
      any(!bundle$markers$gene %in% genes)
  ) {
    stop("The marker table is invalid.", call. = FALSE)
  }
  pathway_columns <- c(
    "pathway_id",
    "label",
    "source",
    "direction",
    "description",
    "p_value",
    "p_adjust",
    "score",
    "gene_count"
  )
  if (
    !is.data.frame(bundle$pathways) ||
      !identical(names(bundle$pathways), pathway_columns) ||
      nrow(bundle$pathways) == 0L ||
      anyNA(bundle$pathways$pathway_id) ||
      any(!nzchar(as.character(bundle$pathways$pathway_id))) ||
      anyDuplicated(as.character(bundle$pathways$pathway_id)) ||
      anyNA(bundle$pathways$label) ||
      any(!nzchar(as.character(bundle$pathways$label))) ||
      anyDuplicated(as.character(bundle$pathways$label)) ||
      any(!as.character(bundle$pathways$source) %in% c("GSEA", "ORA")) ||
      any(
        !as.character(bundle$pathways$direction) %in% c("Control", "E-Stim")
      ) ||
      !is.numeric(bundle$pathways$p_adjust) ||
      anyNA(bundle$pathways$p_adjust) ||
      any(!is.finite(bundle$pathways$p_adjust)) ||
      any(bundle$pathways$p_adjust < 0 | bundle$pathways$p_adjust > 1) ||
      !is.numeric(bundle$pathways$score) ||
      anyNA(bundle$pathways$score) ||
      any(!is.finite(bundle$pathways$score)) ||
      !is.numeric(bundle$pathways$gene_count) ||
      anyNA(bundle$pathways$gene_count) ||
      any(!is.finite(bundle$pathways$gene_count)) ||
      any(bundle$pathways$gene_count < 0) ||
      !is.data.frame(bundle$pathway_genes) ||
      !identical(names(bundle$pathway_genes), c("pathway_id", "gene"))
  ) {
    stop("The featured pathway tables are invalid.", call. = FALSE)
  }
  if (
    any(!bundle$pathway_genes$pathway_id %in% bundle$pathways$pathway_id) ||
      any(!bundle$pathway_genes$gene %in% genes)
  ) {
    stop(
      "The featured pathway genes are outside the public schema.",
      call. = FALSE
    )
  }

  normalization <- bundle$normalization
  normalization_names <- c(
    "method", "target", "log1p", "centered", "sparse", "center", "k",
    "alpha", "package_version", "package_remote_sha"
  )
  if (
    !is.list(normalization) ||
      !identical(names(normalization), normalization_names) ||
      !identical(normalization$method, "scclrR::normalize_matrix") ||
      !identical(normalization$target, "auto") ||
      !isTRUE(normalization$log1p) ||
      !isTRUE(normalization$centered) ||
      !inherits(normalization$sparse, "dgCMatrix") ||
      !identical(dim(normalization$sparse), rev(dim(bundle$counts))) ||
      !identical(rownames(normalization$sparse), genes) ||
      !identical(colnames(normalization$sparse), bundle$cells$cell_id) ||
      any(!is.finite(normalization$sparse@x)) ||
      !is.numeric(normalization$center) ||
      length(normalization$center) != nrow(bundle$cells) ||
      !identical(names(normalization$center), bundle$cells$cell_id) ||
      any(!is.finite(normalization$center)) ||
      !is.numeric(normalization$k) ||
      length(normalization$k) != 1L ||
      !is.finite(normalization$k) ||
      normalization$k <= 0 ||
      !is.numeric(normalization$alpha) ||
      length(normalization$alpha) != 1L ||
      !is.finite(normalization$alpha) ||
      normalization$alpha <= 0 ||
      !is.character(normalization$package_version) ||
      length(normalization$package_version) != 1L ||
      !nzchar(normalization$package_version) ||
      !is.character(normalization$package_remote_sha) ||
      length(normalization$package_remote_sha) != 1L ||
      !grepl("^[0-9a-f]{40}$", normalization$package_remote_sha)
  ) {
    stop("The exported scclrR normalization state is invalid.", call. = FALSE)
  }
  if (!is.list(bundle$featured_gene_sets)) {
    stop("Featured gene sets must be stored as a named list.", call. = FALSE)
  }
  featured_genes <- unlist(bundle$featured_gene_sets, use.names = FALSE)
  if (any(!featured_genes %in% genes)) {
    stop("A featured gene set contains an unknown gene.", call. = FALSE)
  }

  actual <- list(
    genes = length(genes),
    cells = nrow(bundle$cells),
    clusters = length(unique(as.character(bundle$cells$cluster))),
    de_rows = nrow(bundle$primary_de)
  )
  for (field in names(actual)) {
    target <- expected[[field]]
    if (
      !is.null(target) &&
        length(target) > 0L &&
        !is.na(target) &&
        actual[[field]] != as.integer(target)
    ) {
      stop(
        sprintf(
          "The bundle %s count is %s; expected %s.",
          field,
          actual[[field]],
          target
        ),
        call. = FALSE
      )
    }
  }

  expected_source_hash <- compact_character(expected$source_sha256)
  if (length(expected_source_hash) > 0L) {
    provenance_hash <- compact_character(bundle$provenance$source_sha256)
    nested_hash <- compact_character(
      bundle$provenance$inputs$source_object$sha256
    )
    if (
      length(provenance_hash) == 0L ||
        length(nested_hash) == 0L ||
        !identical(
          tolower(provenance_hash[[1L]]),
          tolower(expected_source_hash[[1L]])
        ) ||
        !identical(
          tolower(nested_hash[[1L]]),
          tolower(expected_source_hash[[1L]])
        )
    ) {
      stop(
        "The frozen source checksum does not match the data manifest.",
        call. = FALSE
      )
    }
  }
  if (
    !is.null(expected$reduction) &&
      !identical(as.character(bundle$provenance$reduction), expected$reduction)
  ) {
    stop(
      "The bundle reduction does not match the fixed production contract.",
      call. = FALSE
    )
  }
  if (
    !is.null(expected$cluster_column) &&
      !identical(
        as.character(bundle$provenance$cluster_column),
        expected$cluster_column
      )
  ) {
    stop(
      "The bundle cluster column does not match the fixed production contract.",
      call. = FALSE
    )
  }
  expected_inputs <- expected$input_sha256 %||% list()
  if (length(expected_inputs) > 0L) {
    for (input_name in names(expected_inputs)) {
      expected_input_hash <- compact_character(expected_inputs[[input_name]])
      provenance_input_hash <- compact_character(
        bundle$provenance$inputs[[input_name]]$sha256
      )
      if (
        length(expected_input_hash) > 0L &&
          (length(provenance_input_hash) == 0L ||
            !identical(
              tolower(provenance_input_hash[[1L]]),
              tolower(expected_input_hash[[1L]])
            ))
      ) {
        stop(
          "An imported result checksum does not match the data manifest.",
          call. = FALSE
        )
      }
    }
  }

  serialized_strings <- unlist(bundle$provenance, use.names = FALSE)
  serialized_strings <- serialized_strings[is.character(serialized_strings)]
  if (
    any(grepl(
      "/Users/|Box-Box|/Box/|Library/CloudStorage|file://|\\\\Box\\\\",
      serialized_strings,
      fixed = FALSE
    ))
  ) {
    stop(
      "The bundle provenance contains a private storage path.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

download_bundle_asset <- function(url, data_version) {
  if (!grepl("^https://", url)) {
    stop("The bundle URL must use HTTPS.", call. = FALSE)
  }
  cache_dir <- file.path(tempdir(), "espiviz", safe_filename(data_version))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(cache_dir, basename(sub("[?].*$", "", url)))
  if (!file.exists(destination)) {
    temporary <- paste0(destination, ".partial")
    on.exit(unlink(temporary), add = TRUE)
    status <- utils::download.file(url, temporary, mode = "wb", quiet = TRUE)
    if (!identical(status, 0L)) {
      stop("The data bundle could not be downloaded.", call. = FALSE)
    }
    if (!file.rename(temporary, destination)) {
      stop("The downloaded data bundle could not be cached.", call. = FALSE)
    }
  }
  destination
}

load_bundle <- function(
  path,
  manifest,
  expected = NULL,
  expected_sha256 = NULL
) {
  if (is.character(manifest) && length(manifest) == 1L) {
    manifest <- read_data_manifest(manifest)
  }
  if (!is.list(manifest)) {
    stop("The data manifest is invalid.", call. = FALSE)
  }
  path <- path.expand(path)
  hash_override <- compact_character(expected_sha256)
  hash_environment <- compact_character(Sys.getenv("ESPIVIZ_DATA_SHA256", ""))
  expected_hash <- if (length(hash_override) > 0L) {
    hash_override
  } else if (length(hash_environment) > 0L) {
    hash_environment
  } else {
    compact_character(manifest_sha256(manifest))
  }
  if (
    length(expected_hash) == 0L ||
      !grepl("^[0-9a-fA-F]{64}$", expected_hash[[1L]])
  ) {
    stop(
      "The data manifest does not contain a valid SHA-256 checksum.",
      call. = FALSE
    )
  }
  if (!identical(sha256_file(path), tolower(expected_hash[[1L]]))) {
    stop("The data bundle checksum does not match the manifest.", call. = FALSE)
  }

  bundle <- readRDS(path)
  validate_bundle(bundle, expected %||% manifest_expected(manifest))
  attr(bundle, "bundle_path") <- normalizePath(path, mustWork = TRUE)
  attr(bundle, "data_manifest") <- manifest
  bundle
}

load_public_bundle <- function(
  manifest_path = file.path("app", "data-manifest.json"),
  path = NULL,
  expected_sha256 = NULL
) {
  manifest <- read_data_manifest(manifest_path)
  path <- compact_character(path %||% Sys.getenv("ESPIVIZ_DATA_PATH", ""))
  if (length(path) > 0L) {
    path <- path.expand(path[[1L]])
  } else {
    url <- public_bundle_url(manifest)
    if (is.na(url)) {
      stop("The data manifest does not provide a bundle URL.", call. = FALSE)
    }
    path <- download_bundle_asset(url, manifest$data_version %||% "current")
  }
  load_bundle(
    path = path,
    manifest = manifest,
    expected = production_manifest_expected(manifest),
    expected_sha256 = expected_sha256
  )
}

load_bundle_once <- local({
  cache <- NULL
  cache_key <- NULL
  function(manifest_path, path = NULL, expected_sha256 = NULL) {
    key <- paste(
      normalizePath(manifest_path, mustWork = FALSE),
      path %||% Sys.getenv("ESPIVIZ_DATA_PATH", ""),
      expected_sha256 %||% Sys.getenv("ESPIVIZ_DATA_SHA256", ""),
      sep = "|"
    )
    if (is.null(cache) || !identical(cache_key, key)) {
      cache <<- load_public_bundle(
        path = path,
        manifest_path = manifest_path,
        expected_sha256 = expected_sha256
      )
      cache_key <<- key
    }
    cache
  }
})
