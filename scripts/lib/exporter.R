`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

espiviz_abort <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

espiviz_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    espiviz_abort(
      "Package '%s' is required to build the ESPIviz data bundle.",
      package
    )
  }
}

espiviz_contract <- function() {
  list(
    schema_version = "1.1.0",
    data_version = "1.1.0",
    reduction = "umap_pflog_mg_selected_no_filter_cc_dims20",
    cluster_column = "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5",
    genes = 38394L,
    cells = 3456L,
    clusters = as.character(seq_len(8L)),
    conditions = c("p27CKO", "p27CKO +EStim"),
    primary_de_rows = 24601L,
    primary_de_design = "primary_unpaired_condition",
    primary_de_contrast = "estim_vs_control",
    top_markers_per_cluster = 25L
  )
}

espiviz_sha256 <- function(path) {
  espiviz_require("digest")
  if (!file.exists(path) || dir.exists(path)) {
    espiviz_abort("Input file does not exist: %s", path)
  }
  digest::digest(file = path, algo = "sha256")
}

espiviz_is_absolute_path <- function(path) {
  grepl("^/", path) || grepl("^[A-Za-z]:[/\\\\]", path)
}

espiviz_read_source_manifest <- function(path) {
  espiviz_require("jsonlite")
  if (!file.exists(path)) {
    espiviz_abort("Source manifest does not exist: %s", path)
  }

  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  required <- c(
    "source_object",
    "primary_de",
    "markers",
    "ora_up",
    "ora_down",
    "gsea",
    "featured_genes_config",
    "featured_pathways_config"
  )
  allowed <- c("manifest_version", required)

  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0L) {
    espiviz_abort(
      "Source manifest is missing required entries: %s",
      paste(missing, collapse = ", ")
    )
  }

  extra <- setdiff(names(manifest), allowed)
  if (length(extra) > 0L) {
    espiviz_abort(
      "Source manifest contains unsupported entries: %s",
      paste(extra, collapse = ", ")
    )
  }

  manifest_version <- manifest$manifest_version %||% ""
  if (!identical(manifest_version, "1.0.0")) {
    espiviz_abort("Source manifest_version must be '1.0.0'.")
  }

  for (name in required) {
    item <- manifest[[name]]
    item_path <- item$path %||% ""
    if (
      !is.character(item_path) || length(item_path) != 1L || !nzchar(item_path)
    ) {
      espiviz_abort("Source manifest entry '%s' must contain one path.", name)
    }
    item_path <- path.expand(item_path)
    if (!espiviz_is_absolute_path(item_path)) {
      espiviz_abort("Source manifest path '%s' must be absolute.", name)
    }
    if (!file.exists(item_path) || dir.exists(item_path)) {
      espiviz_abort(
        "Source manifest path '%s' does not exist: %s",
        name,
        item_path
      )
    }
    manifest[[name]]$path <- normalizePath(item_path, mustWork = TRUE)
  }

  source_hash <- manifest$source_object$sha256 %||% ""
  if (
    !is.character(source_hash) ||
      length(source_hash) != 1L ||
      !grepl("^[0-9a-fA-F]{64}$", source_hash)
  ) {
    espiviz_abort("source_object.sha256 must be one complete SHA-256 hash.")
  }
  manifest$source_object$sha256 <- tolower(source_hash)
  manifest
}

espiviz_input_provenance <- function(manifest) {
  input_names <- c(
    "source_object",
    "primary_de",
    "markers",
    "ora_up",
    "ora_down",
    "gsea",
    "featured_genes_config",
    "featured_pathways_config"
  )
  stats::setNames(
    lapply(input_names, function(name) {
      list(
        sha256 = espiviz_sha256(manifest[[name]]$path)
      )
    }),
    input_names
  )
}

espiviz_read_table <- function(path, format = c("tsv", "csv")) {
  format <- match.arg(format)
  reader <- if (identical(format, "csv")) utils::read.csv else utils::read.delim
  data <- reader(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", ""),
    quote = if (identical(format, "tsv")) "" else "\"",
    comment.char = ""
  )
  if (nrow(data) > 0L) {
    data <- data[rowSums(!is.na(data)) > 0L, , drop = FALSE]
    rownames(data) <- NULL
  }
  data
}

espiviz_require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    espiviz_abort(
      "%s is missing required columns: %s",
      label,
      paste(missing, collapse = ", ")
    )
  }
  invisible(data)
}

espiviz_pflog_normalization <- function(counts) {
  espiviz_require("Matrix")
  espiviz_require("scclrR")
  if (!inherits(counts, "sparseMatrix")) {
    espiviz_abort("PFlog normalization requires a sparse count matrix.")
  }

  normalized <- scclrR::normalize_matrix(
    counts,
    target = "auto",
    log1p = TRUE,
    center = TRUE
  )
  if (
    !is.list(normalized) ||
      !inherits(normalized$sparse, "sparseMatrix") ||
      !identical(dim(normalized$sparse), dim(counts)) ||
      !identical(dimnames(normalized$sparse), dimnames(counts)) ||
      length(normalized$center) != ncol(counts) ||
      any(!is.finite(normalized$center)) ||
      length(normalized$alpha) != 1L ||
      !is.finite(normalized$alpha) ||
      normalized$alpha <= 0 ||
      length(normalized$k) != 1L ||
      !is.finite(normalized$k) ||
      normalized$k <= 0
  ) {
    espiviz_abort("scclrR returned an invalid PFlog normalization result.")
  }

  description <- utils::packageDescription("scclrR")
  remote_sha <- description$RemoteSha %||%
    description$GithubSHA1 %||%
    NA_character_
  list(
    method = "scclrR::normalize_matrix",
    target = "auto",
    log1p = TRUE,
    centered = TRUE,
    sparse = normalized$sparse,
    center = unname(as.numeric(normalized$center)),
    k = unname(as.numeric(normalized$k)),
    alpha = unname(as.numeric(normalized$alpha)),
    package_version = as.character(utils::packageVersion("scclrR")),
    package_remote_sha = as.character(remote_sha)
  )
}

espiviz_extract_source <- function(object, contract = espiviz_contract()) {
  espiviz_require("Matrix")
  espiviz_require("SeuratObject")

  if (!inherits(object, "Seurat")) {
    espiviz_abort("source_object must contain a Seurat object.")
  }
  if (!"RNA" %in% SeuratObject::Assays(object)) {
    espiviz_abort("source_object is missing the RNA assay.")
  }
  if (!contract$reduction %in% SeuratObject::Reductions(object)) {
    espiviz_abort(
      "source_object is missing reduction '%s'.",
      contract$reduction
    )
  }

  object_dimensions <- dim(object)
  if (
    !identical(as.integer(object_dimensions), c(contract$genes, contract$cells))
  ) {
    espiviz_abort(
      "source_object dimensions are %d genes x %d cells; expected %d x %d.",
      object_dimensions[[1L]],
      object_dimensions[[2L]],
      contract$genes,
      contract$cells
    )
  }

  assay <- object[["RNA"]]
  layers <- SeuratObject::Layers(assay)
  if (!"counts" %in% layers) {
    espiviz_abort("RNA assay must contain a counts layer.")
  }
  counts <- SeuratObject::LayerData(assay, layer = "counts")
  if (!inherits(counts, "sparseMatrix")) {
    espiviz_abort("RNA counts must be a sparse matrix.")
  }
  if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
    espiviz_abort("RNA counts contain non-finite or negative values.")
  }
  if (any(abs(counts@x - round(counts@x)) > 1e-8)) {
    espiviz_abort("RNA counts contain non-integer values.")
  }

  genes <- rownames(counts)
  source_cells <- colnames(counts)
  if (
    is.null(genes) ||
      is.null(source_cells) ||
      anyNA(genes) ||
      anyNA(source_cells) ||
      any(!nzchar(genes)) ||
      any(!nzchar(source_cells)) ||
      anyDuplicated(genes) ||
      anyDuplicated(source_cells)
  ) {
    espiviz_abort(
      "RNA counts must have complete, unique feature and cell names."
    )
  }
  if (anyDuplicated(toupper(genes))) {
    espiviz_abort("Gene names are not unique under case-insensitive lookup.")
  }
  if (
    !identical(genes, rownames(object)) ||
      !identical(source_cells, colnames(object))
  ) {
    espiviz_abort("RNA counts do not match the Seurat feature and cell order.")
  }

  umap <- SeuratObject::Embeddings(object, reduction = contract$reduction)
  if (!identical(rownames(umap), source_cells)) {
    espiviz_abort("UMAP cell order does not match the RNA counts cell order.")
  }
  if (!identical(dim(umap), c(contract$cells, 2L))) {
    espiviz_abort(
      "Reduction '%s' must contain exactly %d cells and two coordinates.",
      contract$reduction,
      contract$cells
    )
  }
  if (any(!is.finite(umap))) {
    espiviz_abort("UMAP coordinates contain non-finite values.")
  }

  metadata <- object[[]]
  if (!identical(rownames(metadata), source_cells)) {
    espiviz_abort(
      "Metadata cell order does not match the RNA counts cell order."
    )
  }
  espiviz_require_columns(
    metadata,
    c(contract$cluster_column, "Condition", "Mouse", "sample_id"),
    "source_object metadata"
  )
  cluster <- as.character(metadata[[contract$cluster_column]])
  observed_clusters <- sort(unique(cluster))
  if (!identical(observed_clusters, contract$clusters)) {
    espiviz_abort(
      "Cluster column '%s' contains [%s]; expected [%s].",
      contract$cluster_column,
      paste(observed_clusters, collapse = ", "),
      paste(contract$clusters, collapse = ", ")
    )
  }
  condition <- as.character(metadata[["Condition"]])
  if (!identical(sort(unique(condition)), sort(contract$conditions))) {
    espiviz_abort(
      "Condition metadata contains [%s]; expected [%s].",
      paste(sort(unique(condition)), collapse = ", "),
      paste(sort(contract$conditions), collapse = ", ")
    )
  }

  library_size <- as.numeric(Matrix::colSums(counts))
  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    espiviz_abort("Cell library sizes must be finite and greater than zero.")
  }
  normalization <- espiviz_pflog_normalization(counts)

  internal_ids <- sprintf(
    paste0("cell_%0", nchar(as.character(contract$cells)), "d"),
    seq_len(contract$cells)
  )
  colnames(normalization$sparse) <- internal_ids
  names(normalization$center) <- internal_ids
  cells <- data.frame(
    cell_id = internal_ids,
    umap_1 = unname(umap[, 1L]),
    umap_2 = unname(umap[, 2L]),
    cluster = as.integer(cluster),
    condition = condition,
    Mouse = as.character(metadata[["Mouse"]]),
    sample = as.character(metadata[["sample_id"]]),
    library_size = library_size,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(cells) <- NULL

  cell_gene_counts <- Matrix::t(counts)
  if (!inherits(cell_gene_counts, "dgCMatrix")) {
    cell_gene_counts <- methods::as(cell_gene_counts, "dgCMatrix")
  }
  dimnames(cell_gene_counts) <- list(internal_ids, genes)

  list(
    cells = cells,
    genes = data.frame(
      gene = genes,
      gene_index = seq_along(genes),
      stringsAsFactors = FALSE
    ),
    counts = cell_gene_counts,
    normalization = normalization
  )
}

espiviz_read_primary_de <- function(
  path,
  genes,
  contract = espiviz_contract()
) {
  data <- espiviz_read_table(path, "tsv")
  required <- c(
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
  espiviz_require_columns(data, required, "Primary DE table")
  if (nrow(data) != contract$primary_de_rows) {
    espiviz_abort(
      "Primary DE table has %d rows; expected %d.",
      nrow(data),
      contract$primary_de_rows
    )
  }
  if (anyNA(data$gene) || any(!nzchar(data$gene)) || anyDuplicated(data$gene)) {
    espiviz_abort("Primary DE genes must be complete and unique.")
  }
  missing_genes <- setdiff(data$gene, genes)
  if (length(missing_genes) > 0L) {
    espiviz_abort(
      "Primary DE table contains genes absent from the source object: %s",
      paste(utils::head(missing_genes, 10L), collapse = ", ")
    )
  }
  if (!identical(unique(data$design), contract$primary_de_design)) {
    espiviz_abort(
      "Primary DE design must be '%s' only.",
      contract$primary_de_design
    )
  }
  if (!identical(unique(data$contrast), contract$primary_de_contrast)) {
    espiviz_abort(
      "Primary DE contrast must be '%s' only.",
      contract$primary_de_contrast
    )
  }
  if (any(grepl("(^|_)paired($|_)", data$design, ignore.case = TRUE))) {
    espiviz_abort("Paired-sensitivity DE is not allowed in the public bundle.")
  }
  rownames(data) <- NULL
  data
}

espiviz_read_markers <- function(path, genes, contract = espiviz_contract()) {
  data <- espiviz_read_table(path, "csv")
  required <- c(
    "gene",
    "cluster",
    "rank_within_cluster",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj"
  )
  espiviz_require_columns(data, required, "Marker table")
  cluster <- as.character(data$cluster)
  if (!identical(sort(unique(cluster)), contract$clusters)) {
    espiviz_abort(
      "Marker table does not contain exactly the eight final clusters."
    )
  }
  if (anyNA(data$gene) || any(!data$gene %in% genes)) {
    espiviz_abort("Marker table contains missing or unknown genes.")
  }
  rank <- as.integer(data$rank_within_cluster)
  if (anyNA(rank) || any(rank < 1L)) {
    espiviz_abort("Marker ranks must be positive integers.")
  }
  if (anyDuplicated(paste(cluster, rank, sep = "::"))) {
    espiviz_abort("Marker ranks must be unique within each cluster.")
  }

  keep <- unlist(
    lapply(contract$clusters, function(current_cluster) {
      index <- which(cluster == current_cluster)
      index[order(rank[index])][seq_len(min(
        length(index),
        contract$top_markers_per_cluster
      ))]
    }),
    use.names = FALSE
  )
  markers <- data.frame(
    cluster = as.integer(cluster[keep]),
    gene = data$gene[keep],
    avg_log2FC = as.numeric(data$avg_log2FC[keep]),
    pct.1 = as.numeric(data$pct.1[keep]),
    pct.2 = as.numeric(data$pct.2[keep]),
    p_val_adj = as.numeric(data$p_val_adj[keep]),
    rank = rank[keep],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  markers <- markers[order(markers$cluster, markers$rank), , drop = FALSE]
  rownames(markers) <- NULL
  markers
}

espiviz_read_featured_gene_sets <- function(path, genes) {
  config <- espiviz_read_table(path, "csv")
  required <- c("set_id", "set_label", "gene", "display_order")
  espiviz_require_columns(config, required, "Featured gene config")
  if (nrow(config) == 0L) {
    espiviz_abort("Featured gene config must contain at least one gene.")
  }
  if (
    anyNA(config$set_id) ||
      anyNA(config$set_label) ||
      anyNA(config$gene) ||
      any(!nzchar(config$set_id)) ||
      any(!nzchar(config$set_label)) ||
      any(!nzchar(config$gene))
  ) {
    espiviz_abort("Featured gene config contains blank values.")
  }
  if (any(!grepl("^[a-z0-9_]+$", config$set_id))) {
    espiviz_abort("Featured gene set_id values must use lower snake case.")
  }
  if (anyDuplicated(paste(config$set_id, config$gene, sep = "::"))) {
    espiviz_abort("Featured gene config contains duplicate genes within a set.")
  }
  unknown <- setdiff(config$gene, genes)
  if (length(unknown) > 0L) {
    espiviz_abort(
      "Featured gene config contains genes absent from the source object: %s",
      paste(unknown, collapse = ", ")
    )
  }

  set_metadata <- unique(config[c("set_id", "set_label", "display_order")])
  if (anyDuplicated(set_metadata$set_id)) {
    espiviz_abort(
      "Each featured gene set_id must have one label and display order."
    )
  }
  if (anyDuplicated(set_metadata$set_label)) {
    espiviz_abort("Featured gene set labels must be unique.")
  }
  set_metadata <- set_metadata[order(as.numeric(set_metadata$display_order)), ]
  sets <- lapply(set_metadata$set_id, function(id) {
    config$gene[config$set_id == id]
  })
  stats::setNames(sets, set_metadata$set_label)
}

espiviz_map_entrez_symbols <- function(entrez_ids) {
  espiviz_require("AnnotationDbi")
  espiviz_require("org.Mm.eg.db")
  database <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
  mapped <- suppressMessages(
    AnnotationDbi::mapIds(
      database,
      keys = unique(entrez_ids),
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
  )
  unname(mapped[match(entrez_ids, names(mapped))])
}

espiviz_read_featured_pathways <- function(
  config_path,
  result_paths,
  genes
) {
  config <- espiviz_read_table(config_path, "csv")
  required <- c(
    "pathway_id",
    "label",
    "source",
    "result_id",
    "direction",
    "display_order"
  )
  espiviz_require_columns(config, required, "Featured pathway config")
  if (nrow(config) == 0L) {
    espiviz_abort("Featured pathway config must contain at least one term.")
  }
  if (anyDuplicated(config$pathway_id) || anyDuplicated(config$label)) {
    espiviz_abort("Featured pathway IDs and labels must be unique.")
  }
  if (any(!grepl("^[a-z0-9_]+$", config$pathway_id))) {
    espiviz_abort("Featured pathway_id values must use lower snake case.")
  }
  allowed_sources <- c("ora_up", "ora_down", "gsea")
  if (any(!config$source %in% allowed_sources)) {
    espiviz_abort("Featured pathway source must be ora_up, ora_down, or gsea.")
  }
  if (any(!config$direction %in% c("E-Stim", "Control"))) {
    espiviz_abort("Featured pathway direction must be E-Stim or Control.")
  }

  result_tables <- lapply(allowed_sources, function(source) {
    data <- espiviz_read_table(result_paths[[source]], "tsv")
    common <- c("ID", "Description", "pvalue", "p.adjust")
    extra <- if (identical(source, "gsea")) {
      c("NES", "setSize", "core_enrichment")
    } else {
      c("FoldEnrichment", "Count", "geneID", "direction")
    }
    espiviz_require_columns(
      data,
      c(common, extra),
      sprintf("%s pathway table", source)
    )
    if (anyDuplicated(data$ID)) {
      espiviz_abort("%s pathway table contains duplicate term IDs.", source)
    }
    data
  })
  names(result_tables) <- allowed_sources

  config <- config[order(as.numeric(config$display_order)), , drop = FALSE]
  pathway_rows <- vector("list", nrow(config))
  pathway_gene_rows <- vector("list", nrow(config))

  for (index in seq_len(nrow(config))) {
    item <- config[index, , drop = FALSE]
    source <- item$source[[1L]]
    source_table <- result_tables[[source]]
    source_row <- source_table[
      source_table$ID == item$result_id[[1L]],
      ,
      drop = FALSE
    ]
    if (nrow(source_row) != 1L) {
      espiviz_abort(
        "Configured pathway '%s' was not found exactly once in %s.",
        item$result_id[[1L]],
        source
      )
    }
    if (is.na(source_row$p.adjust[[1L]]) || source_row$p.adjust[[1L]] >= 0.05) {
      espiviz_abort(
        "Configured pathway '%s' is not significant at adjusted P < 0.05.",
        item$result_id[[1L]]
      )
    }

    if (identical(source, "gsea")) {
      score <- as.numeric(source_row$NES[[1L]])
      expected_direction <- if (score > 0) "E-Stim" else "Control"
      entrez <- strsplit(source_row$core_enrichment[[1L]], "/", fixed = TRUE)[[
        1L
      ]]
      pathway_genes <- espiviz_map_entrez_symbols(entrez)
      gene_count <- as.integer(source_row$setSize[[1L]])
      analysis <- "GSEA"
    } else {
      score <- as.numeric(source_row$FoldEnrichment[[1L]])
      expected_direction <- if (identical(source, "ora_up")) {
        "E-Stim"
      } else {
        "Control"
      }
      pathway_genes <- strsplit(source_row$geneID[[1L]], "/", fixed = TRUE)[[
        1L
      ]]
      gene_count <- as.integer(source_row$Count[[1L]])
      analysis <- "ORA"
    }
    if (!identical(item$direction[[1L]], expected_direction)) {
      espiviz_abort(
        "Configured direction for pathway '%s' does not match its result.",
        item$result_id[[1L]]
      )
    }
    pathway_genes <- unique(pathway_genes[!is.na(pathway_genes)])
    pathway_genes <- pathway_genes[pathway_genes %in% genes]
    if (length(pathway_genes) == 0L) {
      espiviz_abort(
        "Configured pathway '%s' has no genes in the source object.",
        item$result_id[[1L]]
      )
    }

    pathway_rows[[index]] <- data.frame(
      pathway_id = item$pathway_id[[1L]],
      label = item$label[[1L]],
      source = analysis,
      direction = item$direction[[1L]],
      description = source_row$Description[[1L]],
      p_value = as.numeric(source_row$pvalue[[1L]]),
      p_adjust = as.numeric(source_row$p.adjust[[1L]]),
      score = score,
      gene_count = gene_count,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    pathway_gene_rows[[index]] <- data.frame(
      pathway_id = item$pathway_id[[1L]],
      gene = pathway_genes,
      stringsAsFactors = FALSE
    )
  }

  list(
    pathways = do.call(rbind, pathway_rows),
    pathway_genes = do.call(rbind, pathway_gene_rows)
  )
}

espiviz_palette <- function() {
  list(
    expression = c(low = "#2166AC", mid = "#BDBDBD", high = "#E31A8C"),
    condition = c("p27CKO" = "#2166AC", "p27CKO +EStim" = "#E31A8C"),
    cluster = c(
      "1" = "#3B5B92",
      "2" = "#5F8F77",
      "3" = "#C07A4A",
      "4" = "#876B9E",
      "5" = "#C0A23E",
      "6" = "#4D8E9B",
      "7" = "#B65C70",
      "8" = "#6F6F6F"
    )
  )
}

espiviz_software_versions <- function() {
  packages <- c(
    "Matrix",
    "SeuratObject",
    "scclrR",
    "AnnotationDbi",
    "org.Mm.eg.db"
  )
  package_versions <- stats::setNames(
    lapply(packages, function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }),
    packages
  )
  c(list(R = as.character(getRversion())), package_versions)
}

espiviz_collect_characters <- function(value) {
  if (is.character(value)) {
    return(value)
  }
  if (is.data.frame(value)) {
    return(unlist(lapply(value, espiviz_collect_characters), use.names = FALSE))
  }
  if (is.list(value)) {
    return(unlist(lapply(value, espiviz_collect_characters), use.names = FALSE))
  }
  character()
}

espiviz_validate_public_bundle <- function(
  bundle,
  contract = espiviz_contract(),
  enforce_frozen = TRUE
) {
  expected_names <- c(
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
  if (!identical(names(bundle), expected_names)) {
    espiviz_abort(
      "Bundle fields must be exactly: %s",
      paste(expected_names, collapse = ", ")
    )
  }
  if (!identical(bundle$schema_version, contract$schema_version)) {
    espiviz_abort("Bundle schema_version does not match the public contract.")
  }
  if (!identical(bundle$data_version, contract$data_version)) {
    espiviz_abort("Bundle data_version does not match the public contract.")
  }

  expected_cell_columns <- c(
    "cell_id",
    "umap_1",
    "umap_2",
    "cluster",
    "condition",
    "Mouse",
    "sample",
    "library_size"
  )
  if (!identical(names(bundle$cells), expected_cell_columns)) {
    espiviz_abort("Bundle cells contain unsupported or missing columns.")
  }
  if (!identical(names(bundle$genes), c("gene", "gene_index"))) {
    espiviz_abort("Bundle genes must contain gene and gene_index only.")
  }
  if (!inherits(bundle$counts, "dgCMatrix")) {
    espiviz_abort("Bundle counts must be a cells x genes dgCMatrix.")
  }
  if (
    !identical(dim(bundle$counts), c(nrow(bundle$cells), nrow(bundle$genes)))
  ) {
    espiviz_abort("Bundle counts dimensions do not match cells and genes.")
  }
  if (!identical(rownames(bundle$counts), bundle$cells$cell_id)) {
    espiviz_abort("Bundle counts row order does not match internal cell IDs.")
  }
  if (!identical(colnames(bundle$counts), bundle$genes$gene)) {
    espiviz_abort("Bundle counts column order does not match genes.")
  }
  expected_ids <- sprintf(
    paste0("cell_%0", nchar(as.character(nrow(bundle$cells))), "d"),
    seq_len(nrow(bundle$cells))
  )
  if (!identical(bundle$cells$cell_id, expected_ids)) {
    espiviz_abort("Bundle cell IDs must be sequential internal IDs.")
  }
  if (!identical(bundle$genes$gene_index, seq_len(nrow(bundle$genes)))) {
    espiviz_abort("Bundle gene_index must match gene order.")
  }
  if (anyDuplicated(toupper(bundle$genes$gene))) {
    espiviz_abort(
      "Bundle gene names are not unique for case-insensitive lookup."
    )
  }
  if (
    any(!is.finite(bundle$cells$umap_1)) ||
      any(!is.finite(bundle$cells$umap_2)) ||
      any(!is.finite(bundle$cells$library_size)) ||
      any(bundle$cells$library_size <= 0)
  ) {
    espiviz_abort("Bundle cell coordinates and library sizes must be finite.")
  }
  if (
    anyNA(bundle$cells$condition) ||
      anyNA(bundle$cells$Mouse) ||
      anyNA(bundle$cells$sample) ||
      any(!nzchar(as.character(bundle$cells$condition))) ||
      any(!nzchar(as.character(bundle$cells$Mouse))) ||
      any(!nzchar(as.character(bundle$cells$sample)))
  ) {
    espiviz_abort(
      "Bundle condition, Mouse, and sample metadata must be complete."
    )
  }
  if (
    any(!is.finite(bundle$counts@x)) ||
      any(bundle$counts@x < 0) ||
      any(abs(bundle$counts@x - round(bundle$counts@x)) > 1e-8)
  ) {
    espiviz_abort("Bundle counts must contain finite, non-negative raw counts.")
  }
  if (
    !isTRUE(all.equal(
      as.numeric(Matrix::rowSums(bundle$counts)),
      bundle$cells$library_size,
      tolerance = 0
    ))
  ) {
    espiviz_abort("Bundle library sizes do not equal raw-count row sums.")
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
      !identical(rownames(normalization$sparse), bundle$genes$gene) ||
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
    espiviz_abort("Bundle scclrR normalization state is invalid.")
  }

  espiviz_require_columns(
    bundle$primary_de,
    c(
      "gene",
      "baseMean",
      "design",
      "contrast",
      "log2FoldChange",
      "lfcSE",
      "pvalue",
      "padj",
      "mean_count_control",
      "mean_count_estim"
    ),
    "Bundle primary_de"
  )
  if (
    !identical(unique(bundle$primary_de$design), contract$primary_de_design)
  ) {
    espiviz_abort("Bundle contains a non-primary DE design.")
  }
  if (
    !identical(unique(bundle$primary_de$contrast), contract$primary_de_contrast)
  ) {
    espiviz_abort("Bundle contains a non-primary DE contrast.")
  }
  if (
    any(grepl("(^|_)paired($|_)", bundle$primary_de$design, ignore.case = TRUE))
  ) {
    espiviz_abort("Bundle contains paired-sensitivity DE.")
  }
  if (any(!bundle$primary_de$gene %in% bundle$genes$gene)) {
    espiviz_abort("Bundle primary_de contains unknown genes.")
  }
  if (anyNA(bundle$primary_de$gene) || anyDuplicated(bundle$primary_de$gene)) {
    espiviz_abort("Bundle primary_de genes must be complete and unique.")
  }
  if (
    !is.numeric(bundle$primary_de$baseMean) ||
      anyNA(bundle$primary_de$baseMean) ||
      any(!is.finite(bundle$primary_de$baseMean)) ||
      any(bundle$primary_de$baseMean < 0)
  ) {
    espiviz_abort("Bundle primary_de baseMean values must be finite and non-negative.")
  }
  if (
    !is.numeric(bundle$primary_de$log2FoldChange) ||
      anyNA(bundle$primary_de$log2FoldChange) ||
      any(!is.finite(bundle$primary_de$log2FoldChange))
  ) {
    espiviz_abort("Bundle primary_de log2FoldChange values must be finite.")
  }
  if (
    !is.numeric(bundle$primary_de$lfcSE) ||
      anyNA(bundle$primary_de$lfcSE) ||
      any(!is.finite(bundle$primary_de$lfcSE)) ||
      any(bundle$primary_de$lfcSE < 0)
  ) {
    espiviz_abort("Bundle primary_de lfcSE values must be finite and non-negative.")
  }
  for (column in c("pvalue", "padj")) {
    value <- bundle$primary_de[[column]]
    if (
      !is.numeric(value) ||
        any(!is.na(value) & (!is.finite(value) | value < 0 | value > 1))
    ) {
      espiviz_abort(paste0(
        "Bundle primary_de ",
        column,
        " values must be missing or finite probabilities from zero to one."
      ))
    }
  }
  for (column in c("mean_count_control", "mean_count_estim")) {
    value <- bundle$primary_de[[column]]
    if (
      !is.numeric(value) ||
        anyNA(value) ||
        any(!is.finite(value)) ||
        any(value < 0)
    ) {
      espiviz_abort(paste0(
        "Bundle primary_de ",
        column,
        " values must be finite and non-negative."
      ))
    }
  }

  if (
    !identical(
      names(bundle$markers),
      c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "rank")
    )
  ) {
    espiviz_abort("Bundle marker schema is invalid.")
  }
  if (any(!bundle$markers$gene %in% bundle$genes$gene)) {
    espiviz_abort("Bundle markers contain unknown genes.")
  }
  if (
    !identical(
      names(bundle$pathways),
      c(
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
    )
  ) {
    espiviz_abort("Bundle pathway schema is invalid.")
  }
  if (!identical(names(bundle$pathway_genes), c("pathway_id", "gene"))) {
    espiviz_abort("Bundle pathway_genes schema is invalid.")
  }
  if (
    nrow(bundle$pathways) == 0L ||
      anyNA(bundle$pathways$pathway_id) ||
      any(!nzchar(as.character(bundle$pathways$pathway_id))) ||
      anyDuplicated(as.character(bundle$pathways$pathway_id)) ||
      anyNA(bundle$pathways$label) ||
      any(!nzchar(as.character(bundle$pathways$label))) ||
      anyDuplicated(as.character(bundle$pathways$label))
  ) {
    espiviz_abort(
      "Bundle pathways must contain complete, unique IDs and labels."
    )
  }
  if (
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
      any(bundle$pathways$gene_count < 0)
  ) {
    espiviz_abort("Bundle pathway statistics are invalid.")
  }
  if (any(!bundle$pathway_genes$pathway_id %in% bundle$pathways$pathway_id)) {
    espiviz_abort("Bundle pathway_genes contain unknown pathways.")
  }
  if (any(!bundle$pathway_genes$gene %in% bundle$genes$gene)) {
    espiviz_abort("Bundle pathway_genes contain unknown genes.")
  }
  if (
    !is.list(bundle$featured_gene_sets) ||
      is.null(names(bundle$featured_gene_sets)) ||
      any(
        !unlist(bundle$featured_gene_sets, use.names = FALSE) %in%
          bundle$genes$gene
      )
  ) {
    espiviz_abort("Bundle featured gene sets are invalid.")
  }

  if (isTRUE(enforce_frozen)) {
    if (
      nrow(bundle$cells) != contract$cells ||
        nrow(bundle$genes) != contract$genes
    ) {
      espiviz_abort(
        "Bundle dimensions do not match the frozen source contract."
      )
    }
    if (length(unique(bundle$cells$cluster)) != length(contract$clusters)) {
      espiviz_abort("Bundle does not contain exactly eight clusters.")
    }
    if (
      !identical(
        sort(unique(as.character(bundle$cells$cluster))),
        contract$clusters
      )
    ) {
      espiviz_abort(
        "Bundle cluster values do not match the frozen source contract."
      )
    }
    if (
      !identical(
        sort(unique(as.character(bundle$cells$condition))),
        sort(contract$conditions)
      )
    ) {
      espiviz_abort(
        "Bundle condition values do not match the frozen source contract."
      )
    }
    if (nrow(bundle$primary_de) != contract$primary_de_rows) {
      espiviz_abort("Bundle does not contain the complete primary DE table.")
    }
  }
  source_hash <- bundle$provenance$source_sha256 %||% ""
  input_source_hash <- bundle$provenance$inputs$source_object$sha256 %||% ""
  if (
    !grepl("^[0-9a-f]{64}$", source_hash) ||
      !identical(source_hash, input_source_hash)
  ) {
    espiviz_abort(
      "Bundle source hash does not match its public input provenance."
    )
  }

  public_characters <- espiviz_collect_characters(bundle[
    names(bundle) != "counts"
  ])
  forbidden_paths <- c("/Users/", "Box-Box", "Library/CloudStorage", "file://")
  for (pattern in forbidden_paths) {
    if (any(grepl(pattern, public_characters, fixed = TRUE))) {
      espiviz_abort("Bundle contains a forbidden local path or Box reference.")
    }
  }
  if (any(vapply(bundle, inherits, logical(1L), what = "Seurat"))) {
    espiviz_abort("Bundle must not contain a Seurat object.")
  }
  invisible(bundle)
}

espiviz_build_bundle <- function(manifest_path, contract = espiviz_contract()) {
  manifest <- espiviz_read_source_manifest(manifest_path)
  input_provenance <- espiviz_input_provenance(manifest)
  observed_source_hash <- input_provenance$source_object$sha256
  if (!identical(observed_source_hash, manifest$source_object$sha256)) {
    espiviz_abort(
      "source_object SHA-256 is '%s'; expected '%s'.",
      observed_source_hash,
      manifest$source_object$sha256
    )
  }

  object <- readRDS(manifest$source_object$path)
  source <- espiviz_extract_source(object, contract = contract)
  rm(object)
  invisible(gc())

  primary_de <- espiviz_read_primary_de(
    manifest$primary_de$path,
    genes = source$genes$gene,
    contract = contract
  )
  markers <- espiviz_read_markers(
    manifest$markers$path,
    genes = source$genes$gene,
    contract = contract
  )
  featured_gene_sets <- espiviz_read_featured_gene_sets(
    manifest$featured_genes_config$path,
    genes = source$genes$gene
  )
  featured_pathways <- espiviz_read_featured_pathways(
    config_path = manifest$featured_pathways_config$path,
    result_paths = lapply(
      manifest[c("ora_up", "ora_down", "gsea")],
      function(item) item$path
    ),
    genes = source$genes$gene
  )
  final_input_provenance <- espiviz_input_provenance(manifest)
  if (!identical(input_provenance, final_input_provenance)) {
    espiviz_abort("An approved input changed while the bundle was being built.")
  }

  bundle <- list(
    schema_version = contract$schema_version,
    data_version = contract$data_version,
    provenance = list(
      source_sha256 = observed_source_hash,
      reduction = contract$reduction,
      cluster_column = contract$cluster_column,
      assay = "RNA",
      count_layer = "counts",
      normalization = paste(
        "scclrR PFlog (target = 'auto', log1p = TRUE, center = TRUE);",
        "dense expression = shifted sparse value - cell center"
      ),
      normalization_check = list(
        method = source$normalization$method,
        target = source$normalization$target,
        log1p = source$normalization$log1p,
        centered = source$normalization$centered,
        alpha = source$normalization$alpha,
        k = source$normalization$k,
        cell_count = nrow(source$cells),
        gene_count = nrow(source$genes),
        center_min = min(source$normalization$center),
        center_max = max(source$normalization$center),
        package_version = source$normalization$package_version,
        package_remote_sha = source$normalization$package_remote_sha
      ),
      dimensions = list(
        genes = nrow(source$genes),
        cells = nrow(source$cells),
        clusters = length(unique(source$cells$cluster)),
        primary_de_rows = nrow(primary_de)
      ),
      inputs = input_provenance,
      software = espiviz_software_versions()
    ),
    palette = espiviz_palette(),
    cells = source$cells,
    genes = source$genes,
    counts = source$counts,
    normalization = source$normalization,
    primary_de = primary_de,
    markers = markers,
    pathways = featured_pathways$pathways,
    pathway_genes = featured_pathways$pathway_genes,
    featured_gene_sets = featured_gene_sets
  )
  espiviz_validate_public_bundle(
    bundle,
    contract = contract,
    enforce_frozen = TRUE
  )
  bundle
}

build_espiviz_bundle <- espiviz_build_bundle
read_source_manifest <- espiviz_read_source_manifest
validate_public_bundle <- espiviz_validate_public_bundle

espiviz_write_bundle <- function(
  bundle,
  path,
  overwrite = FALSE,
  contract = espiviz_contract()
) {
  espiviz_validate_public_bundle(
    bundle,
    contract = contract,
    enforce_frozen = TRUE
  )
  path <- path.expand(path)
  if (file.exists(path) && !isTRUE(overwrite)) {
    espiviz_abort(
      "Output already exists; pass --overwrite to replace it: %s",
      path
    )
  }
  output_directory <- dirname(path)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_directory)) {
    espiviz_abort("Could not create output directory: %s", output_directory)
  }

  temporary <- tempfile(
    pattern = ".espiviz-data-",
    tmpdir = output_directory,
    fileext = ".rds"
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(bundle, temporary, compress = "xz", version = 3L)
  if (!file.rename(temporary, path)) {
    espiviz_abort("Could not atomically move the data bundle to: %s", path)
  }
  list(
    path = normalizePath(path, mustWork = TRUE),
    sha256 = espiviz_sha256(path),
    bytes = unname(file.info(path)$size)
  )
}

write_espiviz_bundle <- espiviz_write_bundle

espiviz_write_public_data_manifest <- function(
  bundle,
  bundle_path,
  manifest_path,
  asset_url,
  contract = espiviz_contract(),
  overwrite = FALSE
) {
  espiviz_require("jsonlite")
  espiviz_validate_public_bundle(
    bundle,
    contract = contract,
    enforce_frozen = TRUE
  )
  if (!grepl("^https://", asset_url)) {
    espiviz_abort("Public data asset URL must use HTTPS.")
  }
  if (!file.exists(bundle_path)) {
    espiviz_abort("Data bundle does not exist: %s", bundle_path)
  }
  if (file.exists(manifest_path) && !isTRUE(overwrite)) {
    espiviz_abort(
      "Public data manifest already exists; pass --overwrite to replace it: %s",
      manifest_path
    )
  }

  input_hashes <- vapply(
    bundle$provenance$inputs,
    function(input) input$sha256,
    character(1L)
  )
  public_manifest <- list(
    schema_version = bundle$schema_version,
    data_version = bundle$data_version,
    asset_url = asset_url,
    asset_sha256 = espiviz_sha256(bundle_path),
    asset_bytes = unname(file.info(bundle_path)$size),
    source_sha256 = bundle$provenance$source_sha256,
    dimensions = bundle$provenance$dimensions,
    input_sha256 = as.list(input_hashes)
  )
  output_directory <- dirname(manifest_path)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  jsonlite::write_json(
    public_manifest,
    path = manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  invisible(public_manifest)
}

write_public_data_manifest <- espiviz_write_public_data_manifest
