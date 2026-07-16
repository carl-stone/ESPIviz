normalization_state <- function(bundle) {
  state <- bundle$normalization
  if (
    !is.list(state) ||
      !inherits(state$sparse, "sparseMatrix") ||
      !is.numeric(state$alpha) ||
      length(state$alpha) != 1L ||
      !is.finite(state$alpha) ||
      state$alpha <= 0 ||
      !is.numeric(state$center) ||
      length(state$center) != nrow(bundle$cells) ||
      any(!is.finite(state$center)) ||
      !identical(dim(state$sparse), rev(dim(bundle$counts)))
  ) {
    stop("The exported scclrR normalization state is invalid.", call. = FALSE)
  }
  state
}

center_shifted_expression <- function(shifted, center) {
  center <- as.numeric(center)
  values <- if (inherits(shifted, "sparseMatrix")) {
    shifted@x
  } else {
    as.numeric(shifted)
  }
  if (any(!is.finite(values)) || any(!is.finite(center))) {
    stop("Exported scclrR values and centers must be finite.", call. = FALSE)
  }
  if (is.null(dim(shifted))) {
    if (length(shifted) != length(center)) {
      stop(
        "Shifted values and scclrR centers must have the same length.",
        call. = FALSE
      )
    }
    return(as.numeric(shifted) - center)
  }
  if (nrow(shifted) != length(center)) {
    stop("Shifted-matrix rows must match the scclrR center vector.", call. = FALSE)
  }
  sweep(as.matrix(shifted), 1L, center, "-")
}

parse_gene_input <- function(text, universe, max_genes = Inf) {
  universe <- as.character(universe)
  tokens <- unlist(
    strsplit(
      paste(as.character(text %||% ""), collapse = "\n"),
      "[,;[:space:]]+"
    ),
    use.names = FALSE
  )
  tokens <- sub("^[\"']", "", tokens)
  tokens <- sub("[\"']$", "", tokens)
  tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
  tokens <- tokens[!duplicated(casefold_key(tokens))]
  universe_keys <- casefold_key(universe)
  matched <- match(casefold_key(tokens), universe_keys)
  genes <- universe[matched[!is.na(matched)]]
  missing <- tokens[is.na(matched)]

  truncated <- character()
  if (is.finite(max_genes) && length(genes) > as.integer(max_genes)) {
    truncated <- genes[(as.integer(max_genes) + 1L):length(genes)]
    genes <- genes[seq_len(as.integer(max_genes))]
  }
  list(genes = genes, missing = missing, truncated = truncated)
}

read_gene_upload <- function(upload) {
  if (is.null(upload) || !file.exists(upload$datapath)) {
    return(character())
  }
  extension <- tolower(tools::file_ext(upload$name %||% ""))
  values <- tryCatch(
    {
      if (identical(extension, "csv")) {
        unlist(
          utils::read.csv(
            upload$datapath,
            header = FALSE,
            stringsAsFactors = FALSE,
            check.names = FALSE
          ),
          use.names = FALSE
        )
      } else if (identical(extension, "tsv")) {
        unlist(
          utils::read.delim(
            upload$datapath,
            header = FALSE,
            stringsAsFactors = FALSE,
            check.names = FALSE
          ),
          use.names = FALSE
        )
      } else {
        readLines(upload$datapath, warn = FALSE)
      }
    },
    error = function(error) readLines(upload$datapath, warn = FALSE)
  )
  values <- as.character(values)
  values <- sub("^\\ufeff", "", values)
  header_names <- c("gene", "genes", "symbol", "gene_symbol", "gene symbol")
  values[!casefold_key(trimws(values)) %in% header_names]
}

paginate_genes <- function(genes, page, page_size = 25L) {
  genes <- as.character(genes)
  page_size <- max(1L, as.integer(page_size))
  total <- length(genes)
  pages <- max(1L, ceiling(total / page_size))
  page <- min(max(1L, as.integer(page %||% 1L)), pages)
  start <- if (total == 0L) 0L else (page - 1L) * page_size + 1L
  end <- if (total == 0L) 0L else min(total, page * page_size)
  list(
    genes = if (total == 0L) character() else genes[start:end],
    page = page,
    pages = pages,
    total = total,
    start = start,
    end = end
  )
}

match_bundle_genes <- function(bundle, genes) {
  parsed <- parse_gene_input(genes, bundle_gene_names(bundle))
  parsed
}

expression_matrix <- function(bundle, genes, cell_ids = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "The Matrix package is required for sparse expression lookup.",
      call. = FALSE
    )
  }
  parsed <- match_bundle_genes(bundle, genes)
  if (is.null(cell_ids)) {
    cell_index <- seq_len(nrow(bundle$cells))
  } else {
    cell_index <- match(cell_ids, bundle$cells$cell_id)
    cell_index <- cell_index[!is.na(cell_index)]
  }
  gene_index <- match(parsed$genes, bundle_gene_names(bundle))
  if (length(gene_index) == 0L) {
    result <- matrix(numeric(), nrow = length(cell_index), ncol = 0L)
    rownames(result) <- bundle$cells$cell_id[cell_index]
    attr(result, "missing_genes") <- parsed$missing
    return(result)
  }
  state <- normalization_state(bundle)
  shifted <- Matrix::t(state$sparse[gene_index, cell_index, drop = FALSE])
  result <- center_shifted_expression(
    shifted,
    state$center[cell_index]
  )
  rownames(result) <- bundle$cells$cell_id[cell_index]
  colnames(result) <- parsed$genes
  attr(result, "missing_genes") <- parsed$missing
  result
}

group_expression_summary <- function(expression, counts, groups, group_name) {
  if (nrow(expression) == 0L || ncol(expression) == 0L) {
    empty <- data.frame(
      gene = character(),
      group = character(),
      cell_count = integer(),
      mean_expression = numeric(),
      median_expression = numeric(),
      detected_n = integer(),
      detected_pct = numeric(),
      stringsAsFactors = FALSE
    )
    names(empty)[[2L]] <- group_name
    return(empty)
  }
  groups <- as.character(groups)
  group_levels <- unique(groups)
  rows <- vector("list", length(group_levels) * ncol(expression))
  cursor <- 0L
  for (group in group_levels) {
    index <- which(groups == group)
    for (gene_index in seq_len(ncol(expression))) {
      cursor <- cursor + 1L
      values <- expression[index, gene_index]
      raw <- counts[index, gene_index]
      rows[[cursor]] <- data.frame(
        gene = colnames(expression)[[gene_index]],
        group = group,
        cell_count = length(index),
        mean_expression = mean(values),
        median_expression = stats::median(values),
        detected_n = sum(raw > 0),
        detected_pct = 100 * mean(raw > 0),
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  names(result)[[2L]] <- group_name
  rownames(result) <- NULL
  result
}

summarize_selection <- function(
  bundle,
  selected_cell_ids,
  genes,
  include_splits = TRUE,
  chunk_size = 128L
) {
  state <- normalization_state(bundle)
  parsed <- match_bundle_genes(bundle, genes)
  all_ids <- as.character(bundle$cells$cell_id)
  selected_cell_ids <- intersect(
    all_ids,
    as.character(selected_cell_ids %||% character())
  )
  if (length(selected_cell_ids) == 0L) {
    selected_cell_ids <- all_ids
  }
  selected_index <- which(all_ids %in% selected_cell_ids)
  remaining_index <- setdiff(seq_along(all_ids), selected_index)
  gene_index <- match(parsed$genes, bundle_gene_names(bundle))

  comparison <- data.frame(
    gene = parsed$genes,
    selected_n = integer(length(gene_index)),
    remaining_n = integer(length(gene_index)),
    selected_mean = numeric(length(gene_index)),
    selected_median = numeric(length(gene_index)),
    selected_detected_n = integer(length(gene_index)),
    selected_detected_pct = numeric(length(gene_index)),
    remaining_mean = numeric(length(gene_index)),
    remaining_median = numeric(length(gene_index)),
    remaining_detected_n = integer(length(gene_index)),
    remaining_detected_pct = numeric(length(gene_index)),
    mean_difference = numeric(length(gene_index)),
    detection_pp_difference = numeric(length(gene_index)),
    detection_ratio = numeric(length(gene_index)),
    stringsAsFactors = FALSE
  )

  by_condition_chunks <- list()
  by_cluster_chunks <- list()
  by_sample_chunks <- list()
  if (length(gene_index) > 0L) {
    chunk_size <- max(1L, as.integer(chunk_size))
    chunks <- split(
      seq_along(gene_index),
      ceiling(seq_along(gene_index) / chunk_size)
    )
    for (chunk_number in seq_along(chunks)) {
      columns <- chunks[[chunk_number]]
      all_counts <- as.matrix(bundle$counts[,
        gene_index[columns],
        drop = FALSE
      ])
      all_shifted <- Matrix::t(state$sparse[
        gene_index[columns],
        ,
        drop = FALSE
      ])
      all_expression <- center_shifted_expression(
        all_shifted,
        state$center
      )
      colnames(all_counts) <- parsed$genes[columns]
      colnames(all_expression) <- parsed$genes[columns]

      summarize_rows <- function(index) {
        if (length(index) == 0L) {
          return(list(
            n = 0L,
            mean = rep(NA_real_, length(columns)),
            median = rep(NA_real_, length(columns)),
            detected_n = integer(length(columns)),
            detected_pct = rep(NA_real_, length(columns))
          ))
        }
        expression <- all_expression[index, , drop = FALSE]
        raw <- all_counts[index, , drop = FALSE]
        list(
          n = length(index),
          mean = colMeans(expression),
          median = apply(expression, 2L, stats::median),
          detected_n = colSums(raw > 0),
          detected_pct = 100 * colMeans(raw > 0)
        )
      }

      selected <- summarize_rows(selected_index)
      remaining <- summarize_rows(remaining_index)
      comparison$selected_n[columns] <- selected$n
      comparison$remaining_n[columns] <- remaining$n
      comparison$selected_mean[columns] <- selected$mean
      comparison$selected_median[columns] <- selected$median
      comparison$selected_detected_n[columns] <- selected$detected_n
      comparison$selected_detected_pct[columns] <- selected$detected_pct
      comparison$remaining_mean[columns] <- remaining$mean
      comparison$remaining_median[columns] <- remaining$median
      comparison$remaining_detected_n[columns] <- remaining$detected_n
      comparison$remaining_detected_pct[columns] <- remaining$detected_pct

      if (isTRUE(include_splits)) {
        selected_cells <- bundle$cells[selected_index, , drop = FALSE]
        selected_expression <- all_expression[selected_index, , drop = FALSE]
        selected_counts <- all_counts[selected_index, , drop = FALSE]
        by_condition_chunks[[chunk_number]] <- group_expression_summary(
          selected_expression,
          selected_counts,
          selected_cells$condition,
          "condition"
        )
        by_cluster_chunks[[chunk_number]] <- group_expression_summary(
          selected_expression,
          selected_counts,
          selected_cells$cluster,
          "cluster"
        )
        by_sample_chunks[[chunk_number]] <- group_expression_summary(
          selected_expression,
          selected_counts,
          selected_cells$sample,
          "sample"
        )
      }
    }
    comparison$mean_difference <- comparison$selected_mean -
      comparison$remaining_mean
    comparison$detection_pp_difference <-
      comparison$selected_detected_pct - comparison$remaining_detected_pct
    comparison$detection_ratio <- safe_ratio(
      comparison$selected_detected_pct,
      comparison$remaining_detected_pct
    )
  }

  if (length(by_condition_chunks) > 0L) {
    by_condition <- do.call(rbind, by_condition_chunks)
    by_cluster <- do.call(rbind, by_cluster_chunks)
    by_sample <- do.call(rbind, by_sample_chunks)
    rownames(by_condition) <- NULL
    rownames(by_cluster) <- NULL
    rownames(by_sample) <- NULL

    sample_lookup <- unique(bundle$cells[c("sample", "condition")])
    by_sample$condition <- as.character(sample_lookup$condition[
      match(as.character(by_sample$sample), as.character(sample_lookup$sample))
    ])
    sample_levels <- unique(as.character(bundle$cells$sample))
    by_sample$sample <- factor(
      as.character(by_sample$sample),
      levels = sample_levels
    )
  } else {
    by_condition <- group_expression_summary(
      matrix(numeric(), nrow = 0L, ncol = 0L),
      matrix(numeric(), nrow = 0L, ncol = 0L),
      character(),
      "condition"
    )
    by_cluster <- group_expression_summary(
      matrix(numeric(), nrow = 0L, ncol = 0L),
      matrix(numeric(), nrow = 0L, ncol = 0L),
      character(),
      "cluster"
    )
    by_sample <- group_expression_summary(
      matrix(numeric(), nrow = 0L, ncol = 0L),
      matrix(numeric(), nrow = 0L, ncol = 0L),
      character(),
      "sample"
    )
    by_sample$condition <- character()
  }

  list(
    comparison = comparison,
    selected_by_condition = by_condition,
    selected_by_cluster = by_cluster,
    selected_by_sample = by_sample,
    selected_cell_ids = selected_cell_ids,
    remaining_cell_ids = all_ids[remaining_index],
    missing_genes = parsed$missing
  )
}

selection_summary_export <- function(bundle, selected_cell_ids, genes) {
  summarize_selection(
    bundle,
    selected_cell_ids,
    genes,
    include_splits = FALSE
  )$comparison
}

selection_expression_export <- function(bundle, selected_cell_ids, genes) {
  state <- normalization_state(bundle)
  parsed <- match_bundle_genes(bundle, genes)
  all_ids <- as.character(bundle$cells$cell_id)
  selected_cell_ids <- intersect(
    all_ids,
    as.character(selected_cell_ids %||% character())
  )
  if (length(selected_cell_ids) == 0L) {
    selected_cell_ids <- all_ids
  }
  cell_index <- match(selected_cell_ids, all_ids)
  gene_index <- match(parsed$genes, bundle_gene_names(bundle))
  if (length(gene_index) == 0L || length(cell_index) == 0L) {
    return(data.frame(
      cell_id = character(),
      gene = character(),
      count = numeric(),
      normalized_expression = numeric(),
      condition = character(),
      cluster = character(),
      Mouse = character(),
      sample = character(),
      stringsAsFactors = FALSE
    ))
  }
  counts <- as.matrix(bundle$counts[cell_index, gene_index, drop = FALSE])
  shifted <- Matrix::t(state$sparse[gene_index, cell_index, drop = FALSE])
  normalized <- center_shifted_expression(
    shifted,
    state$center[cell_index]
  )
  rows <- rep(seq_along(cell_index), times = length(gene_index))
  columns <- rep(seq_along(gene_index), each = length(cell_index))
  metadata <- bundle$cells[cell_index[rows], , drop = FALSE]
  result <- data.frame(
    cell_id = metadata$cell_id,
    gene = parsed$genes[columns],
    count = as.vector(counts),
    normalized_expression = as.vector(normalized),
    condition = metadata$condition,
    cluster = metadata$cluster,
    Mouse = metadata$Mouse,
    sample = metadata$sample,
    stringsAsFactors = FALSE
  )
  attr(result, "missing_genes") <- parsed$missing
  result
}

write_selection_expression_export <- function(
  bundle,
  selected_cell_ids,
  genes,
  file
) {
  state <- normalization_state(bundle)
  parsed <- match_bundle_genes(bundle, genes)
  all_ids <- as.character(bundle$cells$cell_id)
  selected_cell_ids <- intersect(
    all_ids,
    as.character(selected_cell_ids %||% character())
  )
  if (length(selected_cell_ids) == 0L) {
    selected_cell_ids <- all_ids
  }
  cell_index <- match(selected_cell_ids, all_ids)
  gene_index <- match(parsed$genes, bundle_gene_names(bundle))
  if (
    length(gene_index) == ncol(bundle$counts) &&
      setequal(gene_index, seq_len(ncol(bundle$counts)))
  ) {
    gene_index <- seq_len(ncol(bundle$counts))
    parsed$genes <- bundle_gene_names(bundle)
  }
  shifted <- methods::as(
    state$sparse[gene_index, cell_index, drop = FALSE],
    "dgCMatrix"
  )
  centers <- state$center[cell_index]
  names(centers) <- bundle$cells$cell_id[cell_index]
  export <- list(
    schema_version = ESPIVIZ_SCHEMA_VERSION,
    normalization = paste(
      "scclrR PFlog (target = 'auto', log1p = TRUE, center = TRUE);",
      "dense expression = shifted sparse value - cell center"
    ),
    cells = bundle$cells[cell_index, , drop = FALSE],
    genes = data.frame(
      gene = parsed$genes,
      gene_index = seq_along(parsed$genes),
      stringsAsFactors = FALSE
    ),
    normalized_expression = list(
      sparse = shifted,
      center = centers,
      k = state$k,
      alpha = state$alpha
    ),
    raw_counts = "Available in the complete ESPIviz public data bundle.",
    missing_genes = parsed$missing
  )
  saveRDS(export, file, compress = "gzip", version = 3L)
  invisible(file)
}
