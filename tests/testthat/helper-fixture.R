synthetic_bundle <- function(extra_genes = 0L) {
  core_genes <- c("Glul", "EGFP", "ZeroGene", "Mcm2", "Other")
  additional_genes <- if (extra_genes > 0L) {
    sprintf("Gene%03d", seq_len(extra_genes))
  } else {
    character()
  }
  genes <- c(core_genes, additional_genes)

  core_counts <- rbind(
    c(10, 0, 0, 0, 90),
    c(0, 20, 0, 0, 80),
    c(5, 5, 0, 10, 80),
    c(0, 0, 0, 50, 50),
    c(100, 0, 0, 0, 0),
    c(0, 1, 0, 9, 90)
  )
  if (length(additional_genes) > 0L) {
    extra_counts <- matrix(
      rep(seq_along(additional_genes) %% 4L, each = nrow(core_counts)),
      nrow = nrow(core_counts)
    )
    raw_counts <- cbind(core_counts, extra_counts)
  } else {
    raw_counts <- core_counts
  }

  cell_ids <- paste0("cell_", seq_len(nrow(raw_counts)))
  dimnames(raw_counts) <- list(cell_ids, genes)
  counts <- methods::as(Matrix::Matrix(raw_counts, sparse = TRUE), "dgCMatrix")

  cells <- data.frame(
    cell_id = cell_ids,
    umap_1 = c(-2, -1, 0, 1, 2, 0.5),
    umap_2 = c(0, 1, -1, 1, 0, -0.5),
    cluster = factor(c("0", "0", "1", "1", "2", "2")),
    condition = factor(
      c("p27CKO", "p27CKO", "p27CKO", rep("p27CKO +EStim", 3)),
      levels = c("p27CKO", "p27CKO +EStim")
    ),
    Mouse = c("10", "10", "33", "3", "3", "30"),
    sample = c("10_control", "10_control", "33_control", "3_estim", "3_estim", "30_estim"),
    library_size = as.numeric(Matrix::rowSums(counts)),
    stringsAsFactors = FALSE
  )

  primary_de <- data.frame(
    gene = c("Glul", "EGFP", "Mcm2", "Other"),
    baseMean = c(100, 20, 75, 500),
    log2FoldChange = c(1.2, 0.3, -1.5, 0.1),
    lfcSE = c(0.2, 0.4, 0.3, 0.2),
    stat = c(6, 0.75, -5, 0.5),
    pvalue = c(1e-6, 0.4, 3e-5, 0.6),
    padj = c(4e-5, 0.5, 4e-4, 0.7),
    mean_count_control = c(55, 18, 120, 490),
    mean_count_estim = c(130, 22, 42, 510),
    contrast = rep("estim_vs_control", 4),
    design = rep("primary_unpaired_condition", 4),
    stringsAsFactors = FALSE
  )

  markers <- data.frame(
    cluster = c("0", "0", "1", "1", "2", "2"),
    gene = c("Glul", "EGFP", "Mcm2", "Other", "Glul", "Mcm2"),
    avg_log2FC = c(1.8, 1.1, 2.2, 1.0, 2.4, 1.3),
    pct.1 = c(0.8, 0.5, 0.9, 0.6, 0.95, 0.7),
    pct.2 = c(0.2, 0.1, 0.25, 0.3, 0.3, 0.2),
    p_val_adj = c(1e-8, 1e-4, 1e-9, 1e-3, 1e-10, 1e-5),
    rank = rep(1:2, 3),
    stringsAsFactors = FALSE
  )

  pathways <- data.frame(
    pathway_id = c("oxidative_stress", "dna_replication"),
    label = c("Oxidative stress", "DNA replication"),
    source = c("GSEA", "GSEA"),
    direction = c("E-Stim", "Control"),
    description = c("Oxidative stress", "DNA replication"),
    p_value = c(1e-4, 2e-5),
    p_adjust = c(1e-3, 3e-4),
    score = c(1.8, -2.1),
    gene_count = c(2L, 2L),
    stringsAsFactors = FALSE
  )
  pathway_genes <- data.frame(
    pathway_id = c(
      "oxidative_stress",
      "oxidative_stress",
      "dna_replication",
      "dna_replication"
    ),
    gene = c("Glul", "EGFP", "Mcm2", "Other"),
    stringsAsFactors = FALSE
  )

  list(
    schema_version = "1.0.0",
    data_version = "1.0.0",
    provenance = list(
      source_sha256 = paste(rep("a", 64), collapse = ""),
      reduction = "umap_pflog_mg_selected_no_filter_cc_dims20",
      cluster_column = "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5",
      inputs = list(
        source_object = list(sha256 = paste(rep("a", 64), collapse = "")),
        primary_de = list(sha256 = paste(rep("b", 64), collapse = ""))
      )
    ),
    palette = list(
      expression = c(low = "#2166AC", mid = "#BDBDBD", high = "#E31A8C"),
      condition = c("p27CKO" = "#2166AC", "p27CKO +EStim" = "#E31A8C")
    ),
    cells = cells,
    genes = data.frame(
      gene = genes,
      gene_index = seq_along(genes),
      stringsAsFactors = FALSE
    ),
    counts = counts,
    primary_de = primary_de,
    markers = markers,
    pathways = pathways,
    pathway_genes = pathway_genes,
    featured_gene_sets = list(
      study_highlights = c("Glul", "EGFP", "Mcm2"),
      all_synthetic = genes
    )
  )
}

synthetic_expectations <- function(bundle = synthetic_bundle()) {
  list(
    schema_version = "1.0.0",
    cells = nrow(bundle$cells),
    genes = nrow(bundle$genes),
    clusters = length(unique(as.character(bundle$cells$cluster))),
    de_rows = nrow(bundle$primary_de)
  )
}
