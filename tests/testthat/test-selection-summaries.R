comparison_row <- function(summary, gene) {
  summary$comparison[summary$comparison$gene == gene, , drop = FALSE]
}

test_that("single- and multi-cell selections match direct calculations", {
  expect_app_helper("summarize_selection")
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[c(1, 2)]
  expression <- as.numeric(bundle$normalization$sparse["Glul", ]) -
    bundle$normalization$center
  raw_counts <- as.numeric(bundle$counts[, "Glul"])
  selected_index <- bundle$cells$cell_id %in% selected

  observed <- summarize_selection(bundle, selected, "Glul")
  row <- comparison_row(observed, "Glul")

  selected_detected <- raw_counts[selected_index] > 0
  remaining_detected <- raw_counts[!selected_index] > 0
  selected_pct <- mean(selected_detected) * 100
  remaining_pct <- mean(remaining_detected) * 100

  expect_equal(row$selected_n, 2L)
  expect_equal(row$remaining_n, 4L)
  expect_equal(row$selected_mean, mean(expression[selected_index]))
  expect_equal(row$selected_median, stats::median(expression[selected_index]))
  expect_equal(row$selected_detected_n, sum(selected_detected))
  expect_equal(row$selected_detected_pct, selected_pct)
  expect_equal(row$remaining_mean, mean(expression[!selected_index]))
  expect_equal(row$remaining_median, stats::median(expression[!selected_index]))
  expect_equal(row$remaining_detected_n, sum(remaining_detected))
  expect_equal(row$remaining_detected_pct, remaining_pct)
  expect_equal(row$mean_difference, row$selected_mean - row$remaining_mean)
  expect_equal(row$detection_pp_difference, selected_pct - remaining_pct)
  expect_equal(row$detection_ratio, selected_pct / remaining_pct)
  expect_false("mean_ratio" %in% names(observed$comparison))
  expect_false(any(grepl(
    "p[_.]?value|confidence|(^|_)ci($|_)|rank",
    names(observed$comparison),
    ignore.case = TRUE
  )))

  one <- summarize_selection(bundle, selected[[1]], "Glul")
  expect_equal(comparison_row(one, "Glul")$selected_n, 1L)
  expect_equal(comparison_row(one, "Glul")$remaining_n, 5L)
})

test_that("no selection and all-cell selection keep exploration available", {
  expect_app_helper("summarize_selection")
  bundle <- synthetic_bundle()

  no_selection <- summarize_selection(bundle, character(), "Glul")
  all_selection <- summarize_selection(bundle, bundle$cells$cell_id, "Glul")

  for (summary in list(no_selection, all_selection)) {
    row <- comparison_row(summary, "Glul")
    expect_equal(row$selected_n, nrow(bundle$cells))
    expect_equal(row$remaining_n, 0L)
    expect_true(is.na(row$remaining_mean))
    expect_true(is.na(row$mean_difference))
    expect_true(is.na(row$detection_pp_difference))
    expect_true(is.na(row$detection_ratio))
  }
})

test_that("zero-count genes retain PFlog centering and zero detection", {
  expect_app_helper("summarize_selection")
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[1:3]
  expression <- expression_matrix(bundle, "ZeroGene")[, 1L]
  selected_index <- bundle$cells$cell_id %in% selected

  observed <- summarize_selection(bundle, selected, "ZeroGene")
  row <- comparison_row(observed, "ZeroGene")

  expect_equal(row$selected_mean, mean(expression[selected_index]))
  expect_equal(row$remaining_mean, mean(expression[!selected_index]))
  expect_equal(row$selected_detected_n, 0L)
  expect_equal(row$remaining_detected_n, 0L)
  expect_equal(
    row$mean_difference,
    mean(expression[selected_index]) - mean(expression[!selected_index])
  )
  expect_equal(row$detection_pp_difference, 0)
  expect_true(is.na(row$detection_ratio))
})

test_that("selected-cell summaries split by condition and cluster", {
  expect_app_helper("summarize_selection")
  bundle <- synthetic_bundle()
  selected <- bundle$cells$cell_id[c(1, 2, 4, 6)]

  observed <- summarize_selection(bundle, selected, c("Glul", "EGFP"))

  expect_setequal(
    unique(as.character(observed$selected_by_condition$condition)),
    c("p27CKO", "p27CKO +EStim")
  )
  expect_setequal(
    unique(as.character(observed$selected_by_cluster$cluster)),
    c("0", "1", "2")
  )
  expect_setequal(observed$selected_by_condition$gene, c("Glul", "EGFP"))
  expect_setequal(observed$selected_by_cluster$gene, c("Glul", "EGFP"))

  condition_totals <- stats::aggregate(
    cell_count ~ gene,
    data = observed$selected_by_condition,
    FUN = sum
  )
  cluster_totals <- stats::aggregate(
    cell_count ~ gene,
    data = observed$selected_by_cluster,
    FUN = sum
  )
  names(condition_totals)[names(condition_totals) == "cell_count"] <- "total"
  names(cluster_totals)[names(cluster_totals) == "cell_count"] <- "total"
  expect_equal(condition_totals$total, rep(length(selected), 2))
  expect_equal(cluster_totals$total, rep(length(selected), 2))

  control_glul <- observed$selected_by_condition[
    as.character(observed$selected_by_condition$condition) == "p27CKO" &
      observed$selected_by_condition$gene == "Glul",
    ,
    drop = FALSE
  ]
  control_ids <- selected[bundle$cells$condition[match(selected, bundle$cells$cell_id)] == "p27CKO"]
  control_index <- match(control_ids, bundle$cells$cell_id)
  expected_expression <-
    as.numeric(bundle$normalization$sparse["Glul", control_index]) -
    bundle$normalization$center[control_index]
  expect_equal(control_glul$mean_expression, mean(expected_expression))
  expect_equal(control_glul$median_expression, stats::median(expected_expression))
  expect_equal(control_glul$detected_n, sum(bundle$counts[control_index, "Glul"] > 0))
})
