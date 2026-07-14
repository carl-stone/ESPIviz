test_that("normalization reproduces log1p counts per ten thousand", {
  expect_app_helper("normalize_gene_expression")

  counts <- c(10, 0, 25, 100)
  library_size <- c(100, 100, 250, 1000)

  observed <- normalize_gene_expression(counts, library_size)
  expected <- log1p(10000 * counts / library_size)

  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(observed[[2]], 0)
})

test_that("normalization rejects invalid library sizes", {
  expect_app_helper("normalize_gene_expression")

  expect_error(normalize_gene_expression(1, 0), "library|positive", ignore.case = TRUE)
  expect_error(normalize_gene_expression(c(1, 2), 100), "length", ignore.case = TRUE)
})
