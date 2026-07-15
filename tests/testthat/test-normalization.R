test_that("PFlog applies the shifted centered log-ratio formula", {
  expect_app_helper("normalize_gene_expression")

  counts <- c(10, 0, 25, 100)
  alpha <- 0.125
  center <- c(0.2, 0.4, 0.6, 0.8)

  observed <- normalize_gene_expression(counts, alpha, center)
  expected <- log1p(4 * alpha * counts) - center

  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(observed[[2]], -center[[2]])
})

test_that("PFlog rejects invalid alpha and cell centers", {
  expect_app_helper("normalize_gene_expression")

  expect_error(
    normalize_gene_expression(1, 0, 0),
    "alpha|positive",
    ignore.case = TRUE
  )
  expect_error(
    normalize_gene_expression(c(1, 2), 0.1, 0),
    "center|length",
    ignore.case = TRUE
  )
})

test_that("PFlog state and expression are numerically identical to scclrR", {
  testthat::skip_if_not_installed("scclrR")
  expect_app_helper("compute_pflog_state")
  expect_app_helper("normalize_gene_expression")

  bundle <- synthetic_bundle()
  state <- compute_pflog_state(bundle$counts)
  reference <- scclrR::normalize_matrix(
    Matrix::t(bundle$counts),
    target = "auto"
  )
  reference_dense <- t(as.matrix(reference$sparse))
  reference_dense <- sweep(reference_dense, 1L, reference$center, "-")
  observed <- normalize_gene_expression(
    bundle$counts,
    state$alpha,
    state$center
  )

  expect_equal(state$alpha, reference$alpha, tolerance = 1e-12)
  expect_equal(state$scale, 4 * reference$alpha, tolerance = 1e-12)
  expect_equal(state$k, reference$k, tolerance = 1e-12)
  expect_equal(unname(state$center), reference$center, tolerance = 1e-12)
  expect_equal(observed, reference_dense, tolerance = 1e-12)

  zero_index <- as.matrix(bundle$counts) == 0
  expected_zeros <- -state$center[row(zero_index)[zero_index]]
  expect_equal(
    unname(observed[zero_index]),
    unname(expected_zeros),
    tolerance = 1e-12
  )
})

test_that("sparse expression lookup loads its matrix runtime in a fresh process", {
  bundle_path <- tempfile(fileext = ".rds")
  script_path <- tempfile(fileext = ".R")
  saveRDS(synthetic_bundle(), bundle_path)

  script <- c(
    sprintf("repo <- %s", deparse(repo_root)),
    sprintf("bundle_path <- %s", deparse(bundle_path)),
    "app_files <- sort(list.files(file.path(repo, 'app', 'R'), pattern = '[.]R$', full.names = TRUE))",
    "for (app_file in app_files) source(app_file)",
    "bundle <- readRDS(bundle_path)",
    "stopifnot(!'Matrix' %in% loadedNamespaces())",
    "stopifnot(!'scclrR' %in% loadedNamespaces())",
    "observed <- expression_matrix(bundle, 'Glul')",
    "stopifnot(identical(dim(observed), c(6L, 1L)))",
    "stopifnot(!'scclrR' %in% loadedNamespaces())"
  )
  writeLines(script, script_path)

  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", script_path),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  expect_identical(status, 0L, info = paste(output, collapse = "\n"))
})
