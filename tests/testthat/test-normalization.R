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

  expect_error(
    normalize_gene_expression(1, 0),
    "library|positive",
    ignore.case = TRUE
  )
  expect_error(
    normalize_gene_expression(c(1, 2), 100),
    "length",
    ignore.case = TRUE
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
    "observed <- expression_matrix(bundle, 'Glul')",
    "stopifnot(identical(dim(observed), c(6L, 1L)))"
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
