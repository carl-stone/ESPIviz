test_that("runtime expression only applies scclrR's exported centers", {
  expect_app_helper("center_shifted_expression")

  shifted <- c(1.2, 0, 2.5, 4.1)
  center <- c(0.2, 0.4, 0.6, 0.8)

  observed <- center_shifted_expression(shifted, center)
  expected <- shifted - center

  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(observed[[2]], -center[[2]])
})

test_that("runtime expression rejects invalid exported centers", {
  expect_app_helper("center_shifted_expression")

  expect_error(
    center_shifted_expression(c(1, 2), 0),
    "center|length",
    ignore.case = TRUE
  )
})

test_that("the synthetic bundle stores values produced by scclrR", {
  testthat::skip_if_not_installed("scclrR")
  expect_app_helper("normalization_state")

  bundle <- synthetic_bundle()
  state <- normalization_state(bundle)
  reference <- scclrR::normalize_matrix(
    Matrix::t(bundle$counts),
    target = "auto",
    alpha = 0.125
  )

  expect_equal(state$sparse, reference$sparse, tolerance = 1e-12)
  expect_equal(unname(state$center), reference$center, tolerance = 1e-12)
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
