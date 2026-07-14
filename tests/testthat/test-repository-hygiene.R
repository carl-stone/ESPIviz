test_that("repository contains no committed manuscript objects", {
  files <- system2(
    "git",
    c("-C", shQuote(repo_root), "ls-files"),
    stdout = TRUE,
    stderr = FALSE
  )
  if (length(files) == 0L) {
    files <- list.files(
      repo_root,
      recursive = TRUE,
      full.names = FALSE,
      all.files = TRUE,
      no.. = TRUE
    )
    files <- files[!grepl(
      "^([.]git|data|dist|release|app/[.]data-cache)(/|$)",
      files
    )]
  }

  expect_false(any(grepl("[.](rds|Rds|rda|RData)$", files)))
})

test_that("production code has no private repository or cloud-path dependency", {
  public_metadata <- c(
    list.files(file.path(repo_root, "config"), recursive = TRUE, full.names = TRUE),
    file.path(repo_root, "app", "data-manifest.json")
  )
  public_metadata <- public_metadata[file.exists(public_metadata)]
  public_metadata <- public_metadata[
    basename(public_metadata) != "source-manifest.local.json"
  ]
  metadata_text <- paste(
    unlist(lapply(public_metadata, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  r_files <- c(
    list.files(file.path(repo_root, "app"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
    list.files(file.path(repo_root, "scripts"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  )
  r_text <- paste(
    unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )

  expect_false(grepl("/Users/", metadata_text, fixed = TRUE))
  expect_false(grepl("Box-Box", metadata_text, fixed = TRUE))
  expect_false(grepl("(path|file)\\s*<-\\s*['\"](/Users/|~/Library/CloudStorage)", r_text, perl = TRUE))
  expect_false(grepl("devtools::load_all", r_text, fixed = TRUE))
  expect_false(grepl("library\\s*\\(\\s*ESPI\\s*\\)", r_text, perl = TRUE))
  expect_false(grepl("ESPI::", r_text, fixed = TRUE))
})

test_that("runtime dependencies exclude Seurat and the analysis package", {
  description <- read.dcf(file.path(repo_root, "DESCRIPTION"))
  imports <- description[[1, "Imports"]]

  expect_false(grepl("Seurat", imports, fixed = TRUE))
  expect_false(grepl("ESPI", imports, fixed = TRUE))
})
