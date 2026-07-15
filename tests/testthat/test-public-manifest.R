test_that("public data manifest pins one immutable release asset", {
  manifest_path <- file.path(repo_root, "app", "data-manifest.json")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  expected_fields <- c(
    "schema_version",
    "data_version",
    "asset_url",
    "asset_sha256",
    "asset_bytes",
    "source_sha256",
    "dimensions",
    "input_sha256"
  )

  expect_identical(names(manifest), expected_fields)
  expect_identical(manifest$schema_version, "1.0.0")
  expect_identical(manifest$data_version, "1.1.1")
  expect_equal(
    manifest$asset_url,
    paste0(
      "https://github.com/carl-stone/ESPIviz/releases/download/",
      "data-v1.1.1/espiviz-data-v1.1.1.rds"
    )
  )
  expect_match(manifest$asset_sha256, "^[0-9a-f]{64}$")
  expect_false(identical(
    manifest$asset_sha256,
    paste(rep("0", 64), collapse = "")
  ))
  expect_gt(manifest$asset_bytes, 0)
  expect_lt(manifest$asset_bytes, 150 * 1024^2)
  expect_match(manifest$source_sha256, "^[0-9a-f]{64}$")
  expect_identical(
    manifest$source_sha256,
    unname(manifest$input_sha256[["source_object"]])
  )
  expect_equal(
    unname(unlist(manifest$dimensions)),
    c(genes = 38394, cells = 3456, clusters = 8, primary_de_rows = 24601) |>
      unname()
  )
  expect_length(manifest$input_sha256, 7L)
  expect_true(all(grepl(
    "^[0-9a-f]{64}$",
    unname(unlist(manifest$input_sha256))
  )))
})

test_that("local release asset matches the public manifest when present", {
  bundle_path <- file.path(repo_root, "release", "espiviz-data-v1.1.1.rds")
  testthat::skip_if_not(
    file.exists(bundle_path),
    "Local release asset is not committed"
  )
  manifest <- jsonlite::read_json(
    file.path(repo_root, "app", "data-manifest.json"),
    simplifyVector = TRUE
  )

  expect_identical(sha256_file(bundle_path), manifest$asset_sha256)
  expect_equal(unname(file.info(bundle_path)$size), manifest$asset_bytes)

  bundle <- readRDS(bundle_path)
  expect_silent(validate_bundle(bundle))
  expect_identical(nrow(bundle$primary_de), 24601L)
  expect_identical(nrow(bundle$pathways), 11089L)
  expect_true(any(bundle$pathways$p_adjust >= 0.05))
  marker_counts <- table(bundle$markers$cluster)
  expect_true(all(marker_counts <= 25L))
  expect_false(any(grepl(
    "(^|_)paired($|_)",
    bundle$primary_de$design,
    ignore.case = TRUE
  )))
  expect_false(any(
    c(
      "barcode",
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "paired_sensitivity"
    ) %in%
      c(names(bundle), names(bundle$cells))
  ))
})

test_that("Connect manifest contains only app runtime files and dependencies", {
  connect_manifest_path <- file.path(repo_root, "app", "manifest.json")
  expect_true(file.exists(connect_manifest_path))
  manifest <- jsonlite::read_json(connect_manifest_path, simplifyVector = FALSE)
  file_names <- names(manifest$files)
  package_names <- names(manifest$packages)

  expect_identical(manifest$platform, "4.6.0")
  expect_identical(manifest$metadata$appmode, "shiny")
  expect_true(all(
    c(
      "app.R",
      "data-manifest.json",
      "R/00-utils.R",
      "R/90-app-shell.R",
      "www/favicon.svg",
      "www/styles.css"
    ) %in%
      file_names
  ))
  expect_false(any(grepl("^(scripts|config|tests)/", file_names)))
  expect_false(any(
    c(
      "Seurat",
      "SeuratObject",
      "scclrR",
      "AnnotationDbi",
      "org.Mm.eg.db",
      "ESPI"
    ) %in%
      package_names
  ))
  expect_true(all(
    c(
      "shiny",
      "bslib",
      "plotly",
      "DT",
      "ggplot2",
      "jsonlite",
      "digest",
      "Matrix"
    ) %in%
      package_names
  ))

  for (file_name in file_names) {
    file_path <- file.path(repo_root, "app", file_name)
    expect_true(file.exists(file_path), info = file_name)
    expect_identical(
      unname(tools::md5sum(file_path)),
      manifest$files[[file_name]]$checksum,
      info = paste("Stale Connect manifest entry:", file_name)
    )
  }
})
