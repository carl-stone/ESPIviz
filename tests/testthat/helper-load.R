repo_root <- normalizePath(
  testthat::test_path("..", ".."),
  winslash = "/",
  mustWork = TRUE
)

app_r_dir <- file.path(repo_root, "app", "R")
app_r_files <- if (dir.exists(app_r_dir)) {
  list.files(app_r_dir, pattern = "[.]R$", full.names = TRUE) |>
    sort()
} else {
  character()
}

for (app_r_file in app_r_files) {
  source(app_r_file, local = TRUE)
}

exporter_file <- file.path(repo_root, "scripts", "lib", "exporter.R")
if (file.exists(exporter_file)) {
  source(exporter_file, local = TRUE)
}

required_app_helpers <- c(
  "compute_pflog_state",
  "pflog_state",
  "normalize_gene_expression",
  "parse_gene_input",
  "paginate_genes",
  "summarize_selection",
  "validate_bundle",
  "load_bundle",
  "selection_expression_export",
  "selection_summary_export",
  "public_bundle_url"
)

expect_app_helper <- function(name) {
  testthat::expect_true(
    exists(name, mode = "function", inherits = TRUE),
    info = paste("Missing app helper:", name)
  )
}
