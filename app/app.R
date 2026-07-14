app_directory <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(app_directory, "R", "00-utils.R"))) {
  source_path <- tryCatch(sys.frame(1)$ofile, error = function(error) NULL)
  if (is.null(source_path)) source_path <- "app.R"
  app_directory <- normalizePath(dirname(source_path), mustWork = TRUE)
}

source_files <- sort(list.files(file.path(app_directory, "R"), pattern = "[.]R$", full.names = TRUE))
invisible(lapply(source_files, sys.source, envir = environment()))

bundle_result <- tryCatch(
  list(
    bundle = load_bundle_once(
      manifest_path = file.path(app_directory, "data-manifest.json")
    ),
    error = NULL
  ),
  error = function(error) {
    message("ESPIviz data load failed: ", conditionMessage(error))
    list(bundle = NULL, error = error)
  }
)

if (is.null(bundle_result$error)) {
  shiny::shinyApp(
    ui = app_ui(bundle_result$bundle),
    server = app_server(bundle_result$bundle)
  )
} else {
  shiny::shinyApp(ui = maintenance_ui(), server = function(input, output, session) {})
}
