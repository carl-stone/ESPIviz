#!/usr/bin/env Rscript

espiviz_script_path <- function() {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_argument) > 0L) {
    candidate <- sub("^--file=", "", file_argument[[1L]])
    if (!identical(candidate, "-") && file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  frames <- sys.frames()
  for (index in rev(seq_along(frames))) {
    source_file <- frames[[index]]$ofile
    if (!is.null(source_file)) {
      return(normalizePath(source_file, mustWork = TRUE))
    }
  }
  stop("Could not determine scripts/build-data.R location.", call. = FALSE)
}

source(file.path(dirname(espiviz_script_path()), "lib", "exporter.R"))

espiviz_cli_help <- function() {
  cat(
    paste0(
      "Build the immutable ESPIviz public data bundle.\n\n",
      "Usage:\n",
      "  Rscript scripts/build-data.R --manifest PATH --dry-run\n",
      "  Rscript scripts/build-data.R --manifest PATH --output PATH [--overwrite]\n",
      "    [--public-manifest PATH --asset-url HTTPS_URL]\n\n",
      "Options:\n",
      "  --manifest PATH         Gitignored local source manifest.\n",
      "  --output PATH           Output espiviz-data-v1.0.0.rds path.\n",
      "  --dry-run               Validate and build in memory without writing.\n",
      "  --overwrite             Replace existing output files.\n",
      "  --public-manifest PATH  Write the small deploy-time JSON manifest.\n",
      "  --asset-url URL         Immutable GitHub Release asset URL.\n",
      "  --help                   Show this message.\n"
    )
  )
}

espiviz_parse_cli <- function(arguments) {
  values <- list(
    manifest = NULL,
    output = NULL,
    public_manifest = NULL,
    asset_url = NULL,
    dry_run = FALSE,
    overwrite = FALSE,
    help = FALSE
  )
  value_options <- c(
    "--manifest" = "manifest",
    "--output" = "output",
    "--public-manifest" = "public_manifest",
    "--asset-url" = "asset_url"
  )
  flag_options <- c(
    "--dry-run" = "dry_run",
    "--overwrite" = "overwrite",
    "--help" = "help",
    "-h" = "help"
  )

  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    matched_value <- names(value_options)[startsWith(
      argument,
      paste0(names(value_options), "=")
    )]
    if (length(matched_value) > 0L) {
      option <- matched_value[[1L]]
      values[[value_options[[option]]]] <- sub(
        paste0("^", option, "="),
        "",
        argument
      )
    } else if (argument %in% names(value_options)) {
      if (index == length(arguments)) {
        espiviz_abort("Option %s requires a value.", argument)
      }
      index <- index + 1L
      values[[value_options[[argument]]]] <- arguments[[index]]
    } else if (argument %in% names(flag_options)) {
      values[[flag_options[[argument]]]] <- TRUE
    } else {
      espiviz_abort("Unknown option: %s", argument)
    }
    index <- index + 1L
  }
  values
}

espiviz_run_cli <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  options <- espiviz_parse_cli(arguments)
  if (isTRUE(options$help)) {
    espiviz_cli_help()
    return(invisible(NULL))
  }
  if (is.null(options$manifest)) {
    espiviz_abort("--manifest is required.")
  }
  if (!isTRUE(options$dry_run) && is.null(options$output)) {
    espiviz_abort("--output is required unless --dry-run is used.")
  }
  if (xor(is.null(options$public_manifest), is.null(options$asset_url))) {
    espiviz_abort(
      "--public-manifest and --asset-url must be supplied together."
    )
  }
  if (isTRUE(options$dry_run) && !is.null(options$public_manifest)) {
    espiviz_abort("--public-manifest cannot be used with --dry-run.")
  }

  message("Validating frozen source and approved results...")
  bundle <- build_espiviz_bundle(options$manifest)
  message(
    sprintf(
      paste0(
        "Validated %s genes, %s cells, %s clusters, %s primary DE rows, ",
        "and %s featured pathways."
      ),
      format(nrow(bundle$genes), big.mark = ","),
      format(nrow(bundle$cells), big.mark = ","),
      length(unique(bundle$cells$cluster)),
      format(nrow(bundle$primary_de), big.mark = ","),
      nrow(bundle$pathways)
    )
  )
  message(
    sprintf(
      "Normalization sampled max absolute error: %.3g",
      bundle$provenance$normalization_check$max_abs_error
    )
  )
  if (isTRUE(options$dry_run)) {
    message("Dry run complete; no files were written.")
    return(invisible(bundle))
  }

  written <- write_espiviz_bundle(
    bundle,
    options$output,
    overwrite = options$overwrite
  )
  message(sprintf("Wrote %s", written$path))
  message(sprintf("SHA-256 %s", written$sha256))
  message(sprintf("Compressed size %.1f MB", written$bytes / 1024^2))

  if (!is.null(options$public_manifest)) {
    write_public_data_manifest(
      bundle = bundle,
      bundle_path = written$path,
      manifest_path = options$public_manifest,
      asset_url = options$asset_url,
      overwrite = options$overwrite
    )
    message(sprintf("Wrote %s", normalizePath(options$public_manifest)))
  }
  invisible(written)
}

if (sys.nframe() == 0L) {
  tryCatch(
    espiviz_run_cli(),
    error = function(error) {
      message("ERROR: ", conditionMessage(error))
      quit(status = 1L, save = "no")
    }
  )
}
