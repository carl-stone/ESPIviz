set shell := ["bash", "-euo", "pipefail", "-c"]

source_manifest := "config/source-manifest.local.json"
data_output := "dist/espiviz-data-v1.1.0.rds"

default:
    @just --list

# Validate every frozen input without writing a bundle.
data-dry-run manifest=source_manifest output=data_output:
    Rscript scripts/build-data.R --manifest "{{manifest}}" --output "{{output}}" --dry-run

# Build the immutable public bundle. Pass overwrite=true only intentionally.
data-build manifest=source_manifest output=data_output overwrite="false":
    #!/usr/bin/env bash
    args=(--manifest "{{manifest}}" --output "{{output}}")
    if [[ "{{overwrite}}" == "true" ]]; then
      args+=(--overwrite)
    fi
    Rscript scripts/build-data.R "${args[@]}"

# Run locally. An optional absolute bundle path bypasses the public download.
app-run data_path="":
    #!/usr/bin/env bash
    if [[ -n "{{data_path}}" ]]; then
      export ESPIVIZ_DATA_PATH="{{data_path}}"
    fi
    exec Rscript -e 'shiny::runApp("app", host = "127.0.0.1", launch.browser = TRUE)'

# Run synthetic contract and behavior tests.
test:
    Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)'

# Generate deployment metadata by scanning the deployable app directory only.
manifest:
    Rscript -e 'if (!requireNamespace("rsconnect", quietly = TRUE)) stop("Install rsconnect first"); rsconnect::writeManifest(appDir = "app", appPrimaryDoc = "app.R")'
