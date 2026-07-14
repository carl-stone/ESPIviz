# ESPIviz agent guide

ESPIviz is a standalone public Shiny companion to the ESPI manuscript. Keep it
small, reviewable, and independent. Never add private repository history,
private repository links, local cloud-storage paths, original cell barcodes, or
analysis-development artifacts.

## Repository boundaries

- `app/` contains the deployable application and production runtime code.
- `scripts/` contains the standalone data exporter.
- `config/` contains human-readable public configuration and a source-manifest
  example. The real source manifest stays gitignored.
- `tests/` contains synthetic data only.
- Generated bundles belong outside Git. Commit only their public manifest.

The app loads one immutable, checksummed bundle. It must not depend on Seurat or
the private analysis package at runtime. The exporter may use SeuratObject, but
it must read every source path from the local source manifest.

## R conventions

- Use `<-` for assignment and the base pipe `|>`.
- Prefer explicit package namespaces.
- Keep pure data-contract and summary functions under `app/R/` so tests can
  source them without starting Shiny.
- Derive normalized expression as
  `log1p(10000 * count / cell_library_size)`.
- Leave ratios missing when their denominator is zero.

Run `just test` after code changes. Run `just manifest` only from `app/`; do not
let exporter-only packages enter the deployment manifest.
