# ESPIviz

ESPIviz is the public interactive companion to the ESPI single-cell manuscript.
It provides final-clustering UMAP exploration, gene and gene-set expression
summaries, the complete primary condition-model differential-expression result,
and a short set of manuscript-aligned pathways.

**ESPIviz — analysis and application by Carl Stone.**

The interface reports values directly and does not calculate differential
expression from ad hoc cell selections. A selection can contain one cell, many
cells, or every cell. Expression uses the same normalization as the frozen
analysis object:

```text
log1p(10000 × raw count / cell library size)
```

## Repository contents

- `app/`: modular Shiny application and its Connect Cloud manifest
- `scripts/`: standalone, manifest-driven bundle exporter
- `config/`: featured genes, featured pathways, and source-manifest example
- `tests/`: synthetic contract and behavior tests

No manuscript object, original barcode, private path, or exhaustive enrichment
output is committed here. The application downloads a versioned public bundle
whose URL and SHA-256 checksum live in `app/data-manifest.json`.

## Local development

Install R 4.6 and the packages declared in `DESCRIPTION`, then use:

```sh
just data-dry-run
just data-build
just app-run
just test
just manifest
```

Copy `config/source-manifest.example.json` to the gitignored
`config/source-manifest.local.json` and fill in explicit local source paths
before building data. `just data-dry-run` validates those inputs without writing
a bundle. To run against a local bundle:

```sh
just app-run /absolute/path/to/espiviz-data-v1.0.0.rds
```

`just manifest` scans `app/` only, keeping Seurat and exporter dependencies out
of production.

## Deployment

Publish the immutable data release before the first application deployment.
Generate `app/manifest.json` with `just manifest`, then configure Posit Connect
Cloud to deploy `app/` from the public repository's `main` branch with automatic
deployment on pushes. Connect Cloud supports public GitHub sources on its
[free plan](https://docs.posit.co/connect-cloud/user/account/plans.html) and R
deployments described in its [R requirements](https://docs.posit.co/connect-cloud/user/platform/r.html).

Application-code pushes may reuse the pinned data release. A data change always
requires a new immutable release asset and a corresponding
`app/data-manifest.json` update.

## Data and licensing

Application code is MIT licensed. The bundle follows the reuse terms of the
eventual public data repository; see `DATA_LICENSE.md`. Manuscript DOI and final
data links will be added after those public records exist.

Do not replace a published bundle in place. Data changes require a new immutable
release, checksum, schema-compatible manifest entry, and application deployment.
