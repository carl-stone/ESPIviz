# ESPIviz

ESPIviz is the public interactive companion to the ESPI single-cell manuscript.
It provides final-clustering UMAP exploration, gene and gene-set expression
summaries across clusters, conditions, and biological samples, descriptive
sample-level cluster composition, marker overviews, two-gene UMAP blends,
cluster violins, per-cell gene-pair scatter plots, the complete primary
condition-model differential-expression result with MA and volcano views, and a
curated short set of manuscript-aligned pathways shown on method-specific score
scales.

**[Open ESPIviz](https://019f6264-f0b5-e432-795c-c2f7e9fd5c95.share.connect.posit.cloud/)**

**ESPIviz — analysis and application by Carl Stone.**

The interface reports values directly and does not calculate differential
expression from ad hoc cell selections. A selection can contain one cell, many
cells, or every cell. Expression is derived from raw counts with
[`scclrR` PFlog](https://github.com/cleartools/scclrR), using its automatic
overdispersion target. For cell *i* and gene *j*:

```text
PFlog[i,j] = log1p(4 × alpha × count[i,j]) - cell_center[i]
```

The center is calculated across the complete gene universe, so a zero count can
have a negative PFlog value. Detection percentages always use raw counts. The
standalone exporter runs `scclrR` and stores its shifted sparse matrix and
per-cell center vector in the bundle. The Shiny runtime only reconstructs dense
values by subtracting those exported centers; it does not implement PFlog or
install the Rust-backed package in the deployment image.

For a two-gene UMAP blend, raw counts determine whether each gene is detected.
Color intensity then increases with that gene's PFlog strength among detected
cells; double-negative cells remain neutral and double-positive cells mix the
two feature colors. Violin and gene-pair scatter plots retain the centered PFlog
scale and report raw-count detection separately.

The pinned v1.0.0 data asset remains byte-for-byte immutable. Bundle schema
v1.1.0 adds the exact normalization values returned by the pinned `scclrR`
revision.

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

The standalone exporter and PFlog equivalence test additionally require Rust
and the pinned `scclrR` source revision:

```sh
Rscript -e 'pak::pak("cleartools/scclrR@6d6378dd41502a8606da14adb01a01032cb75224")'
```

Copy `config/source-manifest.example.json` to the gitignored
`config/source-manifest.local.json` and fill in explicit local source paths
before building data. `just data-dry-run` validates those inputs without writing
a bundle. To run against a local bundle:

```sh
just app-run /absolute/path/to/espiviz-data-v1.1.0.rds
```

`just manifest` scans `app/` only, keeping Seurat, `scclrR`, and other exporter
dependencies out of production.

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
