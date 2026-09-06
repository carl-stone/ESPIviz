# visual review

ESPIviz · reviewed and implemented 5 September 2026

The 14 primary findings, A1/A2, and the copy-edit recommendations have been addressed in the local application. The original review is retained below as the dated baseline; its descriptions of defects and verification limits describe the earlier state.

## Implementation record

| Finding | Implemented resolution | Verification |
| --- | --- | --- |
| 1 · Expression blend | Raw-count detection gates each gene's contribution. Detected values scale independently; constant detected values use 0.5. Legend, help, README, and exports explain the encoding. | Gm22307 contributes zero strength. Both named example pairs leave every double-negative cell neutral. Synthetic negative/constant-value checks pass. |
| 2 · Exported gene scope | Analysis downloads use current + second + saved-set genes across every page. Saved-set TXT remains explicitly separate. Counts appear beside exports. | Browser CSV/RDS downloads include Glul and Rlbp1. RDS reconstruction matches the public PFlog values exactly. Tests cover saved sets spanning pages. |
| 3 · Pathway figures | PNG/PDF have visible size glyphs, matching circles/diamonds, wrapped captions, and spacing based on wrapped term length. A quantitative size key and interpretation help precede the interactive plot. | Actual PNG and PDF downloads inspected visually; long labels, legends, and captions fit. |
| 4 · Unavailable adjusted P | Following the author's review, MA and volcano plots omit genes without adjusted P values. Gene cards state estimated direction separately from FDR status. Figures identify the highlighted gene or explain its omission. | Both plots show the same 17,428 eligible genes. All 7,173 missing adjusted P values remain in the results table and CSV; regression checks preserve tiny and zero P values as plotted results. |
| 5 · Linked result scope | Local scope notes distinguish selected-cell summaries from all-cell context. Saved-set size and additional current genes remain visible. Pathway exploration opens By cluster and identifies the loaded set and retained selection. | Browser pathway navigation and retained 17-cell selection checked. |
| 6 · Reproducible views | Copy view link captures supported settings and data identity. JSON saves support arbitrary selections and larger states; restoration validates the bundle version/checksum. | Fresh-link and uploaded-file restoration checked in the browser; arbitrary-cell state round trips and mismatched-bundle rejection tested. |
| 7 · Large sets | Automatic mode uses dot plots above six genes. Violin pages contain at most six genes; dot pages contain 25. Previous/Next controls appear only where they affect the result. | An 86-gene pathway set defaults to a dot plot; page navigation and irrelevant-tab hiding checked. |
| 8 · Readability | Larger helper text, clearer table spacing, larger composition annotations, responsive menus, and grouped gene statistics complete the visual pass. DE plots resize with the available width, reducing tick density and stacking the legend as needed. Wide pathway plots retain keyboard-accessible horizontal scrolling. | MA and volcano checked at nine viewport widths from 1,280 to 320 pixels, without plot or page-wide horizontal overflow. |
| 9 · Small and absent groups | Groups below ten cells use points and a median. Sample comparisons retain absent samples with zero cells and missing expression/detection percentages. | Cluster-8 sample counts and downloaded summaries checked; regression checks cover absent samples in tables and plot axes. |
| 10 · Plot interpretation | LOESS and its interval are described as cell-level summaries; No trend is available. Density normalization and included cells are explicit. Dot scales offer automatic/shared limits; zero detection has zero area, while missing groups use crosses. | Browser scatter/density controls and downloads checked; scale, detection-area, and summary tests pass. |
| 11 · Scientific tables | Readable labels, numerical formatting, method/direction/FDR filters, reset actions, result counts, and stated default sorting replace generic defaults. Pathway CSVs offer complete and filtered results. | Full pathway CSV has 11,089 rows; ORA filtering and DE missing-value filtering match their controls. All 185 marker rows match chosen-cluster versus remaining-cell detection percentages at the stored precision. |
| 12 · Reusable output | Explore results have dedicated PNG/PDF/data exports. Captions or metadata record genes, selection/model scope, and data version. Formats are named; sparse RDS reconstruction is documented. | Actual download buttons exercised across UMAP, summaries, cell distributions, pair plots, composition, markers, DE, and pathways. Twelve PDF variants contain a data version and have text within page bounds. |
| 13 · Accessibility | Tables receive descriptive accessible names. Selectable rows support Enter/Space; individual cells can be selected through searchable controls. Expanded cards restore focus and announce their state. | Keyboard gene/pathway/cell selection, menu Escape, expanded-card focus trapping/return, and four-page axe checks pass with zero automated violations. |
| 14 · Terminology/provenance | Comparison labels use p27CKO + E-Stim vs p27CKO. Expression axes use log normalized expression. The About Methods section was removed at the author's request; normalization remains documented in the README and export metadata. Data version, citation/deposit status, code license, and provisional data reuse status are visible. The author line is Megan L. Stone, Carl Stone, Joel Jovanovic, Edward M. Levine, in the requested order. | Names verified against the [bioRxiv preprint author record](https://pubmed.ncbi.nlm.nih.gov/41727120/). Copy and public-bundle metadata checked. The requested DE introduction and source-condition-mean disclaimer were removed. Complete result exports retain their original numerical fields. |
| A1/A2 · Working content and control states | Table/export polish follows the scientific contracts above. Empty sets have an actionable message; removal, clearing, empty downloads, and pagination reflect the current state. | Empty-set and mixed valid/invalid pasted-gene flows, menu layouts, and large-set navigation checked. |

## Implementation verification

The follow-up NA check found no zero P values in the public bundle. Of its 7,173 missing adjusted P values, 7,155 have finite raw P values and mean counts below the lowest mean count with an adjusted P value; the other 18 also lack raw P values. This pattern is consistent with DESeq2 filtering, rather than numerical underflow. The source volcano implementation already omits missing adjusted P values. The app now applies that omission to both DE plots, preserving the table and CSV values. [DESeq2 results documentation][deseq2]

`just test` passes. The optional committed-release-asset test is skipped because the generated bundle is kept outside Git; the actual public v1.2.0 bundle was separately loaded through the checksummed app loader. `just manifest` was run from `app/`, and `git diff --check` passes. The immutable data bundle and exporter normalization contract are preserved.

Browser coverage includes all four main pages and all seven Explore tabs, actual figure/data downloads, filter changes, a second gene, single-cell and cluster selection, large-set pagination, copied links, saved-file upload, keyboard interactions, and narrow layouts. CSV and RDS contents were read back. Pathway PNG/PDF and representative other exported figures were inspected visually; all twelve checked PDFs have complete text bounds and a data-version label.

This verifies the local application in Chrome on macOS. Automated accessibility checks and keyboard checks are not a formal screen-reader audit; viewport tests are not physical touch-device tests. The work does not revalidate the upstream study analysis or establish a new public deployment.

## Original review baseline

The original highest priorities were scientific interpretation and reproducibility: an undetected gene could contribute a strong expression-blend color; visible and downloaded gene scopes differed; the pathway export had a broken size legend and clipped caption; and the MA plot grouped unavailable adjusted P values with nonsignificant results.

## Assessment by area

| Area | Assessment | Main opportunity |
| --- | --- | --- |
| Visual style | Remaining concerns are concentrated in dense labels, tables, and exported figures. | Finish readability at the actual display and output sizes. |
| UI/UX | Core controls work, but users must track several kinds of state across a long page. | Put scope, navigation, and feedback beside the result they affect. |
| Scientific plotting | Several interpretation and export defects remain. | Correct misleading encodings, clarify small-group summaries and smoothing, and verify exported figures visually. |
| Data exploration | Broad capabilities, particularly replicate summaries and complete result tables. | Make linked exploration predictable, shareable, and easier to turn into reusable figures or tables. |
| Prose and instructions | Useful caveats coexist with jargon, repetition, inconsistent labels, and generic empty states. | Use short, contextual explanations and consistent scientific terminology. |

## Scope and evidence

The review combined live inspection of the [public application][app], source review at commit [22b146b][revision], and numerical checks using the public v1.2.0 bundle. The bundle checksum was validated through the app loader. Its contents include 3,456 cells, eight clusters, six samples, 38,394 genes, 24,601 primary DE rows, and 11,089 enrichment results. The deployed app's exact commit was not independently established.

This cleanup compares the current local working tree with the review and uses the desktop and narrow-viewport checks performed during the redesign. The deployment has not been rechecked. Pinned source links document the original evidence; links labeled “Current source” refer to the files now under review.

Live coverage included all four primary pages and all seven Explore result tabs; gene and second-gene search; cluster coloring; selection of cluster 8; selection comparisons; violin, dot, scatter, and density views; sample composition; marker summaries; DE search and MA/volcano switching; pathway details and gene-set navigation; a pasted list containing an invalid symbol; an empty gene set; and opening a copied exploration URL in a fresh session. Layout was inspected at 1280 × 720, 390 × 844, and 320 × 844 CSS pixels.

Static MA and pathway figures were generated from the same functions used by downloads. The pathway PNG was also rendered at the configured 9 × 11 inches and 320 dpi. This was a check of export rendering, not an end-to-end test of every download button or PDF.

Evidence labels distinguish **Live**, **Source + bundle**, and **Rendered export** observations. Accessibility findings are a DOM/accessibility-tree and visual spot check, not a formal conformance or screen-reader audit. This cleanup edits the review document only.

## Priority findings

**P1** means resolve before relying on the affected view for scientific interpretation or sharing. **P2** means an important usability or presentation improvement. Suggested changes are recommendations, not implemented behavior.

### 1. P1 — Undetected genes can create an expression-blend gradient

**Evidence: Source + bundle.** The expression blend independently rescales each gene's centered PFlog values across all cells. It does not use the raw-count detection flags to restrict the strength calculation. An undetected gene therefore varies with the per-cell centering offset and can receive a large relative strength. [Implementation][blend]

A concrete reproduction is to pair **Glul** with **Gm22307**. Gm22307 has zero raw counts throughout this public bundle, but its computed blend strength spans **0 to 1**. For **Mki67 + Ascl1**, 2,398 cells detect neither gene; 2,397 nevertheless receive a non-neutral blend color. The display problem is that a rescaled centering offset can be interpreted as gene strength; this finding does not establish an error in PFlog normalization.

The legend's “High [gene]” wording invites an expression or coexpression interpretation. The README also explicitly describes a detection-gated blend in which double-negative cells remain neutral, which differs from the implementation. [README][readme]

**Recommendation:** resolve the intended encoding explicitly. A detection-gated expression blend should assign zero strength to undetected genes and scale expression among detected cells, with defined behavior for constant values. At minimum, prevent an entirely undetected gene from producing a “high” gradient and explain what low/high mean beside the legend. Keep the distinct raw-count detection view. Align the README, About text, hover text, and exported caption with the final behavior.

**Acceptance check:** an entirely undetected gene contributes no apparent high-expression signal; double-negative behavior matches the documentation; independent scaling is disclosed and is not presented as an expression ratio.

### 2. P1 — Downloaded summaries do not necessarily contain the genes visible in the summaries

**Evidence: Live; Source + bundle.** With no active gene set, selecting **Glul** and second plot gene **Rlbp1** produces two-gene plots and comparison-table rows. The Gene set, Selection summary, and Cell expression download handlers use `analysis_genes()`, which returns only Glul in this state. With an active set, visible summaries include the current and second genes as well as the set, but the same downloads use only the set. [Visible gene scope][explore-scope] · [Download handlers][explore-downloads] · [Download scope helper][state-genes]

A user downloading the visible “Selection summary” can reasonably expect both displayed rows to be included.

**Recommendation:** establish one visible export contract. Either export the full displayed analysis gene scope, or make the distinction explicit with actions such as “Download plotted genes” and “Download gene set.” Show the number of genes and cells beside each data download. Preserve all pages of the chosen scope.

**Acceptance check:** the interface states exactly which genes each export contains, and a two-gene comparison can be downloaded without manually reconstructing the missing row.

### 3. P1 — The pathway PNG has an unusable size legend and a clipped caption

**Evidence: Rendered export.** At the download's configured dimensions and resolution, the size legend shows the labels **20, 40, 60** without visible size symbols. The final caption extends beyond the right edge. The static export also uses squares for E-Stim while the interactive plot uses diamonds. [Static pathway export][pathway-static]

**Recommendation:** give the size guide an explicit visible glyph and fill/stroke combination; wrap and position the caption within the full figure width; and use the same direction shapes on screen and in exported figures. Make a readable figure at the final physical size the acceptance criterion.

The interactive chart's explanation of method-specific scales and point size sits below the 1,200-pixel figure. Put essential interpretation guidance near the start of the chart. [Current source][current-pathways]

The interactive chart has a direction legend but no quantitative size key. Add a size key labeled “−log10 adjusted P” or simplify the encoding and expose adjusted P in an aligned value column. “Adjusted-P strength” is too vague as the only explanation of size.

**Acceptance check:** PNG and PDF both have visible size keys, complete captions, consistent direction symbols, and readable long terms. The current review reproduced the PNG defects; PDF rendering still needs a dedicated check.

### 4. P1 — “Not significant” conflates unavailable and evaluated adjusted P values in the MA view

**Evidence: Live; Source + bundle.** **7,173 of 24,601 genes (29.2%)** have no adjusted P value. The MA view places them in “Not significant.” The explanatory paragraph acknowledges this choice, and the volcano view gives them a separate category, but the MA legend and its export do not preserve that distinction. [Classification][de-classification] · [Plot and export][de-plots]

**Recommendation:** use “Adjusted P unavailable” in both views. Label evaluated results that exceed the threshold as “FDR > 0.05,” and show the threshold with the MA legend as well as on the volcano plot. Preserve missing values in data exports. In the current-gene card, distinguish the sign of the estimated effect from its FDR status.

**Acceptance check:** absence of an adjusted P value cannot be mistaken for a tested nonsignificant result, including when a figure is viewed outside the app. The DESeq2 documentation distinguishes unavailable adjusted P values arising from filtering and other result conditions. [DESeq2 documentation][deseq2]

### 5. P2 — Result scope remains ambiguous during linked exploration

**Evidence: Current source; earlier live check.** “Explore gene set” changes the saved set and current gene but preserves the cell selection, optional second gene, and Explore result tab. In the earlier reproduction, loading 86 leading-edge genes after selecting cluster 8 retained the 17-cell selection and Cell-level tab. The header reported “Genes 1–25 of 87,” while Cell-level displayed only the current pair. This behavior remains in the pathway action and result-scope logic. [Current source: pathway actions][current-pathways] · [Current source: Explore scopes][current-explore]

By cluster, By sample, and Pooled condition summaries use selected cells, whereas Cell-level cluster distributions and the composition overview provide all-cell context. Their local headings do not consistently identify that distinction. The saved gene-set count is inside the expandable Gene set menu, while the result header counts the combined analysis genes.

**Recommendation:** label the cell and gene scope of each result, especially views that show all-cell context despite an active selection. After a pathway action, identify the loaded set and retained selection, and open a result tab that displays the set. Distinguish the saved set size from additional current/second genes.

### 6. P2 — A copied URL restores only a small part of the exploration

**Evidence: Live; Source.** Opening the current Explore URL in a new session preserved Mre11a but lost the gene set, second gene, 17-cell selection, color mode, and result tab. The source writes view, current gene, and sometimes pathway to the URL; it does not serialize most exploration settings. It can read a `genes` parameter but does not preserve it in its generated URL. [URL handling][url-state]

**Recommendation:** add an explicit “Copy view link” with a clear description of what it includes. Encode safe, compact settings such as genes, cluster selection, plot type, and grouping; use a local state download for more complex cell selections if needed. Include data version in saved state. Avoid implying that copying the address reproduces an arbitrary selection.

### 7. P2 — Large gene sets still produce long figures and irrelevant pagination

**Evidence: Current source.** A 25-gene violin page now requests **2,990 pixels** of height: 13 rows at 230 pixels per row. Larger sets still default to violins. The Page control is shown whenever the analysis spans multiple pages, including on Cell-level and Cluster markers, whose plots are not driven by that page. [Current source: figure sizing][current-plots] · [Current source: pagination][current-explore]

**Recommendation:** default larger gene sets to a compact dot plot, with violins available on demand. Use “Gene page 1 of 4” with Previous/Next controls, and hide pagination on views where it has no effect.

### 8. P2 — Some dense labels and helper text still need a readability pass

**Evidence: Current source.** Composition tile counts and percentages retain a fixed ggplot text size of 3. Supporting copy, figure notes, and table headings remain 12-pixel CSS text, including long quantitative labels. These dense cases still need checking at their final display size. [Current source: composition plot][current-plots] · [Current source: typography][current-styles]

**Recommendation:** check composition labels and long quantitative table headings at ordinary browser zoom and narrow widths. Increase the remaining small annotations where needed and shorten repeated metric labels after defining their units.

### 9. P2 — Small groups and absent samples need clearer treatment

**Evidence: Current source; earlier bundle check.** Within selected cluster 8, the five represented samples contain **1, 5, 1, 1, and 9 cells**. Singletons are shown as points, but groups with two or more cells still use a smoothed violin regardless of how small the group is. Samples with no selected cells are omitted from these summaries. [Current source: violin grouping and rendering][current-plots]

**Recommendation:** prefer points and a simple summary for very small groups. Keep absent samples visible as “0 selected cells” where the goal is replicate comparison.

### 10. P2 — Several plot encodings need more precise interpretation support

**Evidence: Live; Source.** The gene-pair scatter includes a LOESS line and a 95% confidence ribbon, but the ribbon's meaning is not explained in the visible instructions or legend. Its fit is based on cells, not a model of biological-replicate uncertainty. The density view says “Cell density” without explaining that it is a smoothed probability density. [Pair plotting][pair-plots]

Expression color limits are recalculated from each displayed gene/summary, and the blend rescales genes independently. Dot and marker plots give 0% detection a nonzero glyph size. These choices can be defensible, but comparisons between pages or views require explicit guidance. [Expression dot plot][dot-plot]

**Recommendation:** label the smooth and its interval, identify it as descriptive, and provide a “No trend” option. Explain density normalization and the cells included in it. Add a clear automatic/shared color-scale choice for multi-gene comparisons. Use an unambiguous zero-detection symbol or zero area, distinct from missing data. Preserve the existing useful statement about excluding double-negative cells from pair plots, with the longer centering explanation available in help.

### 11. P2 — Results tables expose implementation terminology and weak default organization

**Evidence: Live; Source + bundle.** Pathway columns include `pathway_id`, `source`, `p_value`, `p_adjust`, `score`, and `gene_count`. The meaning of `score` and `gene_count` varies by method. Marker columns include `avg_log2FC`, `pct.1`, `pct.2`, and `p_val_adj`; their units and comparison groups are not explained locally. The pathway table still includes a description column that duplicates its label for **all 11,089 rows**. [Pathway tables][pathway-tables] · [Marker table][marker-table]

The default DE table begins in source order, with Xkr4, Gm53491, Rp1, Sox17, and Gm37323, rather than a stated scientific ranking. Search works, but generic “All” column filters and horizontal scrolling make common tasks less discoverable. [DE table][de-table]

**Recommendation:** use readable display labels and format numerical columns deliberately; define percentage versus fraction; identify the marker reference group; label NES and fold enrichment distinctly; and remove the duplicated description column from the pathway table. Provide explicit method, direction, and adjusted-P filters for pathways and DE, with result counts and a reset action. State or offer the default sort. Add complete and filtered enrichment-table downloads, parallel to the existing complete DE download.

### 12. P2 — Reusable output is uneven across plot types

**Evidence: Source; Rendered export.** UMAP, DE, and pathway plots have dedicated image-download actions. The violin, dot, marker, and composition views have no equivalent app buttons. Pair plots have a Plotly image action but no matching dedicated PDF/data action. The MA export highlights Glul with an outline but does not name the highlighted gene in a title or caption. Exported figures generally lack bundle version and selection context. [Explore downloads][explore-downloads] · [DE export][de-plots]

**Recommendation:** provide consistent PNG, PDF where appropriate, and underlying-data downloads across plot types. Put gene names, selected-cell scope, data version, and relevant scale/model information into reusable captions or an accompanying metadata file. Label TXT, CSV, and RDS formats on the buttons. Explain that Cell expression is an RDS with a shifted sparse matrix and a separate centering vector, and provide a short reconstruction example alongside it. [Expression export contract][expression-export]

### 13. P2 — Table filters and keyboard workflows need an accessibility pass

**Evidence: Earlier live DOM/accessibility check; current source.** Several table filters are exposed simply as “All,” without identifying the column. Full keyboard and assistive-technology coverage remains unverified. [Current source: enrichment tables][current-pathways] · [Current source: Explore tables][current-explore]

**Recommendation:** give table filters column-specific accessible names. Verify keyboard-only selection, table-row activation, focus order, and return from expanded cards. Keep table alternatives to color-dependent or pointer-dependent plots. Do not treat the viewport and interaction spot checks as a completed accessibility audit.

### 14. P2 — Condition terminology, methods, and provenance need clearer explanation

**Evidence: Current source; earlier live check.** The p27CKO and E-Stim labels still need concise definitions and an explicit mapping to “Control” and “E-Stim” in the inferential views. About contains unexplained “MG-selected,” “Mouse × Condition,” PFlog parameters, and design syntax. GSEA, ORA, NES, and leading-edge genes are not expanded near their first use. [Current source: Explore context][current-explore] · [Current source: About][current-about] · [Current source: pathway details][current-pathways]

The Data card provides the bundle and application-code links, but does not visibly state the bundle version, manuscript/deposit status, citation guidance, or the provisional reuse notice. The fallback sentence about future records is unreachable while the default code link exists. The repository's data-license notice already states that final deposit terms are pending; surface that status without inventing an accession, DOI, or license. [Data links][data-links] · [Data notice][data-license]

**Recommendation:** define the condition names and map the two labeling systems explicitly; use a small glossary or inline help for method names; put implementation details behind a methods disclosure; and add a concise provenance block with the actual data version and available public records.

## Copy-edit recommendations

These are proposed replacements. They should be adjusted to the final behavior where a P1 finding changes an encoding or export contract.

| Location / current text | Suggested text or treatment |
| --- | --- |
| “Gene” / repeated “current gene” | Use “Current gene” consistently for the gene shared across pages. |
| “Second plot gene (optional)” | “Second gene (optional)” followed by a precise statement of which plots, tables, and downloads it affects. |
| “Log normalized expression” repeated in axes, headers, and sentences | Define “PFlog expression” once in nearby help, then use it consistently. Explain that values are centered and may be negative, and that detection means raw count > 0. |
| “No region selected” | “All cells — no selection applied.” |
| “Click a cell or use box or lasso on the UMAP.” | “Click a cell, drag a lasso, or choose a cluster. Use the plot toolbar to switch to box selection.” |
| “Use featured set” | “Replace gene set with featured set,” or a shorter label with an adjacent replacement notice. |
| “Remove chosen” | “Remove selected genes.” |
| “Not found: not_a_gene” | “Added 2 genes. Unrecognized symbol: not_a_gene.” Report newly added, already present, and unrecognized counts when relevant. |
| Empty Gene set tab: “No data available in table” | “No gene set yet. Choose a featured set, paste gene symbols, or add the current gene.” Hide the empty pagination controls. |
| “Page” / “Genes 1–25 of 87” | “Gene page 1 of 4” and “Showing 25 of 87 analysis genes.” Explain any additional current/second genes outside the saved set. |
| “Observed cluster composition” | “Observed cluster composition — all cells.” Keep the existing descriptive-analysis caveat. |
| “Other cells” and “Remaining cells” used interchangeably | Use “Remaining cells” consistently. |
| “Mean difference” / “Detection difference (pp)” | “Selected − remaining mean” and “Detection difference,” with “percentage points” explained in help or a units row. |
| “Mean Count Control” / “Mean Count E-Stim” | Specify whether these are raw or normalized counts and the averaging unit, using the verified export definition. Avoid making the label imply a unit that has not been established. |
| “Pathways” with “Top directional GSEA and ORA results…” | Keep the short tab label; use “GO Biological Process enrichment” as the page description. Expand gene set enrichment analysis (GSEA), over-representation analysis (ORA), and normalized enrichment score (NES) on first use. |
| “Explore gene set” / “Pathway genes” | For GSEA: “Explore 86 leading-edge genes” / “Leading-edge genes (TXT).” For ORA, use “overlapping genes.” Generate the count dynamically. |
| “Cell expression” / “Selection summary” | Include format and scope: for example, “Cell expression (RDS)” and “Selection summary (CSV),” with gene and cell counts nearby. |
| “for 1 genes” in plot alternate text | Use singular/plural-aware wording: “for 1 gene” / “for 2 genes.” |
| About: “frozen final MG-selected analysis object” | “The map shows cells retained for the final Müller glia analysis, using the study's final UMAP coordinates and cluster assignments.” Confirm that “Müller glia” accurately expands MG before applying. |
| Maintenance: “Please try again later.” | Keep the clear unavailable message and provide “Reload” plus a public source/data link when useful. Do not expose internal loader errors. |

## Recommended revision order

1. **Correct interpretation and export reliability:** findings 1–4. Recheck absent-gene blends, unavailable adjusted P values, visible/exported gene scopes, and the pathway figure at its actual output size.
2. **Make exploration state understandable:** findings 5–7. Clarify per-result scope, pathway navigation, saved-view links, and large-set pagination.
3. **Finish readability and reuse:** findings 8–14 and the copy table. Address remaining small annotations, small-group behavior, interpretation guidance, tables, downloads, terminology, and accessibility.

Preserve the immutable public-bundle architecture and the clear distinction between descriptive exploration and the study's fixed inferential results.

## Verification status and limits

`just test` completed successfully. One test was skipped because the optional local release asset was not present in the repository. The public bundle used in this review was independently loaded and checksum-validated through the app's loader. The redesign's final `just test` run also passed with the same optional check skipped. This document-only cleanup did not rerun the code tests; passing tests do not close the remaining review findings.

The review did not validate the underlying study design or upstream statistical analysis, test every gene/pathway, exercise all lasso and single-cell edge cases, test file-upload round trips, audit every download/PDF, or run a full accessibility/browser-compatibility suite. Small-screen observations are viewport tests, not physical touch-device tests. Maintenance behavior was reviewed from source rather than induced in the public application.

[app]: https://019f6264-f0b5-e432-795c-c2f7e9fd5c95.share.connect.posit.cloud/
[revision]: https://github.com/carl-stone/ESPIviz/tree/22b146b5be40c9aeb77fcc4518f76082ed2f4b23
[readme]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/README.md
[blend]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/04-cell-level.R#L58
[explore-scope]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/10-explore.R#L617
[explore-downloads]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/10-explore.R#L1334
[state-genes]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/03-state-plots.R#L33
[pathway-static]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/30-pathways.R#L345
[pathway-tables]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/30-pathways.R#L682
[de-classification]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/20-differential-expression.R#L15
[de-plots]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/20-differential-expression.R#L391
[de-table]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/20-differential-expression.R#L697
[url-state]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/90-app-shell.R#L142
[pair-plots]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/04-cell-level.R#L399
[dot-plot]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/03-state-plots.R#L928
[marker-table]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/10-explore.R#L1322
[expression-export]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/02-expression.R#L453
[data-links]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/app/R/40-about.R#L141
[data-license]: https://github.com/carl-stone/ESPIviz/blob/22b146b5be40c9aeb77fcc4518f76082ed2f4b23/DATA_LICENSE.md
[deseq2]: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## AI slop design review

Remaining issues after the redesign · 5 September 2026

The remaining concerns are unfinished table/export presentation and controls that do not reflect the current task state. “AI slop” here describes generic, insufficiently task-specific design decisions, not evidence about who or what produced the application.

### A1. Working tables and exports still expose unedited defaults

**Evidence:** Enrichment and marker tables retain technical column names and weak numeric formatting. Pathway figures still have the export legend and caption defects documented above, and downloads often omit format and scope from their labels. [Current source: tables and pathway exports][current-pathways] · [Current source: Explore tables and downloads][current-explore]

**Why it reads as slop:** Generic table and export components leave readers to translate implementation details into scientific meaning. The working content needs the same attention as the surrounding interface.

**Change:** finish the scientific labels, units, numeric formatting, export scope, and figure captions described in the remaining primary findings. Verify exported figures at their actual output size.

### A2. Gene-set controls still need to reflect empty and irrelevant states

**Evidence:** The expandable Gene set editor always includes Remove chosen and Clear set, even with no set to edit. The empty Gene set table still relies on its generic empty-table message. For a set spanning multiple pages, Page remains available on result tabs where it does not affect the plot. [Current source][current-explore]

**Why it reads as slop:** Controls are still exposed according to component availability rather than whether the current data and task make them useful.

**Change:** hide or disable removal actions until they apply, provide an actionable empty-set message, and show pagination only on views driven by the gene page. Keep saved-set scope visible when it is active.

### Criteria for the remaining visual work

| Check | Desired result |
| --- | --- |
| Open an empty Gene set | The message explains how to add genes; removal and empty pagination controls do not compete. |
| Browse a large gene set | A readable summary fits the task, and pagination affects the visible result. |
| Read quantitative tables | Labels, units, numeric formatting, and filters are understandable at ordinary zoom. |
| Export a result | The control states format and scope; the file has legible labels, a usable legend, and a complete caption. |

Preserve the scientific specificity: mouse retinal cells, the stimulation contrast, replicate summaries, named genes, GO terms, fixed analysis outputs, and useful interpretation caveats.

[current-explore]: app/R/10-explore.R
[current-plots]: app/R/03-state-plots.R
[current-styles]: app/www/styles.css
[current-pathways]: app/R/30-pathways.R
[current-about]: app/R/40-about.R
