provenance_value <- function(provenance, names) {
  for (name in names) {
    value <- compact_character(provenance[[name]])
    if (length(value) > 0L) return(value[[1L]])
  }
  NA_character_
}

https_url <- function(value) {
  value <- compact_character(value)
  if (length(value) == 0L ||
    !grepl("^https://[^[:space:]]+$", value[[1L]], perl = TRUE)) {
    return(NA_character_)
  }
  value[[1L]]
}

published_bundle_url <- function(bundle) {
  manifest <- attr(bundle, "data_manifest") %||% list()
  hash <- manifest_sha256(manifest)
  url <- https_url(public_bundle_url(manifest))
  if (is.na(url) || is.na(hash) || !grepl("^[0-9a-f]{64}$", hash) ||
    identical(hash, paste(rep("0", 64L), collapse = ""))) {
    return(NA_character_)
  }
  url
}

about_ui <- function(id) {
  ns <- shiny::NS(id)
  htmltools::div(
    class = "about-shell",
    htmltools::tags$section(
      class = "about-hero",
      htmltools::p("About the study", class = "eyebrow"),
      htmltools::h1(
        "Electrical stimulation combined with p27Kip1 inactivation drives proliferative neurogenic reprogramming of Müller glia in the adult mouse retina"
      ),
      htmltools::p(
        "ESPIviz is the interactive companion to the study's final single-cell RNA-seq analysis.",
        class = "lede"
      ),
      htmltools::p(
        "ESPIviz — analysis and application by Carl Stone.",
        class = "credit-line"
      )
    ),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        bslib::card_header("Study"),
        htmltools::p(
          "The experiment compares p27CKO retinal cells with p27CKO cells collected after electrical stimulation. The final explorer contains 3,456 cells across eight clusters."
        ),
        htmltools::p(
          "Use Explore to inspect genes and cell selections, Differential expression to search the primary condition model, and Pathways to open the featured manuscript-aligned gene sets."
        )
      ),
      bslib::card(
        bslib::card_header("Data"),
        shiny::uiOutput(ns("data_links")),
        shiny::downloadButton(ns("download_bundle"), "Download processed app bundle", class = "btn-primary")
      ),
      bslib::card(
        class = "span-12 methods-card",
        bslib::card_header("Methods"),
        htmltools::div(
          class = "methods-grid",
          htmltools::div(
            htmltools::h2("Expression"),
            htmltools::p(
              "Normalized expression is calculated from raw counts for each cell as log1p(10,000 × count / cell library size)."
            )
          ),
          htmltools::div(
            htmltools::h2("Cell map"),
            htmltools::p(
              "The UMAP coordinates and numeric cluster assignments come from the frozen final MG-selected analysis object."
            )
          ),
          htmltools::div(
            htmltools::h2("Condition model"),
            htmltools::p(
              "The differential-expression view contains all 24,601 genes from the primary six-sample Mouse × Condition pseudobulk model with design ~ condition."
            )
          ),
          htmltools::div(
            htmltools::h2("Cell selections"),
            htmltools::p(
              "Selection summaries report expression and detection for selected cells and all remaining cells. Ratios are left blank when their denominator is zero."
            )
          )
        )
      )
    )
  )
}

data_links_ui <- function(bundle) {
  provenance <- bundle$provenance %||% list()
  links <- list()
  bundle_url <- published_bundle_url(bundle)
  data_url <- https_url(provenance_value(provenance, c("data_url", "data_repository_url")))
  code_url <- https_url(provenance_value(provenance, c("code_url", "application_code_url")))
  if (is.na(code_url)) code_url <- "https://github.com/carl-stone/ESPIviz"
  doi <- provenance_value(provenance, c("manuscript_doi", "doi"))
  if (!is.na(bundle_url)) {
    links <- c(links, list(htmltools::tags$li(htmltools::a(
      "Processed app bundle",
      href = bundle_url,
      target = "_blank",
      rel = "noopener noreferrer"
    ))))
  }
  if (!is.na(data_url)) {
    links <- c(links, list(htmltools::tags$li(htmltools::a(
      "Public data repository",
      href = data_url,
      target = "_blank",
      rel = "noopener noreferrer"
    ))))
  }
  if (!is.na(code_url)) {
    links <- c(links, list(htmltools::tags$li(htmltools::a(
      "Application source code",
      href = code_url,
      target = "_blank",
      rel = "noopener noreferrer"
    ))))
  }
  if (!is.na(doi)) {
    doi_url <- if (grepl("^10[.][0-9]{4,9}/[^[:space:]]+$", doi)) {
      paste0("https://doi.org/", doi)
    } else {
      https_url(doi)
    }
    if (!is.na(doi_url)) {
      links <- c(links, list(htmltools::tags$li(htmltools::a(
        "Manuscript",
        href = doi_url,
        target = "_blank",
        rel = "noopener noreferrer"
      ))))
    }
  }
  if (length(links) == 0L) {
    htmltools::p(
      "Public manuscript and repository records will appear here when available.",
      class = "supporting-copy"
    )
  } else {
    htmltools::tags$ul(class = "data-link-list", links)
  }
}

about_server <- function(id, bundle) {
  shiny::moduleServer(id, function(input, output, session) {
    output$data_links <- shiny::renderUI({
      data_links_ui(bundle)
    })

    output$download_bundle <- shiny::downloadHandler(
      filename = function() {
        paste0("espiviz-data-v", as.character(bundle$data_version), ".rds")
      },
      content = function(file) {
        source <- attr(bundle, "bundle_path")
        if (!is.null(source) && file.exists(source)) {
          if (!file.copy(source, file, overwrite = TRUE)) {
            stop("The bundle could not be copied.")
          }
        } else {
          clean_bundle <- bundle
          attributes(clean_bundle) <- attributes(clean_bundle)[
            setdiff(names(attributes(clean_bundle)), c("bundle_path", "data_manifest"))
          ]
          saveRDS(clean_bundle, file, compress = "xz")
        }
      }
    )
  })
}
