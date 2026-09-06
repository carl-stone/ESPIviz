provenance_value <- function(provenance, names) {
  for (name in names) {
    value <- compact_character(provenance[[name]])
    if (length(value) > 0L) return(value[[1L]])
  }
  NA_character_
}

https_url <- function(value) {
  value <- compact_character(value)
  if (
    length(value) == 0L ||
      !grepl("^https://[^[:space:]]+$", value[[1L]], perl = TRUE)
  ) {
    return(NA_character_)
  }
  value[[1L]]
}

published_bundle_url <- function(bundle) {
  manifest <- attr(bundle, "data_manifest") %||% list()
  hash <- manifest_sha256(manifest)
  url <- https_url(public_bundle_url(manifest))
  if (
    is.na(url) ||
      is.na(hash) ||
      !grepl("^[0-9a-f]{64}$", hash) ||
      identical(hash, paste(rep("0", 64L), collapse = ""))
  ) {
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
      htmltools::h1("About the study"),
      htmltools::h2(
        class = "manuscript-title",
        "Electrical stimulation combined with p27Kip1 inactivation drives proliferative neurogenic reprogramming of Müller glia in the adult mouse retina"
      ),
      htmltools::p(
        "Megan L. Stone · Carl Stone · Joel Jovanovic · Edward M. Levine",
        class = "manuscript-authors"
      ),
      htmltools::p(
        "ESPIviz is an interactive companion to the study's single-cell RNA-seq analysis.",
        class = "lede"
      ),
      htmltools::p(
        "ESPIviz application by Carl Stone.",
        class = "credit-line"
      )
    ),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        bslib::card_header("Study"),
        htmltools::p(
          "Comparisons are p27CKO + E-Stim vs p27CKO: p27 inactivation with electrical stimulation versus p27 inactivation alone.",
          class = "condition-key"
        ),
        htmltools::p(
          "The study examines Müller glia in the adult mouse retina after p27 inactivation (p27CKO), with and without electrical stimulation (E-Stim). The final explorer contains 3,238 cells across five clusters and six biological samples. Clustering uses PFlog normalization, no cell-cycle filtering, 20 principal components, and resolution 0.3."
        ),
        htmltools::p(
          paste(
            "Use Explore to inspect genes and cell selections, Differential",
            "expression to search the primary condition model, and Pathways",
            "to search the complete Gene Ontology Biological Process",
            "enrichment results."
          )
        )
      ),
      bslib::card(
        bslib::card_header("Data"),
        shiny::uiOutput(ns("data_links")),
        shiny::downloadButton(
          ns("download_bundle"),
          "Processed app bundle (RDS)",
          class = "btn-primary"
        )
      )
    )
  )
}

data_links_ui <- function(bundle) {
  provenance <- bundle$provenance %||% list()
  links <- list()
  bundle_url <- published_bundle_url(bundle)
  data_url <- https_url(provenance_value(
    provenance,
    c("data_url", "data_repository_url")
  ))
  code_url <- https_url(provenance_value(
    provenance,
    c("code_url", "application_code_url")
  ))
  if (is.na(code_url)) {
    code_url <- "https://github.com/carl-stone/ESPIviz"
  }
  doi <- provenance_value(provenance, c("manuscript_doi", "doi"))
  if (!is.na(bundle_url)) {
    links <- c(
      links,
      list(htmltools::tags$li(htmltools::a(
        "Processed app bundle",
        href = bundle_url,
        target = "_blank",
        rel = "noopener noreferrer"
      )))
    )
  }
  if (!is.na(data_url)) {
    links <- c(
      links,
      list(htmltools::tags$li(htmltools::a(
        "Public data repository",
        href = data_url,
        target = "_blank",
        rel = "noopener noreferrer"
      )))
    )
  }
  if (!is.na(code_url)) {
    links <- c(
      links,
      list(htmltools::tags$li(htmltools::a(
        "Application source code",
        href = code_url,
        target = "_blank",
        rel = "noopener noreferrer"
      )))
    )
  }
  if (!is.na(doi)) {
    doi_url <- if (grepl("^10[.][0-9]{4,9}/[^[:space:]]+$", doi)) {
      paste0("https://doi.org/", doi)
    } else {
      https_url(doi)
    }
    if (!is.na(doi_url)) {
      links <- c(
        links,
        list(htmltools::tags$li(htmltools::a(
          "Manuscript",
          href = doi_url,
          target = "_blank",
          rel = "noopener noreferrer"
        )))
      )
    }
  }
  htmltools::div(
    class = "provenance-block",
    htmltools::span(
      class = "data-version",
      paste("Data version", bundle$data_version)
    ),
    htmltools::tags$ul(class = "data-link-list", links),
    htmltools::p(
      if (is.na(doi)) {
        "Manuscript DOI: not yet provided. Cite the study title, Carl Stone's ESPIviz application, and the data version used; a formal manuscript citation will be added when available."
      } else {
        "Use the linked manuscript citation and report the ESPIviz data version used."
      }
    ),
    htmltools::p(
      if (is.na(data_url)) {
        "Final data deposit and accession: pending. The processed app bundle is available from the versioned public release."
      } else {
        "The final data record is linked above."
      }
    ),
    htmltools::p(
      "Reuse notice: final data-deposit terms are pending. This exploration bundle does not grant broader rights than the eventual repository license. Application code is separately MIT licensed."
    )
  )
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
            setdiff(
              names(attributes(clean_bundle)),
              c("bundle_path", "data_manifest")
            )
          ]
          saveRDS(clean_bundle, file, compress = "xz")
        }
      }
    )
  })
}
