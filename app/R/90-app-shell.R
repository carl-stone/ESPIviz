app_theme <- function() {
  font <- bslib::font_collection(
    "-apple-system",
    "BlinkMacSystemFont",
    "Segoe UI",
    "Helvetica Neue",
    "Arial",
    "sans-serif"
  )
  bslib::bs_theme(
    version = 5,
    bg = "#f1f3f4",
    fg = "#20282c",
    primary = "#14764f",
    secondary = "#58646a",
    base_font = font,
    heading_font = font,
    `border-radius` = "0.3rem"
  )
}

app_assets <- function() {
  htmltools::tagList(
    htmltools::tags$link(
      rel = "icon",
      type = "image/svg+xml",
      href = "favicon.svg"
    ),
    htmltools::tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = "styles.css"
    ),
    htmltools::tags$script(src = "interactions.js", defer = NA)
  )
}

maintenance_ui <- function() {
  bslib::page_fluid(
    lang = "en",
    title = "ESPIviz",
    theme = app_theme(),
    app_assets(),
    htmltools::tags$main(
      class = "maintenance-shell",
      htmltools::div(
        class = "maintenance-card",
        htmltools::p("ESPIviz", class = "maintenance-brand"),
        htmltools::h1("The explorer is temporarily unavailable."),
        htmltools::p("The application data could not be loaded."),
        htmltools::a("Reload", href = "", class = "btn btn-primary"),
        htmltools::a(
          "Public source and data",
          href = "https://github.com/carl-stone/ESPIviz",
          class = "maintenance-link"
        )
      )
    )
  )
}

app_ui <- function(bundle) {
  bslib::page_fluid(
    lang = "en",
    title = "ESPIviz single-cell explorer",
    theme = app_theme(),
    app_assets(),
    htmltools::tags$meta(
      name = "description",
      content = "Interactive explorer for the ESPI single-cell RNA-seq study."
    ),
    htmltools::a(
      "Skip to analysis",
      href = "#analysis-workspace",
      class = "skip-link"
    ),
    htmltools::tags$header(
      class = "app-header",
      htmltools::div(
        class = "brand-name",
        htmltools::HTML('ESPI<span class="brand-accent">viz</span>')
      ),
      htmltools::span("Single-cell explorer", class = "brand-subtitle")
    ),
    htmltools::div(
      class = "app-nav-layout",
      id = "analysis-workspace",
      bslib::navset_pill_list(
        id = "main_nav",
        selected = "Explore",
        well = FALSE,
        widths = c(2, 10),
        bslib::nav_panel(
          "Explore",
          explore_ui("explore", bundle),
          value = "Explore"
        ),
        bslib::nav_panel(
          "Differential expression",
          differential_expression_ui("de"),
          value = "Differential expression"
        ),
        bslib::nav_panel(
          "Pathways",
          pathways_ui("pathways"),
          value = "Pathways"
        ),
        bslib::nav_panel("About", about_ui("about"), value = "About")
      )
    )
  )
}

app_view_label <- function(value) {
  view_map <- c(
    explore = "Explore",
    de = "Differential expression",
    differential_expression = "Differential expression",
    pathways = "Pathways",
    about = "About"
  )
  key <- tolower(as.character(value %||% ""))
  view <- unname(view_map[key])
  if (length(view) == 0L || is.na(view)) {
    return(NULL)
  }
  view
}

app_view_slug <- function(label) {
  switch(
    label %||% "Explore",
    "Differential expression" = "de",
    "Pathways" = "pathways",
    "About" = "about",
    "explore"
  )
}

app_server <- function(bundle) {
  force(bundle)
  function(input, output, session) {
    state <- new_app_state(bundle)
    url_initialized <- shiny::reactiveVal(FALSE)
    pending_view <- shiny::reactiveVal(NULL)
    navigate_explore <- function() {
      bslib::nav_select("main_nav", "Explore", session = session)
    }

    explore_server("explore", bundle, state)
    differential_expression_server("de", bundle, state, navigate_explore)
    pathways_server("pathways", bundle, state, navigate_explore)
    about_server("about", bundle)

    shiny::observeEvent(
      session$clientData$url_search,
      {
        query <- shiny::parseQueryString(session$clientData$url_search %||% "")
        if (!is.null(query$view_state)) {
          tryCatch(
            {
              value <- jsonlite::fromJSON(
                query$view_state,
                simplifyVector = TRUE
              )
              restored <- restore_view_state(value, bundle, state)
              pending_view(app_view_label(restored$view) %||% "Explore")
            },
            error = function(error) {
              shiny::showNotification(
                conditionMessage(error),
                type = "error",
                duration = NULL
              )
            }
          )
        }
        view <- app_view_label(query$view)
        if (!is.null(view)) {
          pending_view(view)
          bslib::nav_select("main_nav", view, session = session)
        }
        if (!is.null(query$gene)) {
          set_state_gene(state, bundle, query$gene)
        }
        if (!is.null(query$genes)) {
          replace_state_gene_set(state, bundle, query$genes)
          if (!is.null(query$gene)) set_state_gene(state, bundle, query$gene)
        }
        if (
          !is.null(query$pathway) &&
            query$pathway %in% bundle$pathways$pathway_id
        ) {
          state$active_pathway(query$pathway)
        }
        url_initialized(TRUE)
      },
      once = TRUE
    )

    shiny::observe({
      shiny::req(url_initialized())
      current_view <- input$main_nav %||% "Explore"
      initial_view <- pending_view()
      if (!is.null(initial_view) && !identical(current_view, initial_view)) {
        return()
      }
      if (!is.null(initial_view)) {
        pending_view(NULL)
      }
      view_slug <- app_view_slug(current_view)
      query <- paste0(
        "?view=",
        utils::URLencode(view_slug, reserved = TRUE),
        "&gene=",
        utils::URLencode(state$active_gene(), reserved = TRUE)
      )
      if (
        identical(view_slug, "pathways") && !is.null(state$active_pathway())
      ) {
        query <- paste0(
          query,
          "&pathway=",
          utils::URLencode(state$active_pathway(), reserved = TRUE)
        )
      }
      shiny::updateQueryString(query, mode = "replace", session = session)
    })
  }
}
