test_that("programmatic gene changes keep server selectize data valid", {
  bundle <- synthetic_bundle()
  state <- new_app_state(bundle)
  registered <- list()
  mock_session <- shiny::MockShinySession$new()
  mock_session$registerDataObj <- function(name, data, filterFunc) {
    registered[[length(registered) + 1L]] <<- list(
      name = name,
      data = data,
      filter = filterFunc
    )
    paste0("mock-data-object-", length(registered))
  }

  shiny::testServer(
    explore_server,
    args = list(bundle = bundle, state = state),
    session = mock_session,
    {
      session$flushReact()
      set_state_gene(state, bundle, "EGFP")
      session$flushReact()
    }
  )

  gene_objects <- Filter(
    function(object) endsWith(object$name, "active_gene"),
    registered
  )
  secondary_gene_objects <- Filter(
    function(object) endsWith(object$name, "secondary_gene"),
    registered
  )
  expect_gt(length(gene_objects), 0L)
  expect_gt(length(secondary_gene_objects), 0L)

  request <- list(
    QUERY_STRING = paste0(
      "query=&field=%5B%5B%22label%22%5D%5D&value=value",
      "&conju=and&maxop=1000"
    )
  )
  for (object in c(gene_objects, secondary_gene_objects)) {
    expect_false(is.null(object$data))
    expect_true(all(c("label", "value") %in% names(object$data)))
    expect_silent(object$filter(object$data, request))
  }
})
