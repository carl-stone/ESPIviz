test_that("gene lookup is case-insensitive and preserves canonical names", {
  expect_app_helper("parse_gene_input")
  universe <- c("Glul", "EGFP", "Mcm2", "ZeroGene")

  parsed <- parse_gene_input(
    "glul, EGFP\nGLUL; mcm2  not_a_gene",
    universe = universe
  )

  expect_equal(parsed$genes, c("Glul", "EGFP", "Mcm2"))
  expect_equal(parsed$missing, "not_a_gene")
  expect_length(parsed$truncated, 0L)
})

test_that("gene parsing removes duplicates without imposing an app-sized cap", {
  expect_app_helper("parse_gene_input")
  universe <- c("EGFP", sprintf("Gene%03d", seq_len(80)))
  requested <- c(universe, "gene001", "EGFP")

  parsed <- parse_gene_input(
    paste(requested, collapse = "\n"),
    universe = universe,
    max_genes = Inf
  )

  expect_equal(parsed$genes, universe)
  expect_length(parsed$genes, 81L)
  expect_length(parsed$truncated, 0L)
})

test_that("parser accepts the complete 38,394-gene universe", {
  expect_app_helper("parse_gene_input")
  universe <- c("EGFP", sprintf("Gene%05d", seq_len(38393L)))

  parsed <- parse_gene_input(
    paste(universe, collapse = "\n"),
    universe = universe,
    max_genes = Inf
  )

  expect_length(parsed$genes, 38394L)
  expect_identical(parsed$genes, universe)
  expect_length(parsed$missing, 0L)
  expect_length(parsed$truncated, 0L)
})

test_that("missing and empty gene input do not block valid genes", {
  expect_app_helper("parse_gene_input")

  mixed <- parse_gene_input("missing, egfp", c("EGFP"))
  empty <- parse_gene_input("  \n, ;", c("EGFP"))

  expect_equal(mixed$genes, "EGFP")
  expect_equal(mixed$missing, "missing")
  expect_length(empty$genes, 0L)
  expect_length(empty$missing, 0L)
})

test_that("visual gene comparisons paginate at 25 without dropping genes", {
  expect_app_helper("paginate_genes")
  genes <- sprintf("Gene%03d", seq_len(63))

  first <- paginate_genes(genes, page = 1L)
  second <- paginate_genes(genes, page = 2L, page_size = 25L)
  last <- paginate_genes(genes, page = 3L, page_size = 25L)

  expect_equal(first$genes, genes[1:25])
  expect_equal(second$genes, genes[26:50])
  expect_equal(last$genes, genes[51:63])
  expect_equal(c(first$total, second$total, last$total), rep(63L, 3))
  expect_equal(c(first$pages, second$pages, last$pages), rep(3L, 3))
  expect_equal(unlist(list(first$genes, second$genes, last$genes)), genes)
})

test_that("an explicit parser cap reports rather than silently loses genes", {
  expect_app_helper("parse_gene_input")
  universe <- sprintf("Gene%03d", seq_len(30))

  parsed <- parse_gene_input(
    paste(universe, collapse = ","),
    universe = universe,
    max_genes = 25L
  )

  expect_equal(parsed$genes, universe[1:25])
  expect_equal(parsed$truncated, universe[26:30])
})
