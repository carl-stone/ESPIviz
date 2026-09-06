`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

compact_character <- function(x) {
  x <- as.character(x %||% character())
  x[!is.na(x) & nzchar(x)]
}

casefold_key <- function(x) {
  tolower(enc2utf8(as.character(x)))
}

condition_label <- function(x) {
  labels <- c(
    "Control" = "p27CKO",
    "control" = "p27CKO",
    "E-Stim" = "p27CKO + E-Stim",
    "estim" = "p27CKO + E-Stim",
    "p27CKOEStim" = "p27CKO + E-Stim",
    "p27CKO +EStim" = "p27CKO + E-Stim"
  )
  values <- as.character(x)
  matched <- match(values, names(labels))
  values[!is.na(matched)] <- unname(labels[matched[!is.na(matched)]])
  values
}

sample_label <- function(x) {
  x <- sub("_p27CKOEStim$", " · p27CKO + E-Stim", as.character(x))
  sub("_p27CKO$", " · p27CKO", x)
}

label_result_conditions <- function(data) {
  for (column in intersect(c("condition", "direction"), names(data))) {
    data[[column]] <- condition_label(data[[column]])
  }
  if ("contrast" %in% names(data)) {
    data$contrast[data$contrast == "estim_vs_control"] <-
      "p27CKO + E-Stim vs p27CKO"
  }
  data
}

safe_ratio <- function(numerator, denominator) {
  ifelse(
    is.finite(denominator) & denominator != 0,
    numerator / denominator,
    NA_real_
  )
}

format_number <- function(x, digits = 3L) {
  ifelse(
    is.na(x) | !is.finite(x),
    "—",
    formatC(x, digits = digits, format = "fg", flag = "#")
  )
}

format_percent <- function(x, digits = 1L) {
  ifelse(
    is.na(x) | !is.finite(x),
    "—",
    paste0(formatC(x, digits = digits, format = "f"), "%")
  )
}

as_named_palette <- function(x, fallback) {
  values <- unlist(x %||% list(), use.names = TRUE)
  values <- values[nzchar(values)]
  if (length(values) == 0L) fallback else values
}

safe_filename <- function(x, fallback = "espiviz") {
  x <- gsub("[^A-Za-z0-9._-]+", "-", as.character(x %||% ""))
  x <- gsub("(^-+|-+$)", "", x)
  if (!nzchar(x)) fallback else x
}
