.valcurva_abort <- function(message) stop(message, call. = FALSE)

.valcurva_name <- function(data, expression) {
  if (is.symbol(expression)) {
    name <- as.character(expression)
    if (name %in% names(data)) return(name)
  }
  if (is.character(expression) && length(expression) == 1L && expression %in% names(data)) return(expression)
  .valcurva_abort("Informe uma coluna existente por nome ou sem aspas.")
}

.valcurva_palette <- function() {
  list(data = "#1F5FBF", line = "#111111", band = "#D62CD1", grid = "#D9D9D9",
       wls = "#E67E22", zero = "#6C8EBF")
}

.has_replicates <- function(data) any(table(data$.x) > 1L)
