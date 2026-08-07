#' Validate analytical calibration data
#'
#' Keeps the original measurements in long format. Replicates must not be
#' averaged before validation because they are needed for variance and
#' lack-of-fit diagnostics.
#'
#' @param dados A data frame.
#' @param concentracao Concentration column, supplied bare or as a string.
#' @param sinal Instrumental response column, supplied bare or as a string.
#' @param replica Optional replicate identifier column.
#' @return A `data.frame` with standard internal columns `.x`, `.y` and
#'   `.replica`.
#' @export
validar_curva <- function(dados, concentracao, sinal, replica = NULL) {
  if (!is.data.frame(dados)) .valcurva_abort("`dados` deve ser um data frame.")
  x_name <- .valcurva_name(dados, substitute(concentracao))
  y_name <- .valcurva_name(dados, substitute(sinal))
  x <- dados[[x_name]]
  y <- dados[[y_name]]
  if (!is.numeric(x) || !is.numeric(y)) .valcurva_abort("Concentracao e sinal devem ser numericos.")
  if (anyNA(x) || anyNA(y) || any(!is.finite(x)) || any(!is.finite(y))) {
    .valcurva_abort("Concentracao e sinal devem ser finitos e sem valores ausentes.")
  }
  if (length(unique(x)) < 3L) .valcurva_abort("Sao necessarios ao menos tres niveis de concentracao.")
  replica_expr <- substitute(replica)
  replica_name <- if (identical(replica_expr, quote(NULL))) NULL else .valcurva_name(dados, replica_expr)
  result <- data.frame(.x = x, .y = y, .replica = if (is.null(replica_name)) seq_along(x) else dados[[replica_name]])
  result <- result[order(result$.x, result$.replica), , drop = FALSE]
  attr(result, "valcurva_columns") <- c(concentracao = x_name, sinal = y_name, replica = replica_name)
  result
}
