#' Validate analytical calibration data
#'
#' Keeps the original measurements in long format. A horizontal table is also
#' accepted: omit `sinal` and every numeric column other than `concentracao`
#' is converted internally to a replicate column. Replicates must not be
#' averaged before validation because they are needed for variance and
#' lack-of-fit diagnostics.
#'
#' @param dados A data frame.
#' @param concentracao Concentration column, supplied bare or as a string.
#' @param sinal Instrumental response column, supplied bare or as a string.
#'   Omit it for horizontal data with one numeric response column per replicate.
#' @param replica Optional replicate identifier column.
#' @return A `data.frame` with standard internal columns `.x`, `.y` and
#'   `.replica`.
#' @export
validar_curva <- function(dados, concentracao, sinal = NULL, replica = NULL) {
  if (!is.data.frame(dados)) .valcurva_abort("`dados` deve ser um data frame.")
  x_name <- .valcurva_name(dados, substitute(concentracao))
  sinal_expr <- substitute(sinal)
  if (missing(sinal) || identical(sinal_expr, quote(NULL))) {
    return(.validar_tabela_horizontal(dados, x_name))
  }
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

.validar_tabela_horizontal <- function(dados, x_name) {
  x <- dados[[x_name]]
  if (!is.numeric(x)) .valcurva_abort("Concentracao deve ser numerica.")
  colunas_resposta <- setdiff(names(dados), x_name)
  colunas_resposta <- colunas_resposta[vapply(dados[colunas_resposta], is.numeric, logical(1))]
  if (!length(colunas_resposta)) .valcurva_abort("A tabela horizontal precisa de ao menos uma coluna numerica de resposta alem da concentracao.")
  longo <- do.call(rbind, lapply(colunas_resposta, function(nome) {
    data.frame(.x = x, .y = dados[[nome]], .replica = nome, stringsAsFactors = FALSE)
  }))
  if (anyNA(longo$.x) || anyNA(longo$.y) || any(!is.finite(longo$.x)) || any(!is.finite(longo$.y)))
    .valcurva_abort("Concentracao e sinais devem ser finitos e sem valores ausentes.")
  if (length(unique(longo$.x)) < 3L) .valcurva_abort("Sao necessarios ao menos tres niveis de concentracao.")
  longo <- longo[order(longo$.x, longo$.replica), , drop = FALSE]
  attr(longo, "valcurva_columns") <- c(concentracao = x_name, sinal = paste(colunas_resposta, collapse = ", "), replica = "nome da coluna")
  longo
}
