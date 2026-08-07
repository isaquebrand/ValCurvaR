#' Fit an analytical calibration curve
#'
#' @param dados Validated data from [validar_curva()] or a data frame with
#'   internal `.x` and `.y` columns.
#' @param metodo One of `"auto"`, `"ols"` or `"wls"`.
#' @param pesos WLS strategy. `"variancia_nivel"` uses the inverse experimental
#'   variance at each concentration; other choices are empirical alternatives.
#' @return A `valcurva_fit` object.
#' @export
ajustar_curva <- function(dados, metodo = c("auto", "ols", "wls"),
                          pesos = c("variancia_nivel", "1/x", "1/x2", "1/y", "1/y2")) {
  metodo <- match.arg(metodo); pesos <- match.arg(pesos)
  if (!all(c(".x", ".y") %in% names(dados))) .valcurva_abort("Use `validar_curva()` antes de ajustar a curva.")
  bf <- .brown_forsythe(dados$.y, dados$.x)
  selected <- if (metodo == "auto") {
    if (isTRUE(bf$disponivel) && is.finite(bf$p_valor) && bf$p_valor < .05) "wls" else "ols"
  } else metodo
  weights <- if (selected == "wls") .make_weights(dados, pesos) else rep(1, nrow(dados))
  model <- if (selected == "wls") {
    stats::lm(.y ~ .x, data = dados, weights = weights)
  } else {
    stats::lm(.y ~ .x, data = dados)
  }
  output <- structure(list(modelo = model, dados = dados, metodo = selected,
                           estrategia_pesos = if (selected == "wls") pesos else "nenhuma",
                           pesos = weights, selecao_automatica = metodo == "auto"),
                      class = "valcurva_fit")
  output$diagnosticos <- diagnosticar_curva(output)
  output
}

.make_weights <- function(data, strategy) {
  x <- data$.x; y <- data$.y
  if (strategy == "variancia_nivel") {
    vars <- tapply(y, x, stats::var)
    if (anyNA(vars) || any(vars <= 0)) .valcurva_abort("Pesos por variancia requerem ao menos duas replicatas e variancia positiva em todos os niveis.")
    return(1 / unname(vars[as.character(x)]))
  }
  base <- switch(strategy, "1/x" = abs(x), "1/x2" = x^2, "1/y" = abs(y), "1/y2" = y^2)
  if (any(base <= 0)) .valcurva_abort("O esquema de pesos escolhido nao aceita concentracao ou sinal igual a zero.")
  1 / base
}

#' @export
print.valcurva_fit <- function(x, ...) {
  co <- stats::coef(x$modelo)
  cat("ValCurvaR", toupper(x$metodo), "fit\n")
  cat(sprintf("sinal = %.8g + %.8g * concentracao\n", co[1], co[2]))
  cat("n =", nrow(x$dados), "| niveis =", length(unique(x$dados$.x)), "\n")
  invisible(x)
}
