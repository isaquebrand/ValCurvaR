#' Evaluate the impact of suspicious observations without deleting them
#'
#' Fits leave-one-out models solely for technical review. It never changes the
#' original data or approves exclusion of a measurement.
#'
#' @param ajuste A `valcurva_fit` object.
#' @param observacoes Optional row numbers. By default, all observations flagged
#'   by Cook's distance, leverage, DFFITS, DFBETAS or studentized residuals are used.
#' @return A data frame comparing each leave-one-out model with the original fit.
#' @export
analisar_sensibilidade <- function(ajuste, observacoes = NULL) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  inf <- ajuste$diagnosticos$influencia
  if (is.null(observacoes)) observacoes <- which(inf$flag_cook | inf$flag_alavancagem | inf$flag_dffits | inf$flag_dfbeta | inf$flag_studentizado)
  observacoes <- sort(unique(as.integer(observacoes)))
  if (!length(observacoes)) return(data.frame())
  if (anyNA(observacoes) || any(observacoes < 1L | observacoes > nrow(ajuste$dados))) .valcurva_abort("`observacoes` deve conter indices de linhas validos.")
  full <- stats::coef(ajuste$modelo); r2 <- summary(ajuste$modelo)$r.squared
  do.call(rbind, lapply(observacoes, function(i) {
    d <- ajuste$dados[-i, , drop = FALSE]
    m <- if (ajuste$metodo == "wls") stats::lm(.y ~ .x, d, weights = .make_weights(d, ajuste$estrategia_pesos)) else stats::lm(.y ~ .x, d)
    co <- stats::coef(m)
    data.frame(observacao = i, concentracao = ajuste$dados$.x[i], sinal = ajuste$dados$.y[i],
               intercepto_sem_observacao = co[1], inclinacao_sem_observacao = co[2],
               mudanca_intercepto_percentual = 100 * (co[1] - full[1]) / ifelse(abs(full[1]) < sqrt(.Machine$double.eps), 1, abs(full[1])),
               mudanca_inclinacao_percentual = 100 * (co[2] - full[2]) / abs(full[2]),
               mudanca_r2 = summary(m)$r.squared - r2, row.names = NULL)
  }))
}
