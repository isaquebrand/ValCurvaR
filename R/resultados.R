#' Compare OLS and WLS calibration models
#'
#' @param ols An OLS `valcurva_fit` object.
#' @param wls A WLS `valcurva_fit` object.
#' @return A data frame intended as decision support, not an automatic approval.
#' @export
comparar_modelos <- function(ols, wls) {
  if (!inherits(ols, "valcurva_fit") || !inherits(wls, "valcurva_fit"))
    .valcurva_abort("Forneca dois objetos produzidos por `ajustar_curva()`.")
  make_row <- function(x) {
    fitted_x <- stats::fitted(x$modelo)
    back <- (x$dados$.y - stats::coef(x$modelo)[1]) / stats::coef(x$modelo)[2]
    data.frame(metodo = toupper(x$metodo), pesos = x$estrategia_pesos,
               intercepto = stats::coef(x$modelo)[1], inclinacao = stats::coef(x$modelo)[2],
               r2 = summary(x$modelo)$r.squared,
               rmse_sinal = sqrt(mean((x$dados$.y - fitted_x)^2)),
               vies_medio_retrocalculado = mean(back - x$dados$.x),
               erro_relativo_medio_percentual = mean(abs((back - x$dados$.x) / x$dados$.x)) * 100)
  }
  rbind(make_row(ols), make_row(wls))
}

#' Back-calculate concentration with GUM uncertainty propagation
#'
#' The measurement model is `x = (y - a) / b`. The covariance of the fitted
#' intercept (`a`) and slope (`b`) is retained explicitly in the combined
#' uncertainty. It must not be discarded or treated as two independent inputs.
#'
#' @param ajuste A `valcurva_fit` object.
#' @param sinal Numeric vector of instrumental responses.
#' @param k Number of replicate signals averaged for each result.
#' @param u_sinal Standard uncertainty of each reported signal. If `NULL`, the
#'   residual standard deviation of the fit divided by `sqrt(k)` is used as a
#'   data-derived estimate.
#' @param nivel Confidence level.
#' @param fator_cobertura Optional coverage factor. If `NULL`, a two-sided
#'   Student t factor using the regression degrees of freedom is used. This
#'   keeps the correlated regression coefficients as one uncertainty block.
#' @return A data frame with concentration, combined and expanded GUM
#'   uncertainty, and the individual variance contributions.
#' @export
retrocalcular_concentracao <- function(ajuste, sinal, k = 1L, u_sinal = NULL,
                                       nivel = .95, fator_cobertura = NULL) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  if (!is.numeric(sinal) || any(!is.finite(sinal)) || length(k) != 1L || k < 1L)
    .valcurva_abort("`sinal` deve ser finito e `k` deve ser maior ou igual a um.")
  if (is.null(u_sinal)) u_sinal <- summary(ajuste$modelo)$sigma / sqrt(k)
  if (!is.numeric(u_sinal) || any(!is.finite(u_sinal)) || any(u_sinal < 0))
    .valcurva_abort("`u_sinal` deve conter incertezas padrao finitas e nao negativas.")
  u_sinal <- rep(u_sinal, length.out = length(sinal))
  co <- stats::coef(ajuste$modelo); a <- co[1]; b <- co[2]
  if (b == 0) .valcurva_abort("A inclinacao da curva nao pode ser zero.")
  xhat <- (sinal - a) / b
  V <- stats::vcov(ajuste$modelo)
  ca <- rep(-1 / b, length(sinal)); cb <- -xhat / b; cy <- rep(1 / b, length(sinal))
  var_sinal <- (cy * u_sinal)^2
  var_intercepto <- (ca^2) * V[1, 1]
  var_inclinacao <- (cb^2) * V[2, 2]
  cov_intercepto_inclinacao <- 2 * ca * cb * V[1, 2]
  u2 <- var_sinal + var_intercepto + var_inclinacao + cov_intercepto_inclinacao
  if (any(u2 < -sqrt(.Machine$double.eps))) .valcurva_abort("A variancia combinada calculada foi negativa; verifique as entradas.")
  u <- sqrt(pmax(u2, 0))
  gl_regressao <- stats::df.residual(ajuste$modelo)
  fator <- if (is.null(fator_cobertura)) stats::qt((1 + nivel) / 2, gl_regressao) else fator_cobertura
  if (!is.numeric(fator) || length(fator) != 1L || !is.finite(fator) || fator <= 0)
    .valcurva_abort("`fator_cobertura` deve ser um numero positivo.")
  data.frame(
    sinal = sinal, concentracao = xhat, u_sinal = u_sinal,
    u_padrao = u, fator_cobertura = fator, U_expandida = fator * u,
    nivel = nivel, gl_regressao = gl_regressao,
    contribuicao_sinal = var_sinal, contribuicao_intercepto = var_intercepto,
    contribuicao_inclinacao = var_inclinacao,
    contribuicao_cov_intercepto_inclinacao = cov_intercepto_inclinacao
  )
}
