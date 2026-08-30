#' Plot an analytical calibration diagram
#'
#' Draws individual data, centroids, fitted line and 95% prediction limits in
#' the visual language of the metrology article supplied for this project.
#'
#' @param ajuste A `valcurva_fit` object.
#' @param nivel Confidence level for prediction limits.
#' @param ... Graphical parameters passed to [graphics::plot()].
#' @return Invisibly, a data frame used to draw the prediction curves.
#' @export
grafico_calibracao <- function(ajuste, nivel = .95, ...) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = graphics::par("mar") + c(0, 0, 0, 3), xpd = NA)
  pal <- .valcurva_palette(); d <- ajuste$dados
  xs <- seq(min(d$.x), max(d$.x), length.out = 200)
  raw_pred <- stats::predict(ajuste$modelo, newdata = data.frame(.x = xs), se.fit = TRUE)
  critical <- stats::qt((1 + nivel) / 2, stats::df.residual(ajuste$modelo))
  half_width <- critical * sqrt(raw_pred$residual.scale^2 + raw_pred$se.fit^2)
  pred <- cbind(fit = raw_pred$fit, lwr = raw_pred$fit - half_width, upr = raw_pred$fit + half_width)
  graphics::plot(d$.x, d$.y, pch = 1, col = pal$data, xlab = "Concentracao", ylab = "Resposta instrumental", ...)
  cent <- stats::aggregate(.y ~ .x, d, mean)
  graphics::points(cent$.x, cent$.y, pch = 3, col = pal$data, lwd = 1.5)
  graphics::lines(xs, pred[, "fit"], col = pal$line, lwd = 2)
  graphics::lines(xs, pred[, "lwr"], col = pal$band, lwd = 1.4, lty = 2)
  graphics::lines(xs, pred[, "upr"], col = pal$band, lwd = 1.4, lty = 2)
  co <- stats::coef(ajuste$modelo)
  graphics::legend("topright", inset = c(.02, 0), bty = "n", cex = .65,
    legend = c("Dados", "Centroide", sprintf("%s: y = %.4g + %.4g x", toupper(ajuste$metodo), co[1], co[2]),
               sprintf("Limites inferior e superior de predicao %.0f%% (tracejados)", 100 * nivel)),
    col = c(pal$data, pal$data, pal$line, pal$band), pch = c(1, 3, NA, NA),
    lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 1.4))
  invisible(data.frame(concentracao = xs, ajuste = pred[, "fit"], inferior = pred[, "lwr"], superior = pred[, "upr"]))
}

#' Plot instrumental-response variability by concentration
#' @param ajuste A `valcurva_fit` object.
#' @return Invisibly, the level summary.
#' @export
grafico_variancia <- function(ajuste) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = graphics::par("mar") + c(0, 0, 0, 3), xpd = NA)
  s <- ajuste$diagnosticos$resumo_niveis; pal <- .valcurva_palette()
  ymax <- max(s$dp_sinal, na.rm = TRUE)
  graphics::plot(s$concentracao, s$dp_sinal, pch = 21, bg = "white", col = pal$data,
                 ylim = c(0, ymax * 1.2 + .Machine$double.eps),
                 xlab = "Concentracao", ylab = "Desvio-padrao do sinal")
  graphics::arrows(s$concentracao, pmax(0, s$dp_sinal - 0.5 * s$dp_sinal), s$concentracao,
                   s$dp_sinal + 0.5 * s$dp_sinal, angle = 90, code = 3, length = .04, col = "#666666")
  if (sum(is.finite(s$dp_sinal)) >= 2L) graphics::abline(stats::lm(dp_sinal ~ concentracao, s), col = pal$data, lwd = 2)
  graphics::legend("topright", inset = c(.02, 0), legend = c("DP por nivel", "Tendencia linear"),
                   pch = c(21, NA), lty = c(NA, 1), col = c(pal$data, pal$data), bty = "n", cex = .8)
  invisible(s)
}

#' Plot standardized residuals against concentration
#' @param ajuste A `valcurva_fit` object.
#' @return Invisibly, residual diagnostic data.
#' @export
grafico_residuos <- function(ajuste) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = graphics::par("mar") + c(0, 0, 0, 3), xpd = NA)
  z <- ajuste$diagnosticos$influencia; pal <- .valcurva_palette()
  graphics::plot(z$concentracao, z$residuo_padronizado, pch = 19, col = pal$data,
                 xlab = "Concentracao", ylab = "Residuo padronizado")
  graphics::abline(h = 0, lty = 2, col = pal$zero)
  graphics::abline(h = c(-2, 2), lty = 3, col = "#888888")
  if (length(unique(z$concentracao)) >= 4L) graphics::lines(stats::lowess(z$concentracao, z$residuo_padronizado), col = pal$data, lty = 3, lwd = 2)
  graphics::legend("topright", inset = c(.02, 0), legend = c("Residuos", "Zero", "Limites +/- 2", "Suavizacao"),
                   pch = c(19, NA, NA, NA), lty = c(NA, 2, 3, 3), col = c(pal$data, pal$zero, "#888888", pal$data), bty = "n", cex = .8)
  invisible(z)
}

#' Plot a normal Q-Q graph with a simulation envelope for residuals
#'
#' @param ajuste A `valcurva_fit` object.
#' @param nivel Envelope coverage level.
#' @param simulacoes Number of normal samples used to build the envelope.
#' @param semente Optional seed for a reproducible envelope.
#' @return Invisibly, data used to construct the graph.
#' @export
grafico_qq <- function(ajuste, nivel = .95, simulacoes = 999L, semente = NULL) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = graphics::par("mar") + c(0, 0, 0, 3), xpd = NA)
  if (!is.null(semente)) set.seed(semente)
  r <- sort(stats::rstandard(ajuste$modelo)); n <- length(r)
  if (n < 3L) .valcurva_abort("O grafico Q-Q requer ao menos tres residuos.")
  p <- stats::ppoints(n); teorico <- stats::qnorm(p)
  sim <- replicate(simulacoes, sort(stats::rnorm(n)))
  alpha <- (1 - nivel) / 2
  env <- t(apply(sim, 1, stats::quantile, probs = c(alpha, 1 - alpha), names = FALSE))
  graphics::plot(teorico, r, pch = 19, col = .valcurva_palette()$data,
                 xlab = "Quantis teoricos normais", ylab = "Residuos padronizados")
  graphics::lines(teorico, env[, 1], lty = 2, col = .valcurva_palette()$band)
  graphics::lines(teorico, env[, 2], lty = 2, col = .valcurva_palette()$band)
  graphics::abline(stats::lm(r ~ teorico), col = .valcurva_palette()$line, lwd = 2)
  graphics::legend("topright", inset = c(.02, 0), legend = c("Residuos", "Envelope", "Reta"),
                   pch = c(19, NA, NA), lty = c(NA, 2, 1), col = c(.valcurva_palette()$data, .valcurva_palette()$band, .valcurva_palette()$line), bty = "n", cex = .8)
  invisible(data.frame(quantil_teorico = teorico, residuo = r, envelope_inferior = env[, 1], envelope_superior = env[, 2]))
}

#' Plot Cook's distance and leverage against their screening limits
#' @param ajuste A `valcurva_fit` object.
#' @return Invisibly, the influence diagnostics.
#' @export
grafico_influencia <- function(ajuste) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  z <- ajuste$diagnosticos$influencia; n <- nrow(z); p <- length(stats::coef(ajuste$modelo))
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(5, 4, 4, 5), xpd = NA)
  graphics::plot(seq_len(n), z$cook, pch = 19, xlab = "Observacao", ylab = "Distancia de Cook")
  graphics::abline(h = 4 / n, lty = 2, col = .valcurva_palette()$band)
  graphics::legend("topright", inset = c(-.07, 0), legend = c("Cook", "Limite 4/n"), pch = c(19, NA), lty = c(NA, 2), col = c("black", .valcurva_palette()$band), bty = "n", cex = .7)
  graphics::plot(seq_len(n), z$alavancagem, pch = 19, xlab = "Observacao", ylab = "Alavancagem")
  graphics::abline(h = 2 * p / n, lty = 2, col = .valcurva_palette()$band)
  graphics::legend("topright", inset = c(-.07, 0), legend = c("Alavancagem", "Limite 2p/n"), pch = c(19, NA), lty = c(NA, 2), col = c("black", .valcurva_palette()$band), bty = "n", cex = .7)
  invisible(z)
}

#' Plot inverse prediction uncertainty over the calibration range
#' @param ajuste A `valcurva_fit` object.
#' @param k Numeric vector of numbers of replicate measurements.
#' @return Invisibly, plotted uncertainty data.
#' @export
grafico_incerteza_predicao <- function(ajuste, k = c(1, 2, 3, 4, 9, 16)) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = graphics::par("mar") + c(0, 0, 0, 4), xpd = NA)
  d <- ajuste$dados; xs <- seq(min(d$.x), max(d$.x), length.out = 200); pal <- .valcurva_palette()
  cols <- grDevices::hcl.colors(length(k), "Dynamic")
  all <- do.call(rbind, lapply(k, function(ki) {
    co <- stats::coef(ajuste$modelo); n <- nrow(d); sxx <- sum((d$.x - mean(d$.x))^2)
    u <- summary(ajuste$modelo)$sigma / abs(co[2]) * sqrt(1 / ki + 1 / n + (xs - mean(d$.x))^2 / sxx)
    data.frame(concentracao = xs, k = ki, u = u)
  }))
  ylim <- range(all$u)
  graphics::plot(NA, xlim = range(xs), ylim = ylim, xlab = "Concentracao", ylab = "Incerteza padrao da concentracao")
  for (i in seq_along(k)) {
    z <- all[all$k == k[i], ]; graphics::lines(z$concentracao, z$u, col = cols[i], lwd = 2)
  }
  graphics::legend("topright", inset = c(.02, 0), legend = paste0("K=", k),
                   col = cols, lty = 1, lwd = 2, bty = "n", xpd = NA, cex = .8)
  invisible(all)
}

#' Create the six-panel calibration figure
#' @param ajuste A `valcurva_fit` object.
#' @param k Replicate counts for the uncertainty panel.
#' @return Invisibly, `ajuste`.
#' @export
painel_calibracao <- function(ajuste, k = c(1, 2, 3, 4, 9, 16)) {
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(3, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 1, 0))
  grafico_calibracao(ajuste, main = sprintf("Curva de calibracao (%s)", toupper(ajuste$metodo)))
  grafico_variancia(ajuste); graphics::title("Variancia por nivel")
  grafico_residuos(ajuste); graphics::title("Diagnostico de residuos")
  grafico_qq(ajuste, simulacoes = 399L, semente = 1); graphics::title("Q-Q com envelope normal")
  grafico_incerteza_predicao(ajuste, k = k); graphics::title("Incerteza de predicao")
  graphics::plot(seq_len(nrow(ajuste$dados)), ajuste$diagnosticos$influencia$cook, pch = 19,
                 xlab = "Observacao", ylab = "Distancia de Cook")
  graphics::abline(h = 4 / nrow(ajuste$dados), lty = 2, col = .valcurva_palette()$band)
  graphics::title("Influencia: Distancia de Cook")
  invisible(ajuste)
}
