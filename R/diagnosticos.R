#' Calculate calibration-curve diagnostics
#'
#' @param ajuste An object returned by [ajustar_curva()].
#' @return A named list of diagnostic results.
#' @export
diagnosticar_curva <- function(ajuste) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  data <- ajuste$dados
  fit <- ajuste$modelo
  n <- nrow(data); k <- length(unique(data$.x))
  residuals <- stats::residuals(fit)
  h <- stats::hatvalues(fit)
  sigma <- sqrt(sum(residuals^2) / stats::df.residual(fit))
  standardized <- residuals / (sigma * sqrt(pmax(1 - h, .Machine$double.eps)))
  cooks <- stats::cooks.distance(fit)
  bf <- .brown_forsythe(data$.y, data$.x)
  lof <- .lack_of_fit(fit, data)
  quadratic <- .mandel_test(fit, data)
  list(
    resumo_niveis = .level_summary(data),
    homocedasticidade = bf,
    falta_ajuste = lof,
    mandel = quadratic,
    influencia = data.frame(concentracao = data$.x, sinal = data$.y,
                             residuo = residuals, residuo_padronizado = standardized,
                             alavancagem = h, cook = cooks),
    normalidade = if (length(residuals) >= 3L && length(residuals) <= 5000L)
      stats::shapiro.test(residuals) else NULL,
    avisos = .diagnostic_messages(bf, lof, quadratic, cooks, n)
  )
}

.level_summary <- function(data) {
  split_y <- split(data$.y, data$.x)
  data.frame(concentracao = as.numeric(names(split_y)), n = lengths(split_y),
             media_sinal = vapply(split_y, mean, numeric(1)),
             dp_sinal = vapply(split_y, stats::sd, numeric(1)),
             cv_percentual = 100 * vapply(split_y, stats::sd, numeric(1)) /
               vapply(split_y, mean, numeric(1)), row.names = NULL)
}

.brown_forsythe <- function(y, x) {
  group <- as.factor(x)
  if (any(table(group) < 2L)) return(list(disponivel = FALSE, p_valor = NA_real_,
                                          mensagem = "Sao necessarias ao menos duas replicatas por nivel."))
  z <- abs(y - stats::ave(y, group, FUN = stats::median))
  fit <- stats::lm(z ~ group)
  tab <- stats::anova(fit)
  list(disponivel = TRUE, estatistica = unname(tab$`F value`[1]),
       gl1 = unname(tab$Df[1]), gl2 = unname(tab$Df[2]), p_valor = unname(tab$`Pr(>F)`[1]))
}

.lack_of_fit <- function(fit, data) {
  counts <- table(data$.x)
  if (!any(counts > 1L)) return(list(disponivel = FALSE, p_valor = NA_real_,
                                      mensagem = "Falta de ajuste requer replicatas."))
  k <- length(counts); n <- nrow(data)
  ss_pe <- sum(vapply(split(data$.y, data$.x), function(z) sum((z - mean(z))^2), numeric(1)))
  ss_res <- sum(stats::residuals(fit)^2)
  ss_lof <- max(0, ss_res - ss_pe)
  df_lof <- k - 2L; df_pe <- n - k
  if (df_lof <= 0L || df_pe <= 0L || ss_pe == 0) return(list(disponivel = FALSE, p_valor = NA_real_,
    mensagem = "Graus de liberdade ou erro puro insuficientes para falta de ajuste."))
  f <- (ss_lof / df_lof) / (ss_pe / df_pe)
  list(disponivel = TRUE, ss_erro_puro = ss_pe, ss_falta_ajuste = ss_lof,
       gl_falta_ajuste = df_lof, gl_erro_puro = df_pe, estatistica = f,
       p_valor = stats::pf(f, df_lof, df_pe, lower.tail = FALSE))
}

.mandel_test <- function(fit, data) {
  if (length(unique(data$.x)) < 4L) return(list(disponivel = FALSE, p_valor = NA_real_))
  quad <- stats::lm(.y ~ .x + I(.x^2), data = data, weights = fit$weights)
  tab <- stats::anova(fit, quad)
  list(disponivel = TRUE, estatistica = unname(tab$F[2]), p_valor = unname(tab$`Pr(>F)`[2]),
       modelo_quadratico = quad)
}

.diagnostic_messages <- function(bf, lof, mandel, cooks, n) {
  result <- character()
  if (isTRUE(bf$disponivel) && is.finite(bf$p_valor) && bf$p_valor < .05)
    result <- c(result, "Ha evidencia de heterocedasticidade; compare OLS e WLS.")
  if (isTRUE(lof$disponivel) && is.finite(lof$p_valor) && lof$p_valor < .05)
    result <- c(result, "Ha evidencia de falta de ajuste linear; investigue a faixa ou a curvatura.")
  if (isTRUE(mandel$disponivel) && is.finite(mandel$p_valor) && mandel$p_valor < .05)
    result <- c(result, "O termo quadratico melhora o ajuste; nao aprove a linearidade sem investigacao.")
  if (any(cooks > 4 / n, na.rm = TRUE))
    result <- c(result, "Ha observacoes influentes; investigue-as antes de excluir qualquer dado.")
  if (!length(result)) "Nenhum alerta automatico foi identificado; confirme visualmente os graficos." else result
}
