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
  standardized <- stats::rstandard(fit)
  studentized <- stats::rstudent(fit)
  cooks <- stats::cooks.distance(fit)
  dffits <- stats::dffits(fit)
  dfbetas <- stats::dfbetas(fit)
  p <- length(stats::coef(fit))
  bf <- .brown_forsythe(data$.y, data$.x)
  bp <- .breusch_pagan(fit)
  gq <- .goldfeld_quandt(fit)
  cochran <- .cochran_test(data$.y, data$.x)
  grubbs <- .grubbs_test(residuals)
  lof <- .lack_of_fit(fit, data)
  quadratic <- .mandel_test(fit, data)
  normality <- .normality_tests(residuals)
  independence <- .independence_tests(fit)
  list(
    resumo_niveis = .level_summary(data),
    homocedasticidade = bf,
    breusch_pagan = bp,
    goldfeld_quandt = gq,
    cochran = cochran,
    grubbs = grubbs,
    falta_ajuste = lof,
    mandel = quadratic,
    influencia = data.frame(concentracao = data$.x, sinal = data$.y,
                             residuo = residuals, residuo_padronizado = standardized,
                             residuo_studentizado = studentized, alavancagem = h,
                             cook = cooks, dffits = dffits,
                             dfbeta_intercepto = dfbetas[, 1], dfbeta_inclinacao = dfbetas[, 2],
                             flag_cook = cooks > 4 / n, flag_alavancagem = h > 2 * p / n,
                             flag_dffits = abs(dffits) > 2 * sqrt(p / n),
                             flag_dfbeta = apply(abs(dfbetas), 1, max) > 2 / sqrt(n),
                             flag_studentizado = abs(studentized) > stats::qt(.975, stats::df.residual(fit))),
    normalidade = normality,
    independencia = independence,
    avisos = .diagnostic_messages(bf, bp, gq, cochran, grubbs, lof, quadratic, cooks, studentized, n)
  )
}

.normality_tests <- function(residuals) {
  n <- length(residuals)
  unavailable <- list(disponivel = FALSE, mensagem = "Amostra insuficiente para o teste.")
  shapiro <- if (n >= 3L && n <= 5000L) stats::shapiro.test(residuals) else unavailable
  ad <- if (n >= 8L) nortest::ad.test(residuals) else unavailable
  ks <- if (n >= 5L) nortest::lillie.test(residuals) else unavailable
  list(shapiro_wilk = shapiro, anderson_darling = ad, kolmogorov_smirnov_lilliefors = ks,
       ryan_joiner = .ryan_joiner_test(residuals))
}

.ryan_joiner_test <- function(x, alpha = .05) {
  n <- length(x)
  if (n < 4L) return(list(disponivel = FALSE, mensagem = "Ryan-Joiner requer ao menos quatro residuos."))
  p <- (seq_len(n) - 3 / 8) / (n + 1 / 4)
  estatistica <- stats::cor(sort(x), stats::qnorm(p))
  crit <- 1.0063 - .1288 / sqrt(n) - .6118 / n + 1.3505 / n^2
  list(disponivel = TRUE, estatistica = estatistica, valor_critico_5_percentual = crit,
       p_valor = .rj_p_value(estatistica, n), rejeita_normalidade = estatistica < crit,
       metodo = "Ryan-Joiner (aproximacao por correlacao normal)")
}

.rj_p_value <- function(r, n) {
  crit <- c(
    `0.10` = 1.0071 - .1371 / sqrt(n) - .3682 / n + .7780 / n^2,
    `0.05` = 1.0063 - .1288 / sqrt(n) - .6118 / n + 1.3505 / n^2,
    `0.01` = .9963 - .0211 / sqrt(n) - 1.4106 / n + 3.1791 / n^2
  )
  if (r > crit[1]) return("> 0.10")
  if (r < crit[3]) return("< 0.01")
  unname(stats::approx(x = rev(crit), y = c(.10, .05, .01), xout = r)$y)
}

.breusch_pagan <- function(fit) {
  z <- tryCatch(lmtest::bptest(fit), error = function(e) e)
  if (inherits(z, "error")) return(list(disponivel = FALSE, mensagem = conditionMessage(z)))
  list(disponivel = TRUE, estatistica = unname(z$statistic), gl = unname(z$parameter), p_valor = z$p.value)
}

.goldfeld_quandt <- function(fit) {
  z <- tryCatch(lmtest::gqtest(fit, fraction = .2), error = function(e) e)
  if (inherits(z, "error")) return(list(disponivel = FALSE, mensagem = conditionMessage(z)))
  list(disponivel = TRUE, estatistica = unname(z$statistic), gl = unname(z$parameter), p_valor = z$p.value)
}

.cochran_test <- function(y, x) {
  groups <- split(y, x); sizes <- lengths(groups)
  if (length(groups) < 2L || any(sizes < 2L) || length(unique(sizes)) != 1L)
    return(list(disponivel = FALSE, mensagem = "Cochran requer ao menos dois niveis com o mesmo numero de replicatas."))
  vars <- vapply(groups, stats::var, numeric(1)); q <- max(vars) / sum(vars)
  k <- length(groups); gl <- sizes[1] - 1L
  critical <- function(alpha) { f <- stats::qf(1 - alpha / k, gl, (k - 1) * gl); f / (f + k - 1) }
  p <- tryCatch(stats::uniroot(function(a) critical(a) - q, c(1e-8, .999999))$root, error = function(e) NA_real_)
  list(disponivel = TRUE, estatistica = q, gl = gl, p_valor_aproximado = p,
       valor_critico_5_percentual = critical(.05), variancias_por_nivel = vars)
}

.grubbs_test <- function(residuals) {
  if (length(residuals) < 3L) return(list(disponivel = FALSE, mensagem = "Grubbs requer ao menos tres residuos."))
  z <- outliers::grubbs.test(residuals, type = 10, opposite = FALSE)
  indice <- which.max(abs(residuals - mean(residuals)))
  list(disponivel = TRUE, estatistica = unname(z$statistic[[1]]), p_valor = unname(z$p.value),
       observacao = indice, residuo = residuals[indice], alternativa = z$alternative,
       nota = "Resultado para investigacao; nao exclua observacoes automaticamente.")
}

.independence_tests <- function(fit) {
  n <- length(stats::residuals(fit))
  dw <- tryCatch(lmtest::dwtest(fit), error = function(e) e)
  bg <- tryCatch(lmtest::bgtest(fit, order = min(2L, max(1L, floor(n / 3L)))), error = function(e) e)
  make <- function(z) if (inherits(z, "error")) list(disponivel = FALSE, mensagem = conditionMessage(z)) else
    list(disponivel = TRUE, estatistica = unname(z$statistic), gl = unname(z$parameter), p_valor = z$p.value, alternativa = z$alternative)
  list(durbin_watson = make(dw), breusch_godfrey = make(bg))
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

.diagnostic_messages <- function(bf, bp, gq, cochran, grubbs, lof, mandel, cooks, studentized, n) {
  result <- character()
  if (isTRUE(bf$disponivel) && is.finite(bf$p_valor) && bf$p_valor < .05)
    result <- c(result, "Ha evidencia de heterocedasticidade; compare OLS e WLS.")
  if (isTRUE(bp$disponivel) && bp$p_valor < .05) result <- c(result, "Breusch-Pagan indica variancia nao constante.")
  if (isTRUE(gq$disponivel) && gq$p_valor < .05) result <- c(result, "Goldfeld-Quandt indica variancia nao constante ao longo da concentracao.")
  if (isTRUE(cochran$disponivel) && is.finite(cochran$p_valor_aproximado) && cochran$p_valor_aproximado < .05) result <- c(result, "Cochran indica uma variancia de nivel desproporcionalmente alta.")
  if (isTRUE(grubbs$disponivel) && is.finite(grubbs$p_valor) && grubbs$p_valor < .05) result <- c(result, "Grubbs indica um residuo extremo; investigue a causa antes de qualquer exclusao.")
  if (isTRUE(lof$disponivel) && is.finite(lof$p_valor) && lof$p_valor < .05)
    result <- c(result, "Ha evidencia de falta de ajuste linear; investigue a faixa ou a curvatura.")
  if (isTRUE(mandel$disponivel) && is.finite(mandel$p_valor) && mandel$p_valor < .05)
    result <- c(result, "O termo quadratico melhora o ajuste; nao aprove a linearidade sem investigacao.")
  if (any(cooks > 4 / n, na.rm = TRUE))
    result <- c(result, "Ha observacoes influentes; investigue-as antes de excluir qualquer dado.")
  if (any(abs(studentized) > stats::qt(.975, max(1, n - 2)), na.rm = TRUE)) result <- c(result, "Ha residuos studentizados extremos; confirme a causa metrologica antes de qualquer decisao.")
  if (!length(result)) "Nenhum alerta automatico foi identificado; confirme visualmente os graficos." else result
}
