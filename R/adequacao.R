#' Evaluate calibration-curve adequacy against explicit criteria
#'
#' No universal acceptance limits are imposed. Supply only the limits applicable
#' to the analytical procedure, matrix, and intended range. A missing criterion
#' is reported as not assessed rather than silently approved.
#'
#' @param ajuste A `valcurva_fit` object.
#' @param criterios A named list. Supported names are `p_falta_ajuste_min`,
#'   `p_mandel_min`, `cv_max_percentual`, and
#'   `erro_relativo_retrocalculo_max_percentual`.
#' @return An object of class `valcurva_adequacao` with a criterion table,
#'   back-calculation table, and conclusion.
#' @export
avaliar_adequacao <- function(ajuste, criterios = list()) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  permitidos <- c("p_falta_ajuste_min", "p_mandel_min", "cv_max_percentual", "erro_relativo_retrocalculo_max_percentual")
  desconhecidos <- setdiff(names(criterios), permitidos)
  if (length(desconhecidos)) .valcurva_abort(paste("Criterios desconhecidos:", paste(desconhecidos, collapse = ", ")))
  if (length(criterios) && any(!vapply(criterios, function(x) is.numeric(x) && length(x) == 1L && is.finite(x), logical(1)))) .valcurva_abort("Cada criterio deve ser um numero finito.")
  d <- ajuste$dados; co <- stats::coef(ajuste$modelo)
  x_back <- (d$.y - co[1]) / co[2]
  back <- data.frame(concentracao = d$.x, sinal = d$.y, concentracao_retrocalculada = x_back, erro_relativo_percentual = 100 * (x_back - d$.x) / d$.x)
  resumo_back <- stats::aggregate(erro_relativo_percentual ~ concentracao, back, function(z) c(vies_percentual = mean(z), erro_absoluto_medio_percentual = mean(abs(z))))
  resumo_back <- data.frame(concentracao = resumo_back$concentracao, vies_percentual = resumo_back$erro_relativo_percentual[, "vies_percentual"], erro_absoluto_medio_percentual = resumo_back$erro_relativo_percentual[, "erro_absoluto_medio_percentual"])
  diag <- ajuste$diagnosticos; niveis <- diag$resumo_niveis
  make_row <- function(nome, valor, limite, sentido, disponivel = TRUE) {
    if (is.null(limite)) return(data.frame(criterio = nome, valor = valor, limite = NA_real_, status = "nao avaliado"))
    if (!disponivel || !is.finite(valor)) return(data.frame(criterio = nome, valor = NA_real_, limite = limite, status = "indisponivel"))
    passou <- if (sentido == "min") valor >= limite else valor <= limite
    data.frame(criterio = nome, valor = valor, limite = limite, status = if (passou) "atende" else "nao atende")
  }
  tab <- rbind(
    make_row("p de falta de ajuste", diag$falta_ajuste$p_valor, criterios$p_falta_ajuste_min, "min", diag$falta_ajuste$disponivel),
    make_row("p do teste de Mandel", diag$mandel$p_valor, criterios$p_mandel_min, "min", diag$mandel$disponivel),
    make_row("CV maximo por nivel (%)", max(niveis$cv_percentual, na.rm = TRUE), criterios$cv_max_percentual, "max"),
    make_row("erro relativo medio maximo por nivel (%)", max(resumo_back$erro_absoluto_medio_percentual), criterios$erro_relativo_retrocalculo_max_percentual, "max")
  )
  avaliados <- tab$status %in% c("atende", "nao atende")
  conclusao <- if (!any(avaliados)) "Sem decisao: informe ao menos um criterio de aceitacao aplicavel." else if (any(tab$status == "nao atende")) "Nao atende aos criterios informados." else if (any(tab$status == "indisponivel")) "Decisao incompleta: ao menos um criterio informado nao pode ser avaliado." else "Atende aos criterios informados; a aprovacao final requer revisao tecnica."
  structure(list(criterios = tab, retrocalculo_por_observacao = back, retrocalculo_por_nivel = resumo_back, conclusao = conclusao), class = "valcurva_adequacao")
}

#' Build an auditable calibration summary
#' @param ajuste A `valcurva_fit` object.
#' @param criterios Optional criteria forwarded to [avaliar_adequacao()].
#' @return A named list with model, diagnostics, adequacy assessment, and data.
#' @export
relatorio_auditoria <- function(ajuste, criterios = list()) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  list(modelo = list(metodo = ajuste$metodo, pesos = ajuste$estrategia_pesos, coeficientes = stats::coef(ajuste$modelo), covariancia_coeficientes = stats::vcov(ajuste$modelo), graus_de_liberdade = stats::df.residual(ajuste$modelo)), dados_por_nivel = ajuste$diagnosticos$resumo_niveis, diagnosticos = ajuste$diagnosticos, adequacao = avaliar_adequacao(ajuste, criterios))
}
